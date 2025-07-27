import subprocess
import os
import pandas as pd
import argparse
import pysam
from collections import defaultdict
import glob
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("genome", help="Path to the genome FASTA file")
parser.add_argument("gene_annotation", help="Path to GFF/GTF annotation file")
parser.add_argument("threads", help="Number of threads")
parser.add_argument("species", help="Species name - the same as output folder")
parser.add_argument("mate1", help="Path to the mate 1 FASTQ")
parser.add_argument("mate2", help="Path to the mate 2 FASTQ")
parser.add_argument("strandness", help="Strandness of the RNA-seq data - either fwd-stranded or rf-stranded")

args = parser.parse_args()

genome = args.genome
gene_annotation = args.gene_annotation
threads = args.threads
species = args.species
mate1 = args.mate1
mate2 = args.mate2
strandness = args.strandness

def stranded_reads(strand_option, bam_f, species, threads):
    if strand_option == "fwd-stranded":
        pair_sense = ["rev", "fwd"]
    elif strand_option == "rf-stranded":
        pair_sense = ["fwd", "rev"]

    for index, sense in enumerate(pair_sense):
        print(f"Creating bed file for {sense} reads...")
        if index == 0:
            ### First strand
            subprocess.call(['samtools', 'view', '-@', str(threads), '-b', '-f', '128', '-F', '16', f"{species}/{bam_f}", '-o', f"{species}/{sense}1_f.bam"]) #, stderr=subprocess.DEVNULL)

            subprocess.call(['samtools', 'view', '-@', str(threads), '-b', '-f', '80', f"{species}/{bam_f}", '-o', f"{species}/{sense}2_f.bam"], stderr=subprocess.DEVNULL)
        elif index == 1:
            ### Second strand
            subprocess.call(['samtools', 'view', '-@', str(threads), '-b', '-f', '144', f"{species}/{bam_f}", '-o', f"{species}/{sense}1_f.bam"], stderr=subprocess.DEVNULL)

            subprocess.call(['samtools', 'view', '-@', str(threads), '-b', '-f', '64', '-F', '16', f"{species}/{bam_f}", '-o', f"{species}/{sense}2_f.bam"], stderr=subprocess.DEVNULL)
        
        ### Combine alignments that originate on the forward strand.
        subprocess.call(['samtools', 'merge', '-@', str(threads), '-f', f"{species}/{sense}.bam", f"{species}/{sense}1_f.bam", f"{species}/{sense}2_f.bam"])

        subprocess.run('bedtools bamtobed -split -i {0}/{1}.bam > {0}/{1}.bed'.format(species, sense), shell = True)

def filter_uniquely_mapped_reads(input_bam, output_bam, threads=4):
    """
    Keeps only reads (single or paired) that align to a unique position in the genome.
    Keeps both mates if both are mapped uniquely.
    Keeps a single mate if it's uniquely mapped and its pair is unmapped.
    Removes reads that have multiple alignments (multi-mappers).
    Requires BAM sorted by read name.
    """
    read_counts = defaultdict(int)

    # First pass: count how many alignments each read name has (including both mates)
    with pysam.AlignmentFile(input_bam, "rb") as bam:
        for read in bam:
            if not read.is_unmapped:
                read_counts[read.query_name] += 1

    # Second pass: write reads where total count for that read name == 1 or 2
    with pysam.AlignmentFile(input_bam, "rb") as bam, \
         pysam.AlignmentFile(output_bam, "wb", template=bam, threads=threads) as outbam:
        for read in bam:
            if read_counts[read.query_name] in (1, 2):
                outbam.write(read)

def get_nonoverlapping_TEs(all_bed, overlap_bed, output_bed):
    cols = ['Chrom', 'Start', 'End', 'Family', 'Score', 'Strand']

    # Load all TEs and overlapped ones
    all_df = pd.read_csv(all_bed, sep='\t', names=cols)
    overlap_df = pd.read_csv(overlap_bed, sep='\t', names=cols)

    # Drop duplicates if any
    overlap_df = overlap_df.drop_duplicates()

    # Use merge with indicator to find rows not in overlap
    merged = all_df.merge(overlap_df, how='left', on=cols, indicator=True)
    non_overlapping = merged[merged['_merge'] == 'left_only'].drop(columns=['_merge'])

    # Write to BED
    non_overlapping.to_csv(output_bed, sep='\t', header=False, index=False)

def get_top_insertions_by_family(bed_file, bam_file, output_bed):
    # Read BED file into DataFrame
    cols = ['chrom', 'start', 'end', 'family', 'score', 'strand']
    df = pd.read_csv(bed_file, sep='\t', names=cols)

    # List to hold results
    results = []

    # Group by TE family
    for family, group in df.groupby('family'):
        family_rows = []

        for _, row in group.iterrows():
            region = f"{row['chrom']}:{row['start']+1}-{row['end']}"  # +1 for 1-based BED to SAM
            cmd = f"samtools view -c {bam_file} {region}"
            try:
                read_count = int(subprocess.check_output(cmd, shell=True).decode().strip())
                row_with_count = row.copy()
                row_with_count['reads'] = read_count
                family_rows.append(row_with_count)
            except subprocess.CalledProcessError as e:
                print(f"Error running samtools for {region}: {e}")
                continue

        # Get top 5 insertions with highest read count
        top5 = sorted(family_rows, key=lambda x: x['reads'], reverse=True)[:5]
        results.extend(top5)

    # Create DataFrame from results and update the score column with read counts
    final_df = pd.DataFrame(results)
    final_df['score'] = final_df['reads']
    final_df = final_df.drop(columns='reads')

    # Save to output bed file
    final_df.to_csv(output_bed, sep='\t', header=False, index=False)

def add_fwd_rev_counts(input_bed, fwd_bam, rev_bam, output_bed):

    # Read the top5 TE insertions BED
    cols = ['chrom', 'start', 'end', 'family', 'score', 'strand']
    df = pd.read_csv(input_bed, sep='\t', names=cols)

    # Initialize lists to store read counts
    fwd_counts = []
    rev_counts = []

    # Loop through each insertion
    for _, row in df.iterrows():
        region = f"{row['chrom']}:{row['start']+1}-{row['end']}"  # BED is 0-based, SAM is 1-based

        # Count in fwd.bam
        try:
            fwd = int(subprocess.check_output(f"samtools view -c {fwd_bam} {region}", shell=True).decode().strip())
        except subprocess.CalledProcessError as e:
            print(f"Error with fwd.bam for {region}: {e}")
            fwd = 0

        # Count in rev.bam
        try:
            rev = int(subprocess.check_output(f"samtools view -c {rev_bam} {region}", shell=True).decode().strip())
        except subprocess.CalledProcessError as e:
            print(f"Error with rev.bam for {region}: {e}")
            rev = 0

        fwd_counts.append(fwd)
        rev_counts.append(rev)

    # Add to dataframe
    df['fwd'] = fwd_counts
    df['rev'] = rev_counts

    # Save to output
    df.to_csv(output_bed, sep='\t', header=False, index=False)

def compute_strand_agreement(df):
    """
    Takes a DataFrame with columns:
    [chrom, start, end, family, total_reads, strand, fwd_reads, rev_reads]

    Returns a DataFrame grouped by family with:
    [family, correct_reads, incorrect_reads, prop_correct, prop_incorrect]
    """
    # Compute correct and incorrect assignments per row
    def strand_match(row):
        if row['strand'] == '+':
            return pd.Series({'correct': row['fwd_reads'], 'incorrect': row['rev_reads']})
        else:
            return pd.Series({'correct': row['rev_reads'], 'incorrect': row['fwd_reads']})

    strand_counts = df.apply(strand_match, axis=1)
    df[['correct', 'incorrect']] = strand_counts

    # Group by family and sum correct/incorrect
    result = df.groupby('family')[['correct', 'incorrect']].sum().reset_index()

    # Calculate proportions
    total = result['correct'] + result['incorrect']
    result['prop_correct'] = (result['correct'] / total).astype(float)
    result['prop_incorrect'] = (result['incorrect'] / total).astype(float)

    return result

def extract_family_sequences(fasta_file, output, family_list1, family_list2):
    """
    Extract sequences from fasta_file and write to output1/output2 based on family lists.
    Applies reverse complement if correct_frame is False.
    """
    set1 = set(family_list1)
    set2 = set(family_list2)

    with open(output, "w") as out_f:
        for record in SeqIO.parse(fasta_file, "fasta"):
            fasta_id_prefix = record.id.split('#')[0]
            if fasta_id_prefix in set1:
                rec = record[:]  # Make a copy to avoid side effects
                rec.seq = rec.seq.reverse_complement()
                rec.description += " [validated_strand]"
                SeqIO.write(rec, out_f, "fasta")
            elif fasta_id_prefix in set2:
                rec = record[:]
                rec.description += " [validated_strand]"
                SeqIO.write(rec, out_f, "fasta")
            else:
                rec = record[:]
                SeqIO.write(rec, out_f, "fasta")


### Create bowtie2 index
os.makedirs(f'{species}/index', exist_ok=True)
print("Creating bowtie2 index")
subprocess.run(f'bowtie2-build {genome} {species}/index/genome_index'.format(genome), shell=True, stdout=subprocess.DEVNULL)

### bowtie2-alignment
print("Performing bowtie2 alignment")
subprocess.run(f'bowtie2 -x {species}/index/genome_index -1 {mate1} -2 {mate2} -S {species}/{species}_aln.sam --very-sensitive -k 2 --threads {threads}', shell=True)
## SAM to BAM
print("Converting SAM to BAM")
subprocess.run(f'samtools view -bS {species}/{species}_aln.sam | samtools view -h -F 4 -o {species}/{species}_aln.bam -@ {threads}', shell=True)
### Sort unique BAM
print("Sorting BAM")
subprocess.run(f'samtools sort -n -@ {threads} -o {species}/{species}_aln_st.bam {species}/{species}_aln.bam ', shell=True)

### Filter uniquely mapped reads
filter_uniquely_mapped_reads(f'{species}/{species}_aln_st.bam', f'{species}/{species}_unique.bam', threads=8)

## Extract stranded 
stranded_reads(strandness, f'{species}_unique.bam', species, threads)
std_files = glob.glob(f'{species}/*_f.bam')
for f in std_files:
    os.remove(f)


### Subset "free-living" TEs
if not os.path.exists(f'{genome}.out'):
    subprocess.run(f'RepeatMasker {genome} -lib {species}/polished_TEs_s5.fa -cutoff 250 -norna -gff -a -s -pa {threads}', shell=True, stdout=subprocess.DEVNULL)

cmd = (f"egrep -v 'Satellite|Simple_repeat|rRNA|Low_complexity|RNA|ARTEFACT' {genome}.out | "
    "tr -s ' ' | "
    "sed 's/^ *//g' | "
    "tr ' ' '\t' | "
    "tail -n +4 | "
    "awk '{Sense=$9;sub(/C/,\"-\",Sense);$9=Sense ;print $5\"\t\"$6\"\t\"$7\"\t\"$10\"\t\"$1\"\t\"$9}' | "
    f"sed 's/ /\\t/g'  > {species}/rm_TEs.bed")
subprocess.run(cmd, shell=True)


### gff to bed with gene regions
cmd = (f"awk '$3 == \"gene\"' {gene_annotation} | "
       "awk -v OFS=\"\t\" '{print $1,$4,$5,$9}' | "
       f"sed 's/ID=gene-//; s/;.*//' > {species}/gene.bed")
subprocess.run(cmd, shell=True)

### Intersect genes with TEs
subprocess.run(f'bedtools intersect -a {species}/rm_TEs.bed -b {species}/gene.bed -wa -wb -v > {species}/rm_TEs_no-exons.bed'.format(species), shell = True)

### TEs not overlapping other TEs
subprocess.run(f"bedtools intersect -a {species}/rm_TEs_no-exons.bed -b {species}/rm_TEs_no-exons.bed -wa -wb -f 0.1 | awk '!(($2 == $8) && ($3 == $9) && ($4 == $10))' | cut -f1-6 | sort | uniq > {species}/TEs_overlap_each_other.bed".format(species), shell=True)
get_nonoverlapping_TEs(all_bed=f"{species}/rm_TEs_no-exons.bed",
    overlap_bed=f"{species}/TEs_overlap_each_other.bed",
    output_bed=f"{species}/TEs_nonoverlapping.bed")

### Sort unique BAM
subprocess.run(f'bedtools bamtobed -i {species}/{species}_unique.bam > {species}/{species}_unique.bed'.format(species), shell=True)
subprocess.run(f'samtools sort {species}/{species}_unique.bam -@ {threads} -o {species}/{species}_unique_st.bam', shell=True)
subprocess.run(f'samtools index {species}/{species}_unique_st.bam', shell=True)

### Get the top5 insertions with more reads
get_top_insertions_by_family(
    bed_file=f'{species}/TEs_nonoverlapping.bed',
    bam_file=f'{species}/{species}_unique_st.bam',
    output_bed=f'{species}/top5_TE_insertions.bed')




### Index stranded bam files
subprocess.run(f'samtools sort -@ {threads} {species}/rev.bam -o {species}/rev_st.bam', shell=True)
subprocess.run(f'samtools index {species}/rev_st.bam', shell=True)
subprocess.run(f'samtools sort -@ {threads} {species}/fwd.bam -o {species}/fwd_st.bam', shell=True)
subprocess.run(f'samtools index {species}/fwd_st.bam', shell=True)

### Count reads per strand
add_fwd_rev_counts(input_bed=f'{species}/top5_TE_insertions.bed',
    fwd_bam=f'{species}/fwd_st.bam',
    rev_bam=f'{species}/rev_st.bam',
    output_bed=f'{species}/top5_TE_insertions_stranded.bed')


### Calculate correct vs incorrect mapping rates
cols = ['chrom', 'start', 'end', 'family', 'total_reads', 'strand', 'fwd_reads', 'rev_reads']
df = pd.read_csv(f"{species}/top5_TE_insertions_stranded.bed", sep='\t', names=cols)
strand_stats = compute_strand_agreement(df)
strand_stats.to_csv(f'{species}/strandess_res.tsv', sep = "\t", index = False)

### Select families with incorrect strand
incorrect_strand = strand_stats[strand_stats["prop_incorrect"] >= 0.75]
incorrect_strand = incorrect_strand[(incorrect_strand["correct"] + incorrect_strand["incorrect"]) > 100]

### Select families with correct strand
correct_strand = strand_stats[strand_stats["prop_correct"] >= 0.75]
correct_strand = correct_strand[(correct_strand["correct"] + correct_strand["incorrect"]) > 100]

                

### Create fasta with corrected consensus
extract_family_sequences(fasta_file=f"{species}/polished_TEs_s5.fa",
    output=f"{species}/polished_TEs_s6.fa",
    family_list1=incorrect_strand["family"].tolist(),
    family_list2=correct_strand["family"].tolist())

total_corrected_strand = len(correct_strand) + len(incorrect_strand)
print(f"{total_corrected_strand} consensus had strandness validated by RNA-seq")