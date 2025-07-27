import pandas as pd
pd.set_option('display.max_columns', None)
import subprocess
import argparse
import glob
import os
from Bio import SeqIO
from io import StringIO

parser = argparse.ArgumentParser()
parser.add_argument("fasta_file", help="Path to the genome FASTA file")
parser.add_argument("TE_database", help="Path to reference TEs in FASTA file")
parser.add_argument("species", help="Species name - the same as output folder")
parser.add_argument("threads", help="Total threads")
args = parser.parse_args()

fasta_file = args.fasta_file
TE_db = args.TE_database
species = args.species
threads = args.threads

### 80% of the custom library must be covered with 80% of identity
blast_result = subprocess.run("blastn -subject {} "
        "-query {} "
        "-outfmt 6 "
        "-qcov_hsp_perc 80 "
        "-perc_identity 80".format(TE_db, fasta_file),
        shell=True,
        capture_output=True,
        text=True)

# Check if the command was successful
if blast_result.returncode != 0:
    print("Error running blastn:", blast_result.stderr)
else:
    # Read the BLAST output into a Pandas DataFrame
    blast_columns = ["query", "subject", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    blast_out_df = pd.read_csv(StringIO(blast_result.stdout), sep="\t", names=blast_columns)

    df_top = blast_out_df.loc[blast_out_df.groupby('query')['bitscore'].idxmax()].reset_index(drop=True)
    df_top["consensus_new_id"] =  (
    df_top["subject"].str.split("#").str[0] + "_" +
    df_top["pident"].round(0).astype(int).astype(str) + "p#" +
    df_top["subject"].str.split("#").str[1])

    df_top.to_csv("top_matches_s5.tmp", sep = "\t", index = False)

    df_top['query'].to_csv("best_queries_s5.tmp", index=False, header=False)

    subprocess.run('python remove_seqs.py {} {} {}'.format(fasta_file, "classif_round1_done.fa_s5.tmp", "best_queries_s5.tmp"), shell=True)

df_top = df_top.rename(columns={"query": "consensus_id"})
df_rename = df_top[["consensus_id", "consensus_new_id"]].copy()
# df_rename.to_csv("rename1.tsv", index = False, header=None, sep ="\t")

### 80% of the database must be covered with 80% of identity
blast_result = subprocess.run("blastn -subject {} "
        "-query {} "
        "-outfmt 6 "
        "-qcov_hsp_perc 80 "
        "-perc_identity 80".format("classif_round1_done.fa_s5.tmp", TE_db),
        shell=True,
        capture_output=True,
        text=True)

# Check if the command was successful
if blast_result.returncode != 0:
    print("Error running blastn:", blast_result.stderr)
else:
    # Read the BLAST output into a Pandas DataFrame
    blast_columns = ["query", "subject", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    blast_out_df = pd.read_csv(StringIO(blast_result.stdout), sep="\t", names=blast_columns)

    df_top = blast_out_df.loc[blast_out_df.groupby('subject')['bitscore'].idxmax()].reset_index(drop=True)
    df_top["consensus_new_id"] =  (
    df_top["query"].str.split("#").str[0] + "_" +
    df_top["pident"].round(0).astype(int).astype(str) + "p#" +
    df_top["query"].str.split("#").str[1])
    
    df_top.to_csv("top_matches_80database_s5.tmp", sep = "\t", index = False)

    df_top['subject'].to_csv("best_queries_db_s5.tmp", index=False, header=False)

    subprocess.run('python remove_seqs.py {} {} {}'.format("classif_round1_done.fa_s5.tmp", "classif_round2_done.fa_s5.tmp", "best_queries_db_s5.tmp"), shell=True)

df_top = df_top.rename(columns={"subject": "consensus_id"})
df_rename = pd.concat([df_rename, df_top[["consensus_id", "consensus_new_id"]]], ignore_index=True)
# df_rename.to_csv("rename2.tsv", index = False, header=None, sep ="\t")


### Remove consensus with less than 200nt
## Create list of IDs with < 200nt
subprocess.run('samtools faidx classif_round2_done.fa_s5.tmp && awk \'$2 < 200\' classif_round2_done.fa_s5.tmp.fai | cut -f1 > cons_less_200nt_s5.tmp', shell=True)
## Remove sequences from fasta
subprocess.run('python remove_seqs.py {} {} {}'.format("classif_round2_done.fa_s5.tmp", "classif_round3_done.fa_s5.tmp", "cons_less_200nt_s5.tmp"), shell=True)
## Mask the remaining consensus with the database
print("Masking consensus... This process may take a while")
subprocess.run('RepeatMasker classif_round3_done.fa_s5.tmp -lib {} -pa {}'.format(TE_db, threads), shell = True) #, stdout=subprocess.DEVNULL)
subprocess.run('samtools faidx classif_round3_done.fa_s5.tmp', shell = True)
fai_df = pd.read_csv("classif_round3_done.fa_s5.tmp.fai", sep="\t", header=None, usecols=[0,1],
                            names=["consensus_id", "seq_length"])


##Define column names as per RepeatMasker's standard .out format
col_names = ["perc_div", "perc_del", "perc_ins",
    "query_sequence", "begin", "end", "left",
    "strand", "matching_repeat", "repeat_class_family",
    "r_begin", "r_end", "r_left", "ID"]

# Read the RepeatMasker .out file (skip the first 3 header lines)

cmd = f"cat classif_round3_done.fa_s5.tmp.out | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | tail -n +4 > rm_out_tab_delim.tsv_s5.tmp"
subprocess.run(cmd, shell=True)
rm_df = pd.read_csv("rm_out_tab_delim.tsv_s5.tmp", 
                 sep="\t", 
                 usecols=[4,5,6,9,10],
                 names=["consensus_id", "start", "end", "family", "class"])

# Ensure numeric types
rm_df["start"] = rm_df["start"].astype(int)
rm_df["end"] = rm_df["end"].astype(int)

# Compute interval length
rm_df["match_length"] = rm_df["end"] - rm_df["start"]

# Step 1: Sum match_length per consensus_id + family + class
family_sums = rm_df.groupby(["consensus_id", "family", "class"], as_index=False)["match_length"].sum()

# Step 2: Add total sequence length from fai_df
family_sums = family_sums.merge(fai_df, on="consensus_id", how="left")

# Step 3: Compute proportion per family
family_sums["family_prop"] = (family_sums["match_length"] / family_sums["seq_length"]) * 100

# Step 4: Compute total length per consensus_id + class
class_sums = family_sums.groupby(["consensus_id", "class"], as_index=False).agg({
    "match_length": "sum",
    "seq_length": "first"  # all rows per consensus_id have the same seq_length
})

# Step 5: Compute class-level proportion
class_sums["class_prop"] = (class_sums["match_length"] / class_sums["seq_length"]) * 100

# Step 6: Merge class proportions back to family_sums
final_df = family_sums.merge(class_sums[["consensus_id", "class", "class_prop"]],
                             on=["consensus_id", "class"], how="left")
final_df["class_prop"] = final_df["class_prop"].round(4)
final_df["family_prop"] = final_df["family_prop"].round(4)

final_df.to_csv("prop_coverage.tsv_s5.tmp", sep = "\t", index = False)


### Select consensus with > 80 with specific family
high_family_prop_df = final_df[final_df["family_prop"] > 80]
highest_family = high_family_prop_df.loc[
    high_family_prop_df.groupby(["consensus_id"])["family_prop"].idxmax()
].reset_index(drop=True)
highest_family["family_prop"] = highest_family["family_prop"].clip(upper=100)
highest_family["consensus_new_id"] = highest_family["family"] + "_" + highest_family["family_prop"].astype(int).astype(str) + "p#" + highest_family["class"] 
highest_family.to_csv("prop_family_80.tsv_s5.tmp", sep = "\t", index = False)
print(highest_family)

df_rename = pd.concat([df_rename, highest_family[["consensus_id", "consensus_new_id"]]], ignore_index=True)

## Subset consensus to check class
high_class_prop_df = final_df[~final_df["consensus_id"].isin(highest_family["consensus_id"])]
high_class_prop_df = high_class_prop_df[high_class_prop_df["class_prop"] > 50]
highest_class = high_class_prop_df.loc[
    high_class_prop_df.groupby(["consensus_id"])["class_prop"].idxmax()
].reset_index(drop=True)
highest_class["class_prop"] = highest_class["class_prop"].clip(upper=100)
highest_class.to_csv("prop_class_80.tsv_s5.tmp", sep = "\t", index = False)
highest_class["consensus_new_id"] = species + "_DSO_" + highest_class["class_prop"].round(4).astype(str) + "p#"  +highest_class["class"] 

df_rename = pd.concat([df_rename, highest_class[["consensus_id", "consensus_new_id"]]], ignore_index=True)

# Remove rows where 'consensus_id' contains any of the unwanted keywords (case-insensitive)
filtered_df = df_rename[~df_rename["consensus_id"].str.contains("satellite|complexity|rrna", case=False, na=False)]
# Save the filtered list
filtered_df["consensus_id"].to_csv("consensus_list_s5.tmp", index=False, header=False)

subprocess.run('python remove_seqs.py {} {} {}'.format("classif_round3_done.fa_s5.tmp", "unclassified_consensus_s5.fa", "consensus_list_s5.tmp"), shell=True)


from Bio import SeqIO
# Load IDs from file
with open("consensus_list_s5.tmp") as f:
    id_set = set(line.strip() for line in f if line.strip())

# Filter and write
with open("polished_TEs_s5.tmp", "w") as out_f:
    records = (record for record in SeqIO.parse(fasta_file, "fasta") if record.id in id_set)
    SeqIO.write(records, out_f, "fasta")

import pandas as pd

# Load the dataframe
id_map = dict(zip(df_rename["consensus_id"], df_rename["consensus_new_id"]))

# Read and rename sequences
with open("polished_TEs_s5.tmp") as infile, open(f"{species}/polished_TEs_s5.fa", "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        if record.id in id_map:
            record.id = id_map[record.id]
            record.name = ""
            record.description = ""
        SeqIO.write(record, outfile, "fasta")

# import shutil
# os.makedirs(species, exist_ok=True)

# files = glob.glob("polished*.fa")
# for f in files:
#     shutil.move(f, species)
# if os.path.exists(rm_out[0]):
#     shutil.move(rm_out[0], species)
# if os.path.exists("unclassified_consensus_s5.fa"):
#     shutil.move("unclassified_consensus_s5.fa", species)