#!/bin/bash

### Example variables for dmoj26 (D. moj. sonorensis)

declare -a libraries=("dmoj26_head_P1"
"dmoj26_head_P2"
"dmoj26_head_P3")

DATA="" ### path to the genome and gff files of D. moj. sonorensis
OUTPUT_STAR="" ### select the path to store the output (.bam files)

mkdir "$OUTPUT_STAR"/dmoj26_star_index "$OUTPUT_STAR"/dmoj26_star_align "$OUTPUT_STAR"/dmoj26_counts "$OUTPUT_STAR"/dmoj26_star_align/stat

STAR --runThreadN 16 --genomeSAindexNbases 12 --runMode genomeGenerate --genomeDir "$OUTPUT_STAR"/dmoj26_star_index --genomeFastaFiles "$DATA"/D_moj_sonorensis_genome.fasta --sjdbGTFfile "$DATA"/D_moj_sonorensis_genes.gff --sjdbOverhang 99

for file in ${libraries[@]}; do
        echo "Alignment "$file""
        STAR --genomeDir "$OUTPUT_STAR"/dmoj26_star_index --runThreadN 16 --readFilesCommand zcat --readFilesIn "$file"_R1.fastq.gz "$file"_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$OUTPUT_STAR"/star_align/"$file"_
        mv "$OUTPUT_STAR"/star_align/*final.out "$OUTPUT_STAR"/star_align/stat
        rm "$OUTPUT_STAR"/star_align/*.out
done
