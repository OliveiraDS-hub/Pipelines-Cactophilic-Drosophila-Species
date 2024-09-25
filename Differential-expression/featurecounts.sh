#!/bin/bash

species_specific_gff="" ### specity the path and file name of the gff file

for bam_file in *bam; do
  out_name=$(basename "$bam_file" | sed 's/.bam//g')
  featureCounts -a $species_specific_gff \
  "$bam_file" \
  -o "$out_name".tsv \
  -p -B  -C -T 20 -t exon -g gene
  cut -f1,7,8 "$out_name".tsv > tmp.tmp && mv tmp.tmp "$out_name".tsv
done
