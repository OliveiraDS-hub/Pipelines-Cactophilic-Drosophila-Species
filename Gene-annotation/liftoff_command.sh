#!/bin/bash

ref_genome="GCF_018153725.1_ASM1815372v1_genomic.fna"
ref_gtf="GCF_018153725.1_ASM1815372v1_genomic.gff"

for genome in *.fasta; do
  out=$(basename "$genome" | sed 's/.fasta//g')
  liftoff "$genome" "$ref_genome" \
  -g "$ref_gtf" \
  -cds -o "$out"_REFncbi_GCF_018153725.gff \
  -p 16 -exclude_partial -a 0.5 -s 0.75
done
