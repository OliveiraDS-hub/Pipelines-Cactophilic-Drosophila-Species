#!/bin/bash

path_to_fastq_files=""

singularity exec -B $PWD:$PWD braker3.sif braker.pl \
--genome=D_koepferae_genome.fasta \
--rnaseq_sets_ids=dkoep_head,dkoep_larv,dkoep_male,dkoep_ovary,dkoep_testes \
--rnaseq_sets_dirs=$path_to_fastq_files \
--prot_seq=dmoj_NCBI_PT.fa \
--threads 30 > braker_dkoep.log 2>&1

mv braker braker_dkoep_FINAL

singularity exec -B $PWD:$PWD braker3.sif braker.pl \
--genome=D_buzzatii_genome.fasta \
--rnaseq_sets_ids=dbuz_egg_rep1,dbuz_head_P1,dbuz_larv_P1,dbuz_pupae_rep1,dbuz_female_rep1,dbuz_male_rep1 \
--rnaseq_sets_dirs=$path_to_fastq_files \
--prot_seq=dmoj_NCBI_PT.fa \
--threads 30 > braker_dbuz.log 2>&1

mv braker braker_dbuz_FINAL
