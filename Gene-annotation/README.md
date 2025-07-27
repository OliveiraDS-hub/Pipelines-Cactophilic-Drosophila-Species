## Running gene annotation

### Required softwares
- liftoff
- singularity
- braker3 singularity container (named braker3.sif)

All of them must be accessible directly from your PATH.

You also need to copy the six genomes to this folder. They can be downloaded from the [Zenodo repository](https://zenodo.org/records/13117512).


### How the pipeline works
Gene annotation with lifftoff will use the publicly available *D. moj. wrigleyi* genome (.fasta) and gene annotation (.gff) from NCBI, version GCF_018153725.1. 

### Running gene annotation

```
bash liftoff_command.sh
```

To run BRAKER, you first need to download the RNA-seq data. You can have the list for *D. buzzatii* and *D. koepferae* in the Supplementary table 11 of the manuscript.
The fastq files **must** have these prefixes:
- dkoep_head, dkoep_larv, dkoep_male, dkoep_ovary, dkoep_testes
- dbuz_egg_rep1, dbuz_head_P1, dbuz_larv_P1, dbuz_pupae_rep1, dbuz_female_rep1, dbuz_male_rep1

In the `braker_command.sh` script, add the path to the folder with all fastq files in the variable `path_to_fastq_files`

Then, run:

```
bash braker_command.sh
```

Each script may take over 24h to finish.
