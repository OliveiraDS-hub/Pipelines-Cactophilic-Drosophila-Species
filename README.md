# Overview
This repository contains the code used to perform genomics and transcriptomics analysis from the manuscript "Oliveira, Daniel Siqueira de et al. Transposable elements as evolutionary driving force to ecological speciation in cactophilic Drosophila species. bioRxiv, p. 2024.03. 27.587021, 2024."

## Gene annotation

The four *D. mojavensis* subspecies and *D. arizonae* genomes were annotated based on the [*D. mojavensis* reference genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_018153725.1/).
The *D. buzzatii* and *D. koepferae* genomes were annotated with a de novo approach.

Dependencies:
  - [liftoff](https://github.com/agshumate/Liftoff)
  - [BRAKER3](https://github.com/Gaius-Augustus/BRAKER)
  - Download both reference genome and its gene annotation (GFF) from NCBI.

We highly recommend using BRAKER3 from a singularity image:
  - [singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)


To use the code available on this repository, make sure you downloaded the following files within the folder "Gene-annotation"
  - The reference genome (.fna) and its annotation (.GFF).
  - The genomes available on zenodo (link here)

**Command-line:**

For the genomes based on the reference genome:

```
bash liftoff_command.sh
```

To annotate the genomes of *D. buzzatii* and *D. koepferae*:

Download the fastq files available on Sup. Table X, using as preffix the file names from `rnaseq_sets_ids`

In the script `braker_command.sh`, set the variable `path_to_fastq_files` with the path where you have the fastq files.

```
bash braker_command.sh
```


## Differential expression analysis

Dependicies list:
  - [STAR](https://github.com/alexdobin/STAR)
  - [featureCounts](https://subread.sourceforge.net/featureCounts.html)
  - [deseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
  - [ggplot2](https://ggplot2.tidyverse.org)

## Permutation test

Permutation tests were perfomed to test enrichemnt and depletion of TEs in both upstream (-2kb) and downstream (+2kb) regions of HLAU genes. 
Dependencies list:
  - [regioneR](https://bioconductor.org/packages/release/bioc/html/regioneR.html)
  - [dplyr](https://dplyr.tidyverse.org)
  - [ggplot2](https://ggplot2.tidyverse.org)

BioCManager may be useful for installing regioneR:
- [BiocManager](https://bioconductor.org/install/)

**Command-line:**

#Example for 2kb upstream of Drosophila arizonae
```
cd Permutation-test/2kb_upstream/dari
```

#Unzip the bed files
```
unzip dari_bed.zip
```

Run the analysis
```
Rscript permut_and_plot.R dari_genes_2000upstream.bed dari_all_TEs.bed dari
```

This command line will generate one histogram per HLAU family, and a dotplot containing all the z-scores named as `dari_enrichment_permut.pdf`
