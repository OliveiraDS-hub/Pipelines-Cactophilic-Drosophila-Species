
## Differential expression analysis

Dependicies list:
  - [STAR](https://github.com/alexdobin/STAR)
  - [featureCounts](https://subread.sourceforge.net/featureCounts.html)
  - [deseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
  - [ggplot2](https://ggplot2.tidyverse.org)

This analysis is composed for three main steps: genome alignment, raw counts, and differential expression. The raw counts matrix produced from the two first steps are provided as csv files in `Differential-expression`


**Genome alignemt**:

In the script `Differential-expression/STAR_alignment.sh` set the preffix of the fastq libraries in the variable `libraries`. For example, if file names is `dmoj26_head_P1_R1.fastq.gz`, then replicate one preffix will be `dmoj26_head_P1`.
Set the variable `DATA` with the path to the folder containing the fasta and gff of the genomes. Finally, set the variable `OUTPUT_STAR` with the path to write the STAR output (.bam files).

Then, run:

```
bash STAR_alignment.sh
```

**Produce raw count matrix**:

In the script `Differential-expresison/featurecounts.sh`, set the variable `species_specific_gff` with gff file for the specific species being analyzed. In the same folder as `featurecounts.sh`, keep the bam files of the same species that you set the gff file.

Run the script for each species:

```
bash featurecounts.sh
```

**Differential expression with DEseq2**:

Four differential expresison analysis were performed, each one with their respective code in R:
  - mojavensis cluster, head tissue: mojavensis_head_DEA.R
  - mojavensis cluster, larvae tissue: mojavensis_larvae_DEA.R
  - buzzatii cluster, head tissue: buzzatii_head_DEA.R
  - buzzatii cluster, larvae tissue: buzzatii_larvae_DEA.R

Each code will produce the complete DEseq data for all genes, a table only with up-regulated genes, and another one with only down-regulated genes. We advise to run them on Rstudio.
