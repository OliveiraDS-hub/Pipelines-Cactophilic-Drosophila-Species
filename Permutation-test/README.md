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

