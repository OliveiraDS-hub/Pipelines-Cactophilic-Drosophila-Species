## Running TFBS identification in OR promoters

### Required softwares

- bedtools
- fimo (MEME suite)

All of them must be accessible directly from your PATH.

You also need to copy the six genomes to this folder. 
They can be downloaded from the [Zenodo repository](https://zenodo.org/records/13117512).

### How the pipeline works
All the following steps are perfomed automatically:
1. TEs overlapping OR promoters are identified in each species with bedtools intersect.
2. The sequences of these TE copies are extracted with bedtools getfasta.
3. TFBSs are identified with fimo, using the models:
   - acj6_MA2188.1.meme
   - fer1_MA2233.1.meme
   - onecut_MA0235.1.meme
   - xbp1_MA2293.1.meme
   - zf30C_UN0798.1.meme
4. The output will be a tab-delimited file with all detected TFBSs

### Running the tests

```
bash TFBS_finder_HLAU_genes.sh
bash TFBS_finder_nonHLAU_genes.sh
```
