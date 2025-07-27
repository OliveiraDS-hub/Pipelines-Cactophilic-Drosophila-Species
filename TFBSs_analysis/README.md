## Running TFBS identification in OR promoters

### Required softwares

- bedtools
- fimo (MEME suite)

All of them must be accessible directly from your PATH.
You also need to copy the six genomes to this folder. They can be downloaded from the [Zenodo repository](https://zenodo.org/records/13117512).
### How the pipeline works

All ortholog HLAU genes can be found in the folders location, acceptance and detox. Each species has a fasta file with the longest transcript of each gene.


### Running the tests

```
bash run_all_families.sh
```
