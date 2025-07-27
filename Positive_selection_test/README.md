## Running selection test with CODEML

### Required softwares

- seqtk
- embo transeq
- mafft
- pal2nal.pl
- codeml

All of them must be accessible directly from your PATH.

### How the pipeline works

All ortholog HLAU genes can be found in the folders location, acceptance and detox. Each species has a fasta file with the longest transcript of each gene.

For each ortholog, the pipeline will automatically perform the following steps:
1. extract the sequence from each species with seqtk
2. convert cds to protein with embo transeq
3. align them with mafft
4. correct to codon-based alignment with pal2nal
5. perform signatures of selection with codeml (branch-site model)
6. perform stastitical analysis in R

The output is a tab-delimited file with the gene ID and the corrected p-value for positive selection. 


### Running the tests

```
bash run_all_families.sh
```
