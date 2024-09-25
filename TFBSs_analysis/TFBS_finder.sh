#!/bin/bash

#$1 = TEs.bed
#$2 = 2kb upstream ORs
#$3 = genome .fasta

BED_TEs="/home/oliveirads/project_phd/2024_data/TE_annotation/REVIEW_july_2024"
BED_2kb_genes="/home/oliveirads/project_phd/2024_data/TFBS_prediction/REVIEW_july_2024/OR_genes"
GENOMES="/home/oliveirads/project_phd/2024_data/TFBS_prediction/REVIEW_july_2024/genomes"

species="dari
dmoj01
dmoj20
dmoj22
dmoj26
dbuz
dkoep"

TF="onecut_MA0235.1.meme
xbp1_MA2293.1.meme
acj6_MA2188.1.meme
fer1_MA2233.1.meme
zf30C_UN0798.1.meme"

while read -r species; do

  bedtools intersect -a "$BED_TEs"/"$species"_all_TEs_noverlap.bed -b "$BED_2kb_genes"/"$species"_or_2kbUP.bed -wa -wb > TEs_int_res.mbed
  cut -f1-6 TEs_int_res.mbed > TEs_upstream.bed
  freq=$(cut -f7-12 TEs_int_res.mbed | sort | uniq | wc -l | awk '{print $1}')
  echo "SPECIES = $species = $freq"

  if [[ "$species" == "dari" ]]; then
   	    bedtools getfasta -fi "$GENOMES"/Darizonae17.FlyeRacon3Medaka.fasta \
   	    -bed TEs_upstream.bed -name+ -s > TEs_upstream.fa
  elif [[ "$species" == "dmoj01" ]]; then
	      bedtools getfasta -fi "$GENOMES"/Dmojavensis01.FlyeRacon3Medaka.fasta \
        -bed TEs_upstream.bed -name+ -s > TEs_upstream.fa
  elif [[ "$species" == "dmoj20" ]]; then
        bedtools getfasta -fi "$GENOMES"/Dmojavensis20_j23.FlyeRacon3Medaka.fasta \
        -bed TEs_upstream.bed -name+ -s > TEs_upstream.fa
  elif [[ "$species" == "dmoj22" ]]; then
        bedtools getfasta -fi "$GENOMES"/Dmojavensis22.FlyeRacon3Medaka.fasta \
        -bed TEs_upstream.bed -name+ -s > TEs_upstream.fa
  elif [[ "$species" == "dmoj26" ]]; then
        bedtools getfasta -fi "$GENOMES"/Dmojavensis26.FlyeRacon3Medaka.fasta \
        -bed TEs_upstream.bed -name+ -s > TEs_upstream.fa
  elif [[ "$species" == "dkoep" ]]; then
        bedtools getfasta -fi "$GENOMES"/Dkoepferae.FlyeRacon3Medaka.fasta \
        -bed TEs_upstream.bed -name+ -s > TEs_upstream.fa
  elif [[ "$species" == "dbuz" ]]; then
        bedtools getfasta -fi "$GENOMES"/Dbuzatti.FlyeRacon3Medaka.fasta \
        -bed TEs_upstream.bed -name+ -s > TEs_upstream.fa
  fi


  while read -r specificTF; do
     if [ -d fimo_out ]; then
	rm -rf fimo_out
     fi
     fimo "$specificTF" TEs_upstream.fa  >> /dev/null 2>&1
     res=$(grep -v '# ' fimo_out/fimo.tsv)
     if [[ ! -z "$res" ]]; then
       TEcopies=$(grep -v '# ' fimo_out/fimo.tsv | grep -v 'motif_id' | cut -f3 | sed 's/::/\t/g; s/:/\t/g; s/(-)/\t/g; s/(+)/\t/g' | awk -v OFS="\t" '{print $1,$2,$3,$4}' | sort | uniq | sed '/^$/d')
       while read -r insertions; do
         if [[ ! -z $insertions ]]; then
           TE_fam=$(cut -f1 <<< "$insertions")
           TEstart=$(cut -f3 <<< "$insertions" | sed 's/-/\t/g' | cut -f1)
           TEend=$(cut -f3 <<< "$insertions" | sed 's/-/\t/g' | cut -f2)
           gene_ID=$(grep -w "$TE_fam" TEs_int_res.mbed | grep -w "$TEstart" | grep -w "$TEend" | cut -f10)
           pvalue=$(grep -w "$TE_fam" fimo_out/fimo.tsv | grep -w "$TEstart" | grep -w "$TEend" | cut -f8 | head -1)
           start_copy=$(grep -w "$TE_fam" fimo_out/fimo.tsv | grep -w "$TEstart" | grep -w "$TEend" | cut -f4 | head -1)
           end_copy=$(grep -w "$TE_fam" fimo_out/fimo.tsv | grep -w "$TEstart" | grep -w "$TEend" | cut -f5 | head -1)
           TFBS=$(grep -w "$TE_fam" fimo_out/fimo.tsv | grep -w "$TEstart" | grep -w "$TEend" | cut -f2 | head -1)
           echo -e "$TFBS\t$TE_fam\t$TEstart\t$TEend\t$start_copy\t$end_copy\t$pvalue\t$gene_ID\t${TE_fam}_${gene_ID}" >> "$species"_TFBSs.tsv
         fi
       done <<< "$TEcopies"
     fi
   done <<< "$TF"

done <<< "$species"
