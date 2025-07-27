#!/bin/bash

BED_TEs="bed_TEs"
BED_2kb_genes="nonHLAU_2kbUP"

species="D_arizonae
D_buzzatii
D_koepferae
D_moj_mojavensis
D_moj_sonorensis
D_moj_wrigleyi"

TF="onecut_MA0235.1.meme
xbp1_MA2293.1.meme
acj6_MA2188.1.meme
fer1_MA2233.1.meme
zf30C_UN0798.1.meme"

gene_fam="serine"

while read -r gene_fam; do
    while read -r species; do
        if [[ ! -z "$species" ]]; then
        bedtools intersect -a "$BED_TEs"/"$species"_TEs.bed -b "$BED_2kb_genes"/"$species"_OR_2kbUP.bed -wa -wb > TEs_int_res_$gene_fam.mbed
        cut -f1-6 TEs_int_res_"$gene_fam".mbed > TEs_upstream_"$gene_fam".bed
        freq=$(cut -f7-12 TEs_int_res_"$gene_fam".mbed | sort | uniq | wc -l | awk '{print $1}')
        echo "SPECIES = $species = $freq"

        if [[ "$species" == "D_arizonae" ]]; then
                bedtools getfasta -fi D_arizonae_genome.fasta \
                -bed TEs_upstream_"$gene_fam".bed -name+ -s > TEs_upstream_"$gene_fam".fa
        elif [[ "$species" == "D_moj_mojavensis" ]]; then
                bedtools getfasta -fi D_moj_mojavensis_genome.fasta \
                -bed TEs_upstream_"$gene_fam".bed -name+ -s > TEs_upstream_"$gene_fam".fa
        elif [[ "$species" == "D_moj_wrigleyi" ]]; then
                bedtools getfasta -fi D_moj_wrigleyi_genome.fasta \
                -bed TEs_upstream_"$gene_fam".bed -name+ -s > TEs_upstream_"$gene_fam".fa
        elif [[ "$species" == "D_moj_sonorensis" ]]; then
                bedtools getfasta -fi D_moj_sonorensis_genome.fasta \
                -bed TEs_upstream_"$gene_fam".bed -name+ -s > TEs_upstream_"$gene_fam".fa
        elif [[ "$species" == "D_koepferae" ]]; then
                bedtools getfasta -fi D_koepferae_genome.fasta \
                -bed TEs_upstream_"$gene_fam".bed -name+ -s > TEs_upstream_"$gene_fam".fa
        elif [[ "$species" == "D_buzzatii" ]]; then
                bedtools getfasta -fi D_buzzatii_genome.fasta \
                -bed TEs_upstream_"$gene_fam".bed -name+ -s > TEs_upstream_"$gene_fam".fa
        fi


        while read -r specificTF; do
            if [ -d fimo_out ]; then
                rm -rf fimo_out
            fi
            fimo "$specificTF" TEs_upstream_"$gene_fam".fa  #>> /dev/null 2>&1
            mv fimo_out "$species"_"$specificTF"_"$gene_fam"_fimo_out
            res=$(grep -v '# ' "$species"_"$specificTF"_"$gene_fam"_fimo_out/fimo.tsv)
            
            if [[ ! -z "$res" ]]; then
                TEcopies=$(grep -v '# ' "$species"_"$specificTF"_"$gene_fam"_fimo_out/fimo.tsv | grep -v 'motif_id' | cut -f3 | sed 's/::/\t/g; s/:/\t/g; s/(-)/\t/g; s/(+)/\t/g' | awk -v OFS="\t" '{print $1,$2,$3,$4}' | sort | uniq | sed '/^$/d')
            
            while read -r insertions; do
                if [[ ! -z $insertions ]]; then
                    TE_fam=$(cut -f1 <<< "$insertions")
                    TEstart=$(cut -f3 <<< "$insertions" | sed 's/-/\t/g' | cut -f1)
                    TEend=$(cut -f3 <<< "$insertions" | sed 's/-/\t/g' | cut -f2)
                    gene_ID=$(grep -w "$TE_fam" TEs_int_res_"$gene_fam".mbed | grep -w "$TEstart" | grep -w "$TEend" | cut -f10)
                    score=$(grep -w "$TE_fam" "$species"_"$specificTF"_"$gene_fam"_fimo_out/fimo.tsv | grep -w "$TEstart" | grep -w "$TEend" | cut -f7 | head -1)
                    pvalue=$(grep -w "$TE_fam" "$species"_"$specificTF"_"$gene_fam"_fimo_out/fimo.tsv | grep -w "$TEstart" | grep -w "$TEend" | cut -f8 | head -1)
                    qvalue=$(grep -w "$TE_fam" "$species"_"$specificTF"_"$gene_fam"_fimo_out/fimo.tsv | grep -w "$TEstart" | grep -w "$TEend" | cut -f9 | head -1)
                    start_copy=$(grep -w "$TE_fam" "$species"_"$specificTF"_"$gene_fam"_fimo_out/fimo.tsv | grep -w "$TEstart" | grep -w "$TEend" | cut -f4 | head -1)
                    end_copy=$(grep -w "$TE_fam" "$species"_"$specificTF"_"$gene_fam"_fimo_out/fimo.tsv | grep -w "$TEstart" | grep -w "$TEend" | cut -f5 | head -1)
                    TFBS=$(grep -w "$TE_fam" "$species"_"$specificTF"_"$gene_fam"_fimo_out/fimo.tsv | grep -w "$TEstart" | grep -w "$TEend" | cut -f2 | head -1)
                    echo -e "$TFBS\t$TE_fam\t$TEstart\t$TEend\t$start_copy\t$end_copy\t$score\t$pvalue\t$qvalue\t$gene_ID\t${TE_fam}_${gene_ID}" >> "$species"_"$gene_fam"_TFBSs.tsv
                fi
            done <<< "$TEcopies"
            fi
        done <<< "$TF"
        fi
    done <<< "$species"

# rm TEs_int_res_"$gene_fam".mbed TEs_upstream_"$gene_fam".bed TEs_upstream_"$gene_fam".fa
done <<< "$gene_fam"