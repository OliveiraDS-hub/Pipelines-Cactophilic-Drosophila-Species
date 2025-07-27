#!/bin/bash
set -e
#$1 folder with CDSs fasta
#$2 HS family name for output

if [[ ! -d output_"$2" ]]; then
	mkdir output_"$2"
fi

cat "$1"/*fa | grep '>' | sed 's/_.*//g; s/>//' | sort | uniq -c | awk '$1 == 6' | awk '{print $2}' > common_geneIDs.lst

#genes () {
while read -r line; do
	for fasta_CDS in "$1"/*.fa; do
		species_name=$(basename "$fasta_CDS" | sed 's/_.*//g')
		grep "$line" "$fasta_CDS" | sed 's/>//g' > spec_ID.txt
		seqtk subseq "$fasta_CDS" spec_ID.txt > "$species_name"_spec_gene.fa
	done

	cat *_spec_gene.fa | sed 's/>.*_/>/g' > merged_genes.fa
	transeq -sequence merged_genes.fa -frame 1 -trim Y -outseq mafft_data_PT.fa
	sed -i 's/_1//g' mafft_data_PT.fa
	mafft --thread 16 mafft_data_PT.fa > mafft_PT.aln #> /dev/null 2>&1

	out=$(sed 's/>//g; s/.fa//g' <<< "$line")
	pal2nal.pl mafft_PT.aln merged_genes.fa -nogap -output paml > p2n_"$out".fa #1> /dev/null
	rm spec_ID.txt *_spec_gene.fa

	if [[ -s p2n_"$out".fa ]]; then
		sed -i "s/seqfile =.*/seqfile = p2n_"$out".fa/" control_file_null.ctrl
                sed -i "s/outfile =.*/outfile = p2n_"$out".codeml/" control_file_null.ctrl

		echo "Testing null hypothesis for $out..."
		codeml control_file_null.ctrl > "$out"_result_null.codeml
                w0=$(egrep -w 'w[0] =' "$out"_result_null.codeml | head -1 | awk '{print $5}')
                w1=$(egrep -w 'w[1] =' "$out"_result_null.codeml | head -1 | awk '{print $5}')
                w2=$(egrep -w 'w[2] =' "$out"_result_null.codeml | head -1 | awk '{print $5}')
                lnL=$(grep -w 'lnL' "$out"_result_null.codeml | head -1 | awk '{print $3}')

                echo -e "$out\t$w0\t$w1\t$w2\t$lnL" >> output_"$2"/results_null_final.tsv
		mv *null.codeml output_"$2"/null_results
		echo "Done"

		echo "Testing alternative hypothesis for $out..."
                sed -i "s/seqfile =.*/seqfile = p2n_"$out".fa/" control_file_alt.ctrl
                sed -i "s/outfile =.*/outfile = p2n_"$out".codeml/" control_file_alt.ctrl
		codeml control_file_alt.ctrl > "$out"_result_alt.codeml
		w0=$(egrep -w 'w[0] =' "$out"_result_alt.codeml | head -1 | awk '{print $5}')
		w1=$(egrep -w 'w[1] =' "$out"_result_alt.codeml | head -1 | awk '{print $5}')
		w2=$(egrep -w 'w[2] =' "$out"_result_alt.codeml | head -1 | awk '{print $5}')
		lnL=$(grep -w 'lnL' "$out"_result_alt.codeml | head -1 | awk '{print $3}')

		echo -e "$out\t$w0\t$w1\t$w2\t$lnL" >> output_"$2"/results_alt_final.tsv
		mv *alt.codeml output_"$2"/alt_results
		echo "Done"
		rm p2n_"$out".fa p2n_"$out".codeml
	else
		rm p2n_$out.fa
	fi
done < common_geneIDs.lst

rm rst rst1 rub lnf 2NG* 4fold.nuc common_geneIDs.lst mafft* merged_genes.fa

paste output_"$2"/results_null_final.tsv output_"$2"/results_alt_final.tsv > output_"$2"/results_selection.tsv

sed -i 1i"gene_id_null\tw0_null\tw1_null\tw2_null\tlnl_null\tgene_id_alt\tw0_alt\tw1_alt\tw2_alt\tlnl_alt" output_"$2"/results_selection.tsv
#}

Rscript --vanilla LRT_test.R ${2}
