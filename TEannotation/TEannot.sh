#!/bin/bash
set -e

function usage() {
   cat << help
Pipeline to polish TE consensus and improve classification of Unknown sequences
Oliveira D.S. Jun 2025

#Mandatory arguments:

  --genome             genome file (.fa) - the same used to build consensus
  --consensus          consensus file (from RepeatModeler2, EDTA, or EarlGrey) (.fa)
  --cds                reference CDS fasta file (.fa) - from your or closest related species
  --database           reference TE consensus (.fa) - i.e. Dfam or Repbase
  --annot              GTF/GFF with gene annotation
  --mate1              R1 of paired-end reads
  --mate2              R2 of paired-end reads
  --strand             RNA-seq strandness: either rf-stranded OR fwd-stranded
  --species            preffix for unclassifed consensus at the family level - i.e. dmel for D. melanogaster

#Optional arguments:
  --coverage           Minimum proportion of similarity between consensus and CDS to be removed (default: 80)
  --threads            threads for processing (default: 6)

help
}

PATH_TO_REPEAT_CRAFT="/home/oliveirads/softwares/repeatcraftp/" ## Set the path to repeatcraft software. From github, it must ends with "repeatcraftp"
COVERAGE="80"
THREADS="6"

while [[ $# -gt 0 ]]; do
    case $1 in
    --genome)
    GENOME=$2
    shift 2
    ;;
    --consensus)
    CONSENSUS=$2
    shift 2
    ;;
    --cds)
    CDS=$2
    shift 2
    ;;
    --database)
    DATABASE=$2
    shift 2
    ;;
    --annot)
    ANNOT=$2
    shift 2
    ;;
    --mate1)
    MATE1=$2
    shift 2
    ;;
    --mate2)
    MATE2=$2
    shift 2
    ;;
    --coverage)
    COVERAGE=$2
    shift 2
    ;;
    --threads)
    THREADS=$2
    shift 2
    ;;
    --strand)
    STRANDNESS=$2
    shift 2
    ;;
    --species)
    SPECIES=$2
    shift 2
    ;;
    -h | --help)
    usage
    exit 1
    ;;
    -*|--*)
    echo "Unknown option $1"
    exit 1
    ;;
    *)
    ARGS+=("$1")
    shift
    ;;
  esac
done

#RepeatModeler2 commandline
#EarlGrey commandline

Remove_nonTE_seqs_s1 () {
  if [[ ! -d "$SPECIES" ]]; then
    mkdir $SPECIES
  fi
  total_seq=$(grep -c '>' $CONSENSUS)
  echo -e "\n$CONSENSUS has $total_seq! Starting filtering..."
  echo -e "============== > STEP1: Removing non-TE sequences (tRNA, rRNA, Satellites, Low complexity)\n"
  python remove_rep_seqs.py $SPECIES $CONSENSUS ## output: only_TEs_consensus.fa
  if [[ -s "$SPECIES/polished_TEs_s1.fa" ]]; then
    total_seq=$(grep -c '^>' $SPECIES/polished_TEs_s1.fa)
    echo -e "Total consensus after STEP 1: $total_seq"
  fi
}

Remove_duplicated_consensus_s2 () {
  [ -s "$SPECIES/polished_TEs_s1.fa" ] || { echo "$SPECIES/polished_TEs_s1.fa missing or empty"; exit 1; }
  echo -e "\n\n============== > STEP2: Removing duplicated consensus\n"

  python remove_dup_cons.py $SPECIES/polished_TEs_s1.fa round1_dedup.fa
  python remove_dup_cons.py round1_dedup.fa $SPECIES/polished_TEs_s2.fa
  rm round1_dedup*
  if [[ -s "$SPECIES/polished_TEs_s2.fa" ]]; then
    total_seq=$(grep -c '^>' $SPECIES/polished_TEs_s2.fa)
    echo -e "Total consensus after STEP 2: $total_seq"
  fi
}

Remove_tandem_rep_consensus_s3 () {
  [ -s "$SPECIES/polished_TEs_s2.fa" ] || { echo "$SPECIES/polished_TEs_s2.fa missing or empty"; exit 1; }
  echo -e "\n\n============== > STEP3: Removing consensus comprised by > 50% tandem repeats...\n"

  ## Mask repeats in each consensus with TRF
  python trf_run.py $SPECIES/polished_TEs_s2.fa

  ## Multi-line fasta to single-line fasta
  awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq = seq $0} END {if (seq) print seq}' polished_TEs_s2.fa.2.5.6.75.20.50.500.mask > $SPECIES/polished_TEs_s2_masked_SL.fa

  ## Compute consensus length
  samtools faidx $SPECIES/polished_TEs_s2_masked_SL.fa

  masked_TEIDs=$(grep '>' $SPECIES/polished_TEs_s2_masked_SL.fa | sed 's/>//g')

  if [[ -f repeats-notTEs.lst ]]; then
    rm repeats-notTEs.lst
  fi

  set +e
  while read -r line; do
    TE_length=$(grep -w "$line" $SPECIES/polished_TEs_s2_masked_SL.fa.fai | cut -f2)
    seq_TE=$(grep -A 1 "$line" $SPECIES/polished_TEs_s2_masked_SL.fa | tail -1)
    masked_freq=$(grep -o 'N' <<< "$seq_TE" | grep -c .)
    if [ "$masked_freq" -gt "0" ]; then
      perc=$(( masked_freq*100/TE_length ))
      if [ "$perc" -gt "50" ]; then
        # echo -e "$line\t$perc"
        echo -e "$line" >> repeats-notTEs.lst
      fi
    fi
  done <<< "$masked_TEIDs"

  if [[ -s repeats-notTEs.lst ]]; then
    total_seq=$(wc -l repeats-notTEs.lst | awk '{print $1}')
    echo -e "$total_seq consensus with over 50% of tandem repeated sequences! Removing from library..."
    python remove_seqs.py $SPECIES/polished_TEs_s2.fa $SPECIES/polished_TEs_s3.fa repeats-notTEs.lst
    rm repeats-notTEs.lst
  else
    echo "0 sequences with over 50% of tandem repeated sequences! Moving to the next filtering..."
    cat $SPECIES/polished_TEs_s2.fa > $SPECIES/polished_TEs_s3.fa
  fi
  rm *.50.500.* $SPECIES/*_masked_SL.fa $SPECIES/*fai

  if [[ -s "$SPECIES/polished_TEs_s3.fa" ]]; then
    total_seq=$(grep -c '^>' $SPECIES/polished_TEs_s3.fa)
    echo -e "Total consensus after STEP 3: $total_seq"
  fi
}

Remove_CDS-like_s4 () {
  [ -s "$SPECIES/polished_TEs_s3.fa" ] || { echo "$SPECIES/polished_TEs_s3.fa missing or empty"; exit 1; }
  echo -e "\n\n============== > STEP4: Removing consensus containing high similarity with CDSs (80% of cons. length)...\n"

  blastn -subject "$CDS" -query $SPECIES/polished_TEs_s3.fa -qcov_hsp_perc 80 -perc_identity 80 -out blast-CDS_s4.tmp -outfmt 6
  cut -f1 blast-CDS_s4.tmp | sort | uniq > cons-from-CDS_s4.tmp

  if [[ -s cons-from-CDS_s4.tmp ]]; then
    total_seq=$(wc -l cons-from-CDS_s4.tmp | awk '{print $1}')
    echo -e "$total_seq consensus with over 80% identity with CDS sequences! Removing from library..."
    python remove_seqs.py $SPECIES/polished_TEs_s3.fa $SPECIES/polished_TEs_s4.fa cons-from-CDS_s4.tmp
  else
    echo "0 consensus with over 80% identity with CDS sequences! Moving to the next filtering..."
    cat $SPECIES/polished_TEs_s3.fa > $SPECIES/polished_TEs_s4.fa
  fi
  rm *s4.tmp

  if [[ -s "polished_TEs_s4.fa" ]]; then
    total_seq=$(grep -c '^>' polished_TEs_s4.fa)
    echo -e "Total consensus after STEP 4: $total_seq"
  fi
}

TE_classification_s5 () {
  [ -s "$SPECIES/polished_TEs_s4.fa" ] || { echo "$SPECIES/polished_TEs_s4.fa missing or empty"; exit 1; }
  echo -e "\n\n============== > STEP5: Classifying consensus based on $DATABASE...\n"

  python TE_classification.py $SPECIES/polished_TEs_s4.fa $DATABASE $SPECIES $THREADS

  ## Clean data
  rm *fai *s5.tmp classif_* *.out unclassified_*

  echo "Masking genome with classified TE library..."
  RepeatMasker "$GENOME" -lib "$SPECIES"/polished_TEs_s5.fa -cutoff 250 -norna -gff -a -s -pa "$THREADS" 1> /dev/null
  echo "Done!"
  rm *.cat *masked *tbl
  
  if [[ -s "polished_TEs_s5.fa" ]]; then
    unk_total_seq=$(grep -c -i 'unknown' polished_TEs_s4.fa)
    unk_total_seq_after=$(grep -c -i 'unknown' polished_TEs_s5.fa)

    echo -e "Total unknown consensus before STEP 5: $unk_total_seq"
    echo -e "Total unknown consensus after STEP 5: $unk_total_seq_after"
  fi
}

TE_strandness_s6 () {
  [ -s "$SPECIES/polished_TEs_s5.fa" ] || { echo "$SPECIES/polished_TEs_s5.fa missing or empty"; exit 1; }
  echo -e "\n\n============== > STEP6: Checking strandness of consensus based on stranded RNA-seq...\n"

  python TE_strandness.py $GENOME $ANNOT $THREADS $SPECIES $MATE1 $MATE2 $STRANDNESS
  if [[ -s "$SPECIES/polished_TEs_s6.fa" ]]; then
    total_seq=$(grep -c '^>' $SPECIES/polished_TEs_s6.fa)
    echo -e "Total consensus after STEP 6: $total_seq"
  fi
  mv *.out* "$SPECIES"
  rm "$GENOME".*
  cd "$SPECIES" && rm -rf *bed *bam *sam *bai classif_round3_done.fa_s5.tmp.out index && cd ..
}

RepeatCraft_s7 () {
  [ -s "$SPECIES/polished_TEs_s6.fa" ] || { echo -e "$SPECIES/polished_TEs_s6.fa missing or empty"; exit 1; }
  echo -e "\n\n============== > STEP7: RepeatCraft -> Merging TE insertions from the same TE family...\n"
  genome_fa=$(basename $GENOME)
  if [[ ! -f "$SPECIES/$genome_fa.out" ]]; then
    echo "Running RepeatMasker... It may take a while"
    RepeatMasker "$GENOME" -lib "$SPECIES"/polished_TEs_s6.fa -cutoff 250 -norna -gff -a -s -pa "$THREADS" -dir $SPECIES 1> /dev/null
    echo "Done!"
  fi
  
  echo "Running LTR finder..."
  LTR_FINDER_parallel -seq "$GENOME" -threads "$THREADS" #1> /dev/null
  LTR_FINDER_out="$genome_fa.finder.combine.gff3"
  echo "Done!"

  if [[ ! -f "$genome_fa.finder.combine.gff3" ]]; then
    echo "$genome_fa.finder.combine.gff3 doesn't exist! Exiting..."
    exit 1
  fi
  echo -e "Parsing TE insertions with RepeatCraft..."
  cat "$PATH_TO_REPEAT_CRAFT"/example/repeatcraft_strict.cfg > config.cfg

  sed -i "s|ltr_finder_gff: None|ltr_finder_gff: "$LTR_FINDER_out"|" config.cfg

  set +e
  python "$PATH_TO_REPEAT_CRAFT"/repeatcraft.py \
  -r "$SPECIES"/"$genome_fa".out.gff \
  -u "$SPECIES"/"$genome_fa".out \
  -c config.cfg -o repcraft.out -m loose 2>/dev/null
  set -e
  egrep -v 'Simple_repeat|Low_complexity|Satellite' repcraft.out.rmerge.gff | sed 's/Tstart.*ID=//; s/;.*//g' > $SPECIES/TEannot_S7.gtf

  if [[ -s "$SPECIES/TEannot_S7.gtf" ]]; then
    total_copies_before=$(tail -n +4 $SPECIES/$genome_fa.out | egrep -v 'Satellite|Simple_repeat|rRNA|Low_complexity|RNA|ARTEFACT' | wc -l)
    echo -e "Total TE copies before STEP 7: $total_copies_before"

    total_copies_after=$(wc -l $SPECIES/TEannot_S7.gtf | awk '{print $1}')
    echo -e "Total TE copies after STEP 7: $total_copies_after"
  fi
  rm config.cfg *.finder.* *list repcraft.out.* ltrfinder_*
}

Short_insertions_s8 () {
  [ -s "$SPECIES/TEannot_S7.gtf" ] || { echo -e"$SPECIES/TEannot_S7.gtf missing or empty"; exit 1; }
  echo "Removing TEs < 80nt"
  awk '$5 - $4 + 1 >= 80' $SPECIES/TEannot_S7.gtf > $SPECIES/TEannot_S8.gtf

  if [[ -s "$SPECIES/TEannot_S8.gtf" ]]; then
    total_copies_before=$(wc -l $SPECIES/TEannot_S7.gtf | awk '{print $1}')
    echo -e "Total TE copies before STEP 8: $total_copies_before"

    total_copies_after=$(wc -l $SPECIES/TEannot_S8.gtf | awk '{print $1}')
    echo -e "Total TE copies after STEP 8: $total_copies_after"
  fi
}

Filter_SSR_s9() {
  [ -s "$SPECIES/TEannot_S8.gtf" ] || { echo -e "$SPECIES/TEannot_S8.gtf missing or empty"; exit 1; }
  ## Add unique identifier to each insertion
  awk 'BEGIN{OFS="\t"} {$9 = $9 "_" NR; print}' $SPECIES/TEannot_S8.gtf | awk '{print $1,$4,$5,$9,$8,$7}' OFS='\t' > TEannot_S8_ID.bed
  echo "Identifying tandem repeats on TE insertions..."
  bedtools getfasta -fi "$GENOME" -bed TEannot_S8_ID.bed -s -nameOnly | sed 's/(-)//g; s/(+)//g' > RMinsertions.fa
  python trf_run.py RMinsertions.fa

  ## Multi-line fasta to single-line fasta
  awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq = seq $0} END {if (seq) print seq}' RMinsertions.fa.2.5.6.75.20.50.500.mask > RMinsertions_masked_SL.fa

  ## Compute consensus length
  samtools faidx RMinsertions_masked_SL.fa

  masked_TEIDs=$(grep '>' RMinsertions_masked_SL.fa | sed 's/>//g')

  if [[ -f repeats-notTEs.lst ]]; then
    rm repeats-notTEs.lst
  fi

  set +e
  while read -r line; do
    TE_length=$(grep -w "$line" RMinsertions_masked_SL.fa.fai | cut -f2)
    seq_TE=$(grep -A 1 "$line" RMinsertions_masked_SL.fa | tail -1)
    masked_freq=$(grep -o 'N' <<< "$seq_TE" | grep -c .)
    if [ "$masked_freq" -gt "0" ]; then
      perc=$(( masked_freq*100/TE_length ))
      if [ "$perc" -gt "50" ]; then
        echo -e "$line" >> repeats-notTEs.lst
      fi
    fi
  done <<< "$masked_TEIDs"

  if [[ -s repeats-notTEs.lst ]]; then
    total_seq=$(wc -l repeats-notTEs.lst | awk '{print $1}')
    echo -e "$total_seq copies with over 50% of tandem repeated sequences! Removing from annotation..."
    grep -v -f repeats-notTEs.lst $SPECIES/TEannot_S8.gtf > $SPECIES/TEannot_POLISHED_tmp.gtf
    ## Remove unique identifiers
    awk 'BEGIN{OFS="\t"} {$9 = gensub(/_[0-9]+$/, "", "g", $9); print}' "$SPECIES"/TEannot_POLISHED_tmp.gtf > "$SPECIES"/TEannot_POLISHED.gtf
    rm repeats-notTEs.lst $SPECIES/TEannot_POLISHED_tmp.gtf RMinsertions.fa TEannot_S8_ID.bed
  else
    echo "0 sequences with over 50% of tandem repeated sequences! Moving to the next filtering..."
    cat "$SPECIES/TEannot_S8.gtf > "$SPECIES/TEannot_POLISHED.gtf
  fi
  rm *.50.500.* *_masked_SL.fa *fai

  if [[ -s "$SPECIES/TEannot_POLISHED.gtf" ]]; then
    total_copies_before=$(wc -l $SPECIES/TEannot_S8.gtf | awk '{print $1}')
    echo -e "Total TE copies before STEP 9: $total_copies_before"

    total_copies_after=$(wc -l $SPECIES/TEannot_POLISHED.gtf | awk '{print $1}')
    echo -e "Total TE copies after STEP 9: $total_copies_after"
  fi

  ### Create GTF with family#order/class format
  awk 'BEGIN{OFS="\t"} {$9 = $9 "#" $3; print}' $SPECIES/TEannot_POLISHED.gtf > $SPECIES/TEannot_POLISHED_fullclass.gtf
  
  ### Create GTF with colors to IGV visualization
  awk 'BEGIN {OFS="\t"}
  {
    color = "#4073FF"  # default color
    if ($3 ~ /RC/)   color = "#ff6600"
    else if ($3 ~ /DNA/)  color = "#ff0000"
    else if ($3 ~ /LINE/) color = "#93C3EB"
    else if ($3 ~ /LTR/)  color = "#008000"
    $9 = "Target=" $9 ";color=" color
    print
  }' "$SPECIES/TEannot_POLISHED_fullclass.gtf" > "$SPECIES/TEannot_POLISHED_IGV.gtf"
}



# RepeatModeler2 function
# EarlGrey function
# Remove_nonTE_seqs_s1  ## OK v2
# Remove_duplicated_consensus_s2 ## OK v2 
# Remove_tandem_rep_consensus_s3 ## OK v2
# Remove_CDS-like_s4    ## OK v2
# TE_classification_s5       ## OK v2
# TE_strandness_s6    ## OK v2
RepeatCraft_s7    ## OK v2
Short_insertions_s8        ## OK v2
Filter_SSR_s9         ## OK 






#