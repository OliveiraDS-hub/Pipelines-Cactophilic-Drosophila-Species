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
    --coverage)
    COVERAGE=$2
    shift 2
    ;;
    --threads)
    THREADS=$2
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
  python remove_rep_seqs.py $CONSENSUS ## output: only_TEs_consensus.fa
}

Remove_duplicated_consensus_s2 () {
  python remove_dup_cons.py polished_TEs_s1.fa "round1_dedup.fa"
  python remove_dup_cons.py round1_dedup.fa "polished_TEs_s2.fa"
  rm round1_dedup*
}

Remove_tandem_rep_consensus_s3 () {
  ## Mask repeats in each consensus with TRF
  python trf_run.py polished_TEs_s2.fa

  ## Multi-line fasta to single-line fasta
  awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq = seq $0} END {if (seq) print seq}' polished_TEs_s2.fa.2.5.6.75.20.50.500.mask > polished_TEs_s2_masked_SL.fa

  ## Compute consensus length
  samtools faidx polished_TEs_s2_masked_SL.fa

  masked_TEIDs=$(grep '>' polished_TEs_s2_masked_SL.fa | sed 's/>//g')

  if [[ -f repeats-notTEs.lst ]]; then
    rm repeats-notTEs.lst
  fi

  set +e
  while read -r line; do
    TE_length=$(grep -w "$line" polished_TEs_s2_masked_SL.fa.fai | cut -f2)
    seq_TE=$(grep -A 1 "$line" polished_TEs_s2_masked_SL.fa | tail -1)
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
    python remove_seqs.py polished_TEs_s2.fa polished_TEs_s3.fa repeats-notTEs.lst
    rm repeats-notTEs.lst
  else
    echo "0 sequences with over 50% of tandem repeated sequences! Moving to the next filtering..."
    cat polished_TEs_s2.fa > polished_TEs_s3.fa
  fi
  rm *.50.500.* *_masked_SL.fa *fai
}

Remove_CDS-like_s4 () {
  echo -e "Removing consensus containing high similarity with CDSs (80% of cons. length)"
  blastn -subject "$CDS" -query polished_TEs_s3.fa -qcov_hsp_perc 80 -perc_identity 80 -out blast-CDS_s4.tmp -outfmt 6
  cut -f1 blast-CDS_s4.tmp | sort | uniq > cons-from-CDS_s4.tmp

  if [[ -s cons-from-CDS_s4.tmp ]]; then
    total_seq=$(wc -l cons-from-CDS_s4.tmp | awk '{print $1}')
    echo -e "$total_seq consensus with over 80% identity with CDS sequences! Removing from library..."
    python remove_seqs.py polished_TEs_s3.fa polished_TEs_s4.fa cons-from-CDS_s4.tmp
  else
    echo "0 consensus with over 80% identity with CDS sequences! Moving to the next filtering..."
    cat polished_TEs_s3.fa > polished_TEs_s4.fa
  fi
  rm *s4.tmp
}

TE_classification_s5 () {
  python TE_classification.py polished_TEs_s4.fa $DATABASE $SPECIES

  ## Clean data
  # rm *fai *s5.tmp #*.out

}

RepeatCraft_s6 () {
  echo "Masking genome with final TE library..."
  # RepeatMasker "$GENOME" -lib "$SPECIES"/polished_TEs_s5.fa -cutoff 250 -norna -gff -a -s -pa "$THREADS" 1> /dev/null
  echo "Done!"
  echo "LTR finder..."
  LTR_FINDER_parallel -seq "$GENOME" -threads "$THREADS" 1> /dev/null
  echo "Done!"
  echo -e "Parsing TE insertions with RepeatCraft..."
  cat "$PATH_TO_REPEAT_CRAFT"/example/repeatcraft_strict.cfg > config.cfg

  LTR_FINDER_out="${GENOME}.finder.combine.gff3"
  sed -i "s|ltr_finder_gff: None|ltr_finder_gff: "$LTR_FINDER_out"|" config.cfg

  set +e
  python "$PATH_TO_REPEAT_CRAFT"/repeatcraft.py \
  -r "$GENOME".out.gff \
  -u "$GENOME".out \
  -c config.cfg -o repcraft.out -m loose 2>/dev/null
  set -e
  egrep -v 'Simple_repeat|Low_complexity|Satellite' repcraft.out.rmerge.gff | sed 's/Tstart.*ID=//; s/;.*//g' > TEannot_RC_raw.gtf
  echo "done"
}

Short_insertions () {
  echo "Removing TEs < 80nt"
  awk 'BEGIN {FS=OFS="\t"} NR >=1 {print $0, $5 - $4}' TEannot_RC_raw.gtf | awk '$10 > 79' | cut -f1-9 > TEannot_RC_tmp.gtf

  echo "Adding TE number to each insertion"
  cut -f9 TEannot_RC_tmp.gtf > col9_tmp.rm
  counter="1"

  while read line; do
     echo -e "$line""@""$counter" >> col9.rm
     ((counter=counter+1))
  done < col9_tmp.rm
  gff_pt1=$(cut -f1-8 TEannot_RC_tmp.gtf)
  paste <(echo "$gff_pt1") col9.rm > TEannot_RC_counted.gtf; rm *.rm

  before=$(cut -f9 TEannot_RC_raw.gtf | sort | uniq | wc -l)
  after=$(cut -f9 TEannot_RC_counted.gtf | sort | uniq | wc -l)
  echo -e "Total consensus before filtering by length = $before\nTotal after filtering by length= $after\n"
  cd ../
}

Filter_SSR() {
  if [ ! -d 6_SSR ]; then
    mkdir 6_SSR; fi
  cd 6_SSR
  set +e
  ln -s ../"$GENOME" .
  ln -s ../5_repcraft/TEannot_RC_counted.gtf .
  echo "Identifying tandem repeats on TE insertions..."
  awk '{print $1,$4,$5,$9,$8,$7}' OFS='\t' TEannot_RC_counted.gtf > RMcounted.bed

  bedtools getfasta -fi "$GENOME" -bed RMcounted.bed -s -nameOnly | sed 's/(-)//g; s/(+)//g' > RMinsertions.fa
  trf RMinsertions.fa 2 5 6 75 20 50 500 -m -h >> /dev/null 2>&1
  echo "Done!"

  echo "Removing TE insertions rich in tandem repeats..."
  trf_file=$(basename *mask)
  echo "Mask file: $trf_file"
  sed -r -i '/^\s*$/d' "$trf_file"
  fasta_formatter -i "$trf_file" -o masked_TEs.fa

  samtools faidx RMinsertions.fa
  masked_TEIDs=$(grep '>' masked_TEs.fa | sed 's/>//g')

  set +e
  while read -r line; do
    TE_length=$(grep -w "$line" RMinsertions.fa.fai | cut -f2)
    seq_TE=$(grep -A 1 "$line" masked_TEs.fa | tail -1)
    masked_freq=$(grep -o 'N' <<< "$seq_TE" | grep -c .)
    if [ "$masked_freq" -gt "0" ]; then
      perc=$(( masked_freq*100/TE_length ))
      if [ "$perc" -gt "50" ]; then
        echo -e "$line" >> repeats-notTEs.lst
      fi
    fi
  done <<< "$masked_TEIDs"

  grep -v -w -f repeats-notTEs.lst TEannot_RC_counted.gtf | sed 's/@.*//g' > TEeasy_RMpolished.gtf

  python ../overlap_removal.py TEeasy_RMpolished.gtf TEeasy_RMpolished_final.gtf
  cd ../
}



#RepeatModeler2 function
#EarlGrey function
# Remove_nonTE_seqs_s1  ## OK v2
# Remove_duplicated_consensus_s2 ## OK v2 
# Remove_tandem_rep_consensus_s3 ## OK v2
# Remove_CDS-like_s4    ## OK v2
TE_classification_s5       ## OK v2
# RepeatCraft_s6    ## OK WORKING
# RepeatCraft        ## OK WORKING
# Short_insertions   ## OK WORKING
# Filter_SSR         ## OK WORKING






#
