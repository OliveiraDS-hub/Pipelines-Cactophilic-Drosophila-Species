#!/bin/bash
set -e

function usage() {
   cat << help
Polishing of transposable elements annotation!
Oliveira D.S. 10/2023

#Mandatory arguments:

  --genome             genome file (.fa)
  --consensus          consensus file (from RepeatModeler2, EDTA) (.fa)
  --cds                reference CDS fasta file (.fa)
  --database           reference TE consensus (.fa)

#Optional arguments:
  --coverage           Minimum proportion of similarity between consensus and CDS to be removed (default: 50)
  --threads            threads for processing (default: 6)

help
}

COVERAGE="50"
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

Remove_CDS-like () {
  if [ ! -d 1_CDS ]; then
    mkdir 1_CDS; fi
  cd 1_CDS
  ln -s ../"$CDS" . ; ln -s ../"$CONSENSUS" .
  echo -e "Removing consensus containing high similarity with CDSs (80% of cons. length)"
  makeblastdb -in "$CDS" -dbtype nucl -out cds-seq 1> /dev/null
  blastn -db cds-seq -query "$CONSENSUS" -qcov_hsp_perc "$COVERAGE" -num_threads "$THREADS" -out blast-CDS.tsv -outfmt 6
  cut -f1 blast-CDS.tsv | sort | uniq > cons-from-CDS.lst
  grep '>' "$CONSENSUS" > consensus_IDs.lst
  grep -w -v -f cons-from-CDS.lst consensus_IDs.lst | sed 's/>//g' > filter1-non-CDS.lst

  seqtk subseq "$CONSENSUS" filter1-non-CDS.lst | awk '{print $1}' > S1_consensus.fa
  n_raw=$(grep -c '>' "$CONSENSUS")
  n_filt1=$(grep -c '>' S1_consensus.fa)
  echo -e "\nConsensus before filtering 1 = $n_raw\nConsensus after filtering 1 = $n_filt1\n"
  cd ../
}

TEclassifier () {
  set +e
  ####old masking function
  if [ ! -d 2_classification ]; then
    mkdir 2_classification; fi
  cd 2_classification
  echo -e "Masking consensus with $DATABASE..."
  ln -s ../1_CDS/S1_consensus.fa .
  ln -s ../"$DATABASE" .
  RepeatMasker S1_consensus.fa -lib "$DATABASE" -cutoff 225 -norna -a -s -par "$THREADS" 1> /dev/null
  echo "Done!"
  echo "Removing simple repeats, satellites and low complexity matches..."
  cat S1_consensus.fa.out | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | tail -n +4 | \
   egrep -i -v 'Simple_repeat|Low_complexity|Satellite' | cut -f5-7,10,11 | \
   awk 'BEGIN {FS=OFS="\t"} NR >=1 {print $0, $3 - $2}' > rm-consensus.tsv
  echo "Done!"
  samtools faidx S1_consensus.fa




  ######old overlaps function
  echo "Removing overlapped regions from TE insertions..."
  if [[ -f insertion_length.tsv ]]; then
    rm insertion_length.tsv; fi

  masked_cons=$(cut -f1 rm-consensus.tsv | sort | uniq)
  while read -r line; do
    ###Correcting overlaps of insertions
    TEcons_length=$(grep -w "$line" S1_consensus.fa.fai | cut -f2)

    freq=$(grep -w "$line" rm-consensus.tsv | cut -f5 | sort | uniq | wc -l) ## Col 5 = class
    if [[ "$freq" == 1 ]]; then
      rm_tab=$(grep -w "$line" rm-consensus.tsv | awk '{printf("%01d %s\n", NR, $0)}' | sed 's/ /\t/g')

      while IFS= read -r match_line; do
        ln_number=$(cut -f1 <<< "$match_line")
        ((ln_number=ln_number+1))
        end_coord=$(cut -f4 <<< "$match_line")
        next_insertion=$(awk -v next_line="$ln_number" '$1==next_line' <<< "$rm_tab")
        if [[ ! -z "$next_insertion" ]]; then
          while read -r compare_with; do
            next_start=$(cut -f3 <<< "$compare_with")
            if [[ "$end_coord" -lt "$next_start" ]]; then
              awk 'BEGIN {FS=OFS="\t"} NR >=1 {print $0, $4 - $3}' <<< "$match_line" | cut -f2-8 >> insertion_length.tsv
            else
              ((diff=end_coord - next_start + 1))
              ((new_end=end_coord - diff))
              start_coord=$(cut -f3 <<< "$match_line")
              ((length=new_end - start_coord + 1))
              echo -e "$match_line\t$length" | cut -f2-8 >> insertion_length.tsv
            fi
          done <<< "$next_insertion"
        else
          awk 'BEGIN {FS=OFS="\t"} NR >=1 {print $0, $4 - $3}' <<< "$match_line" | cut -f2-8  >> insertion_length.tsv
        fi
      done <<< "$rm_tab"
    else ### Multiple classes matching with consensus
      rm_tab=$(grep -w "$line" rm-consensus.tsv)
      n_match=$(grep -w "$line" <<< "$rm_tab" | cut -f5 | sort | uniq -c | awk '$1 > 1')

      if [[ ! -z $n_match ]]; then ###TE class with more than one match, but not multiple classes with single hit
        mult_match=$(grep -w "$line" <<< "$rm_tab" | cut -f5 | sort | uniq -c | awk '$1 > 1' | awk '{print $2}')
        while read -r multiple; do
          multiple_items=$(grep -w "$line" <<< "$rm_tab" | awk -v item="$multiple" '$5==item' | awk '{printf("%01d %s\n", NR, $0)}' | sed 's/ /\t/g')
          while IFS= read -r match_line; do
            ln_number=$(cut -f1 <<< "$match_line")
            ((ln_number=ln_number+1))
            end_coord=$(cut -f4 <<< "$match_line")
            next_insertion=$(awk -v next_line="$ln_number" '$1==next_line' <<< "$multiple_items")
            if [[ ! -z "$next_insertion" ]]; then
              while read -r compare_with; do
                next_start=$(cut -f3 <<< "$compare_with")
                if [[ "$end_coord" -lt "$next_start" ]]; then
                  awk 'BEGIN {FS=OFS="\t"} NR >=1 {print $0, $4 - $3}' <<< "$match_line" | cut -f2-8 >> insertion_length.tsv
                else
                  ((diff=end_coord - next_start + 1))
                  ((new_end=end_coord - diff))
                  start_coord=$(cut -f3 <<< "$match_line")
                  ((length=new_end - start_coord + 1))
                  echo -e "$match_line\t$length" | cut -f2-8 >> insertion_length.tsv
                fi
              done <<< "$next_insertion"
            else
              awk 'BEGIN {FS=OFS="\t"} NR >=1 {print $0, $4 - $3}' <<< "$match_line" | cut -f2-8 >> insertion_length.tsv
            fi
          done <<< "$multiple_items"
        done <<< "$mult_match"
      fi

      single_match=$(grep -w "$line" <<< "$rm_tab" | cut -f5 | sort | uniq -c | awk '$1 == 1' | awk '{print $2}')
      if [ ! -z "$single_match" ]; then #TE class with one hit, or multiple TE classes with single hits each
        while read -r single; do
          single_item=$(grep -w "$line" <<< "$rm_tab" | awk -v item="$single" '$5==item')
          awk 'BEGIN {FS=OFS="\t"} NR >=1 {print $0, $3 - $2}' <<< "$single_item" >> insertion_length.tsv
        done <<< "$single_match"
      fi
    fi
  done <<< "$masked_cons"
  echo "Done!"



  #old classifier function
  if [[ -f insertion_length_complete.tsv ]]; then
    rm insertion_length_complete.tsv; fi
  masked_cons=$(cut -f1 rm-consensus.tsv | sort | uniq)
  echo "Performing class and subfamily classification..."
  cut -f5 insertion_length.tsv | sed 's/_.*//g; s/-I//g; s/-LTR//g; s/JOCKEY2/Jockey-2/; s/Helitron1/Helitron-1/g; s/Gypsy1/Gypsy-1/; s/Gypsy2/Gypsy-2/; s/Gypsy4/Gypsy-4/; s/Gypsy5/Gypsy-5/; s/Gypsy6/Gypsy-6/; s/Gypsy7/Gypsy-7/; s/Copia2/Copia-2/; s/COPIA2I/Copia-2/' > families.lst
  paste insertion_length.tsv families.lst > insertion_length_complete.tsv


  ##### Family/Class level
  while read -r line; do
    TEcons_length=$(grep -w "$line" S1_consensus.fa.fai | cut -f2)
    freq_class=$(grep -w "$line" insertion_length_complete.tsv | cut -f5 | sed 's/_.*-/-/g; s/_.*//g' | sort | uniq | wc -l | awk '{print $1}')
    freq_subfamily=$(grep -w "$line" insertion_length_complete.tsv | cut -f4 | sort | uniq | wc -l | awk '{print $1}')

    if [ "$freq_class" -eq 1 ]; then
      total_masked=$(grep -w "$line" insertion_length_complete.tsv | awk '{sum+=$7;} END{print sum;}')
      class_masked=$(grep -w "$line" insertion_length_complete.tsv | cut -f5 | sort | uniq)
      perc_class=$(( total_masked*100/TEcons_length ))

      if [ "$freq_subfamily" -eq 1 ]; then
        family_masked=$(grep -w "$line" insertion_length_complete.tsv | cut -f4 | sort | uniq)
        total_masked_fam=$(grep -w "$line" insertion_length_complete.tsv | grep -w "$family_masked" | awk '{sum+=$7;} END{print sum;}')
        perc_class_fam=$(( total_masked_fam*100/TEcons_length ))
        echo -e "$line\t$class_masked\t$perc_class\t$family_masked\t$perc_class_fam" >> TE_classif.tsv
      else
        families=$(grep -w "$line" insertion_length_complete.tsv | cut -f4 | sort | uniq)
        while read -r family_TEs; do
          total_masked_fam=$(grep -w "$line" insertion_length_complete.tsv | grep -w "$family_TEs" | awk '{sum+=$7;} END{print sum;}')
          perc_class_fam=$(( total_masked_fam*100/TEcons_length ))
          echo "$family_TEs\t$perc_class_fam" >> perc_class.tsv
        done <<< "$families"
        family_higher_length=$(sort -Vr -k2,2 perc_class.tsv | head -1)
        rm perc_class.tsv
        echo -e "$line\t$class_masked\t$perc_class\t$family_higher_length" >> TE_classif.tsv
      fi

    else
      classes=$(grep -w "$line" insertion_length_complete.tsv | cut -f8 | sort | uniq)
      while read -r classes_TEs; do
        total_masked=$(grep -w "$line" insertion_length_complete.tsv | grep -w "$classes_TEs" | awk '{sum+=$7;} END{print sum;}')
        perc_class=$(( total_masked*100/TEcons_length ))
        echo -e "$classes_TEs\t$perc_class" >> perc_class.tsv
      done <<< "$classes"

      class_higher_length=$(sort -Vr -k2,2 perc_class.tsv | head -1)
      rm perc_class.tsv

      families=$(grep -w "$line" insertion_length_complete.tsv | cut -f4 | sort | uniq)
        while read -r family_TEs; do
          total_masked=$(grep -w "$line" insertion_length_complete.tsv | grep -w "$family_TEs" | awk '{sum+=$7;} END{print sum;}')
          perc_class=$(( total_masked*100/TEcons_length ))
          echo -e "$family_TEs\t$perc_class" >> perc_class.tsv
        done <<< "$families"
      family_higher_length=$(sort -Vr -k2,2 perc_class.tsv | head -1)
      rm perc_class.tsv
      echo -e "$line\t$class_higher_length\t$family_higher_length" >> TE_classif.tsv
      fi
  done <<< "$masked_cons"


  ## old renamer function
  cat S1_consensus.fa > TEclassif_consensus.fa
  new_IDs=$(awk '$3 > 50' TE_classif.tsv | awk -v OFS="#" '{print $4,$2}')
  old_IDs=$(awk '$3 > 50' TE_classif.tsv | cut -f1)
  paste <(echo "$old_IDs") <(echo "$new_IDs") > new_IDs.tsv

  repeated_families=$(cut -f2 new_IDs.tsv | sort | uniq -c | awk '$1 > 1' | grep '#')

  while read -r family_names; do
    freq=$(awk '{print $1}' <<< "$family_names")
    ID=$(awk '{print $2}' <<< "$family_names")
    pc1=$(sed 's/#/\t/' <<< "$ID" | cut -f1)
    pc2=$(sed 's/#/\t/' <<< "$ID" | cut -f2)
    counter="1"
    while [[ "$counter" -le "$freq" ]]; do
        new_ID=$(echo "${pc1}_${counter}#${pc2}")
        awk -v new_ID=$new_ID -v pattern=$ID '!x{x=sub(pattern,new_ID)}1' new_IDs.tsv > new_IDs.tmp && mv new_IDs.tmp new_IDs.tsv
        ((counter = counter + 1))
    done
  done <<< "$repeated_families"

  while read -r line; do
    old_ID=$(cut -f1 <<< "$line")
    new_ID=$(cut -f2 <<< "$line")
    sed -i -e "s|$old_ID|$new_ID|g" TEclassif_consensus.fa
  done < new_IDs.tsv
  n_total=$(wc -l new_IDs.tsv | awk '{print $1}')
  echo "TEclassifier has classifed $n_total TEs!"
  cd ../
}

FL_consensus () {
if [ ! -d 4_FL-insertions ]; then
  mkdir 4_FL-insertions; fi
echo -e "Removing consensus with less than 2 insertions with at least 80% of its length"
cd 4_FL-insertions
ln -s ../"$GENOME" .
ln -s ../2_classification/TEclassif_consensus.fa .
RepeatMasker "$GENOME" -lib TEclassif_consensus.fa -cutoff 250 -norna -gff -a -s -par "$THREADS" 1> /dev/null
cat "$GENOME".out | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | tail -n +4 | egrep -i -v 'Simple_repeat|Low_complexity|Satellite' > rm-consensus.tsv
cut -f10 rm-consensus.tsv | sort | uniq -c | awk '$1 >= 2' | awk '{print $2}' > IDs-post-filtering.tmp

while read -r TE_fam; do
 awk -v TE_family="$TE_fam" '$10==TE_family' rm-consensus.tsv | head -1 | cut -f10,11 | sed 's/\t/#/g' >> IDs_TE_families.tmp
done < IDs-post-filtering.tmp

cut -f10,11 rm-consensus.tsv | sed 's/\t/#/g' > family.tmp
paste rm-consensus.tsv family.tmp > tmp.tmp && mv tmp.tmp rm-consensus.tsv
seqtk subseq TEclassif_consensus.fa IDs_TE_families.tmp | awk '{print $1}' > filt2_consensus.fa.tmp
samtools faidx filt2_consensus.fa.tmp
cut -f1,2 filt2_consensus.fa.tmp.fai > tmp.tmp && mv tmp.tmp filt2_consensus.fa.tmp.fai

IDs=$(awk '{print $1}' filt2_consensus.fa.tmp.fai)
cov="80"
while read -r line; do
  cons_length=$(grep -w "$line" filt2_consensus.fa.tmp.fai | awk '{print $2}')
  grep -w "$line" rm-consensus.tsv > matches.tsv
  i=0

  while read -r match; do
        match_length=$(awk '{ $16 =  $7 - $6 } 1' <<< "$match" | awk '{print $16}')
        perc=$(( match_length*100/cons_length ))

          if [[ $perc -gt $cov ]]; then
            i=$((i+1))
            if [[ $i == "2" ]]; then
              echo "$line" >> IDs_reached.tmp
              break 1
            fi
          fi
  done < matches.tsv
done < IDs_TE_families.tmp

grep -w -f IDs_reached.tmp IDs_TE_families.tmp | sort | uniq > IDs-filt2.tmp
rm IDs_reached.tmp
seqtk subseq filt2_consensus.fa.tmp IDs-filt2.tmp | awk '{print $1}' > TElib_final.fa
n_f1=$(grep -c '>' TEclassif_consensus.fa)
n_f2=$(grep -c '>' TElib_final.fa)
grep -v -w -f IDs-filt2.tmp IDs_TE_families.tmp | sed 's/#.*//g' > IDs2remove.lst
echo -e "\nTotal consensus before filtering 2 = $n_f1\nTotal after filtering 2 = $n_f2\n"
rm *tmp *fai matches.tsv

}

Short_consensus () {
  samtools faidx TElib_final.fa
  sort -k2,2 -V TElib_final.fa.fai | grep -iv 'ltr' | awk '$2 > 200' | cut -f1 > long_non-LTR-cons.lst
  if [[ -f long_non-LTR-cons.lst ]]; then
    grep '>' TElib_final.fa | sed 's/>//g' | grep -i 'ltr' > all-LTR-cons.lst
    cat all-LTR-cons.lst long_non-LTR-cons.lst > long_consensus.lst
    if [[ -s long_consensus.lst ]]; then
      seqtk subseq TElib_final.fa long_consensus.lst > tmp.file && mv tmp.file TElib_final.fa
    fi
  fi
  cd ../
}

RepeatCraft () {
  if [ ! -d 5_repcraft ]; then
    mkdir 5_repcraft; fi
  cd 5_repcraft
  ln -s ../"$GENOME" .
  ln -s ../4_FL-insertions/TElib_final.fa .
  echo "Masking genome with final TE library..."
  RepeatMasker "$GENOME" -lib TElib_final.fa -cutoff 250 -norna -gff -a -s -par "$THREADS" 1> /dev/null
  echo "Done!"
  echo "LTR finder..."
  LTR_FINDER_parallel -seq "$GENOME" -threads 30 1> /dev/null
  echo "Done!"
  echo -e "Parsing TE insertions with RepeatCraft..."
  cat /home/oliveirads/softwares/repeatcraftp/example/repeatcraft_strict.cfg > config.cfg

  LTR_FINDER_out="${GENOME}.finder.combine.gff3"
  sed -i "s|ltr_finder_gff: None|ltr_finder_gff: "$LTR_FINDER_out"|" config.cfg

  set +e
  python /home/oliveirads/softwares/repeatcraftp/repeatcraft.py \
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

  grep -v -w -f repeats-notTEs.lst TEannot_RC_counted.gtf | sed 's/@.*//g' > TEeasy_RMpolished_final.gtf
  cd ../
}

Plot_landscape () {
  echo -e "Generating kimura distance..."
  cd 5_repcraft
  ### I SHOULD CLEAN THE ALIGN FILE!
  perl /home/oliveirads/softwares/RepeatMasker/util/calcDivergenceFromAlign.pl -s landscape.txt "$GENOME".align 2>/dev/null
  tail -n 72 landscape.txt > tmp && mv tmp landscape.txt

  g_size=$(grep 'total' *tbl | awk '{print $3}')
  Rscript --vanilla ../landscape_plot.R ${g_size} ${output}
}

#RepeatModeler2 function
#EarlGrey function
Remove_CDS-like    ## OK WORKING
TEclassifier       ## OK WORKING
FL_consensus       ## OK WORKING
Short_consensus    ## OK WORKING
RepeatCraft        ## OK WORKING
Short_insertions   ## OK WORKING
Filter_SSR         ## OK WORKING
Plot_landscape     ## OK WORKING






#
