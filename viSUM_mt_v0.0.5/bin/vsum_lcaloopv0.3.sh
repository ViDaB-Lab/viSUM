#!/bin/bash

### Classification FILE PATHS ###
chunk=$1
fasta=$2
ct=$3
checkv=$4
vicat=$5
vs2=$6
gnmd=$7
ct3=$8
vitap=$9
genomadannot="${10}"
ICTV="${11}"
THRESHOLD="${12}"
VITAP="${13}"
genomad="${14}"
viCAT="${15}"
cenotetaker3="${16}"
CAT="${17}"
agsc="${18}"

# Debugging and logging setup
DEBUG="${DEBUG:-false}"
LOGDIR="${LOGDIR:-debug_logs}"
mkdir -p "$LOGDIR"

process_sequence() {
  local x="$1"
  local ftax="" frank="" out=()
  local geclass="" cclass="" ct3class="" vtclass="" vcclass=""
  local G=0 C=0 V=0 CT=0 VT=0 VC=0 CV=0
  local chev="" len="" score="" rawscore="" genes="NA"

  local base=$(basename "$chunk")

  if grep -q -w "$x" "$gnmd"; then
    G=1
    geclass=$(grep -w "$x" "$gnmd" | awk -F "," '{print $2","$3","$4","$5","$6","$7","$8","$9","$10}')
  fi

  if grep -q -w "$x" "$ct"; then
    C=1
    cclass=$(grep -w "$x" "$ct" | awk -F "," '{print $2","$3","$4","$5","$6","$7","$8","$9","$10}')
  fi

  if grep -q -w "$x" "$vs2"; then
    V=1
  fi

  if [[ $(grep -wc "$x" "$ct3") -eq 1 ]]; then
    CT=1
    ct3class=$(grep -w "$x" "$ct3" | awk -F "," '{print $2","$3","$4","$5","$6","$7","$8","$9","$10}')
  fi

  if [[ $(grep -wc "$x" "$vitap") -eq 1 ]]; then
    VT=1
    vtclass=$(grep -w "$x" "$vitap" | awk -F "," '{print $2","$3","$4","$5","$6","$7","$8","$9","$10}')
  fi

  if [[ $(grep -wc "$x" "$vicat") -eq 1 ]]; then
    VC=1
    vcclass=$(grep -w "$x" "$vicat" | awk -F "\t" '{print "d__Viruses;"$4}' | sed 's/;/,/g')
  fi

  local line2=$(grep -w "$x" "$checkv")
  chev=$(awk -F '\t' '{print $8}' <<< "$line2")
  len=$(awk -F '\t' '{print $2}' <<< "$line2")

  case "$chev" in
    "Not-determined") CV=0 ;;
    "Medium-quality") CV=0.75 ;;
    "Low-quality") CV=0.5 ;;
    "High-quality"|"Complete") CV=1.0 ;;
  esac

  local max=7
  rawscore=$(echo "scale=4; $G + $VT + $CT + $V + $VC + $C + $CV" | bc)
  score=$(echo "scale=2; $rawscore / $max" | bc)

  if [[ $(echo "$rawscore > 0" | bc -l) -gt 0 ]]; then
    echo "genomad,${geclass}" > tmp-"$x".z
    echo "cenotetaker3,${ct3class}" >> tmp-"$x".z
    echo "VITAP,${vtclass}" >> tmp-"$x".z
    echo "viCAT,${vcclass}" >> tmp-"$x".z
    echo "CAT,${cclass}" >> tmp-"$x".z

    grep -v ",0$" tmp-"$x".z > tmp-"$x".2

    if [[ "$DEBUG" == "true" ]]; then
      cp tmp-"$x".2 "$LOGDIR/${x}_input.tsv"
    fi

    local lines=$(wc -l < tmp-"$x".2)

    if (( lines > 1 )); then
      mapfile -t out < <(
        python3 "$agsc" \
        --input     tmp-"$x".2 \
        --ictv      "$ICTV" \
        --threshold "$THRESHOLD" \
        --weights   VITAP="$VITAP" viCAT="$viCAT" genomad="$genomad" cenotetaker3="$cenotetaker3" CAT="$CAT"
      )
      ftax=${out[0]}
      frank=${out[1]}
    else
      ftax=$(awk -F, '{printf "%s,%s,%s,%s,%s,%s,%s,g__unclassified,s__unclassified\n", $2,$3,$4,$5,$6,$7,$8}' tmp-"$x".2)
      frank="NA"
    fi

    rm -f tmp-"$x".2 tmp-"$x".z
  else
    ftax="0"
    frank="NA"
  fi

  if grep -qE "^${x}_[0-9]+" "$genomadannot"; then
    genes=$(awk -v C="$x" -F'\t' '$1 ~ "^"C"_[0-9]+${"{ gsub(/\t/,"::"); seq = (seq==""?$0:seq ";;"$0) } END{ print seq }' "$genomadannot")
  fi

  if [[ "$DEBUG" == "true" ]]; then
    {
      echo "[$x] Raw support score: $rawscore"
      echo "[$x] Classification scores: G=$G VT=$VT CT=$CT V=$V VC=$VC C=$C CV=$CV"
      echo "[$x] Final taxonomy: $ftax"
      echo "[$x] Final rank confidence: $frank"
      echo "[$x] Genes: $genes"
    } >> "$LOGDIR/$x.log"
  fi

  if [[ $(echo "$rawscore > 0" | bc -l) -gt 0 ]]; then
    echo "$x,$len,$score,$rawscore,$chev,$ftax,$frank,$genes" >> tmp_"$base"_zclassified.csv
  fi
}

export -f process_sequence

while read -r x; do
  process_sequence "$x"
done < "$chunk"

