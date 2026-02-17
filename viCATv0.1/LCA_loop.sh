#!/bin/bash
# LCA.sh <chunk_file>
chunk="$1"

# make sure these env variables or paths are exported before you call parallel:
#   name       ? fasta/orfs prefix
#   dbmeta     ? path to your IMG/VR metadata
#   classify   ? path to your Python classifier
#   f          ? bitscore fraction threshold
dbmeta=""	#database metadata
classify="classify_contig5.py"		#path to classify python script
f="0.5"			#this is the f parameter in the CAT calculation
orfs="$2"
diamond="$3"
best="false"

while IFS=$'\n' read -r z; do
  # count ORFs for this contig
  orf_count=$(grep -c "${z}_" "$orfs")

  # extract all hits for those ORFs
  awk -F'\t' -v c="$z" '$1 ~ ("^"c"_")' "$diamond" > tmp_${z}_hits.tsv

  # skip if no ORFs or no hits
  if [ "$orf_count" -eq 0 ] || [ ! -s tmp_${z}_hits.tsv ]; then
    echo "No ORFs or hits for $z, skipping"
    rm tmp_${z}_*
    continue
  fi

  # count hits
  hit_count=$(wc -l < tmp_${z}_hits.tsv)

  # map IMG/VR IDs ? taxonomy
  awk -F'\t' '{print $5}' tmp_${z}_hits.tsv | cut -d'|' -f1 > tmp_${z}_hits.list
  awk -F, 'FNR==NR{ids[$1]; next} $1 in ids{print $2}' "tmp_${z}_hits.list" "$dbmeta" > tmp_${z}_hits2.list

  # merge, pick top per ORF, classify
  paste "tmp_${z}_hits.tsv" "tmp_${z}_hits2.list" > "tmp_${z}_diamond_tax.tsv"
  if [ "$best" = "true" ]; then
  	sort -t$'\t' -k1,1 -k6,6nr "tmp_${z}_diamond_tax.tsv" | awk '!seen[$1]++' > "tmp_${z}_top_hits.tsv"
	LCA=$(python "$classify" "tmp_${z}_top_hits.tsv" "$f" )
  else 
	sort -t$'\t' -k1,1 -k6,6nr "tmp_${z}_diamond_tax.tsv" > "tmp_${z}_top_hits.tsv"
  	LCA=$(python "$classify" --normalize "tmp_${z}_top_hits.tsv" "$f")
  fi

  LCA=$(python "$classify" "tmp_${z}_top_hits.tsv" "$f")

  # write summary
  prefix=$(basename "$chunk")	
  echo -e "${z}\t${orf_count}\t${hit_count}\t${LCA}" >> "${prefix}_zclassified.tsv"
  rm tmp_${z}_*
done < "$chunk"
