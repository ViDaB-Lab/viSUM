#!/bin/bash

# run by calling script followed by the name of a fasta file
# e.g., bash run_lca_pipeline.sh input.fasta

#####################################################
# ===== USER CONFIGURATION =====
fasta=$1                          # input FASTA file
db=""			#IMR/VR database
dbmeta=""		#IMR/VR metadata
thr=""                           # number of threads
classify="classify_contig2.py"		#path to classify python script
f="0.5"			#this is the f parameter in the CAT calculation
mode="--fast"	#diamond sensitivity mode: --faster;--fast,--mid-sensitive,--sensitive,--more-sensitive,--very-sensitive,--ultra-sensitive
minscore="50" #minimum bitscore for a hit to be reported by diamond
minid="35" #minimum percent ID for alignment to be reported
minqcov="30" #minmum query coverage for alignment
minsubcov="30" #minmum subject coverage for alignment
evalue="1e-5" #evalue for diamond
tp="6"	#top percent of hits kep by diamond -- this is an important parameter when best = false in LCA script
c="4" #index chunks for diamond to increase efficiency; default is 4
b="2" #block size for diamond, default is 2
#####################################################

#SCRIPT BELOW
# Create the list of sequence headers (scaffold names)
grep '^>' "$fasta" | sed 's/^>//' > tmp_seq.list
#confirm the fasta file 
echo "working with $fasta"
name=$(echo $fasta | awk -F ".fasta" '{print $1}')

# Predict ORFs
echo "Predicting ORFs...."
mkdir tmp_chunks
seqkit split2 -p "$thr" -O tmp_chunks $fasta
cd ./tmp_chunks
ls *.fasta | parallel --no-notice -j"$thr" 'pyodigal-gv -p meta -i {} -a {.}.faa -o contigs.gff -q'
cat *.faa >> "$name"_orfs.faa
mv "$name"_orfs.faa ..
cd ..
rm -r tmp_chunks

# Align ORFs to IMG/VR database
echo "Aligning translated ORFs to IMR/VR...."		
diamond blastp --quiet --query "$name"_orfs.faa "$mode" --db $db --out tmp_"$name".diamond.aln --outfmt 6 qseqid qstart qend length sseqid bitscore evalue pident --query-cover $minqcov --subject-cover $minsubcov --id $minid --min-score $minscore --top $tp --evalue $evalue --threads $thr --index-chunks $c --block-size $b --compress 0 --no-self-hits

#Now running LCA classification per scaffold
# full paths to your ORFs and DIAMOND aln
orfs="$(pwd)/${name}_orfs.faa"
diamond="$(pwd)/"tmp_"$name".diamond.aln""
mkdir -p tmp_id_chunks
split -n l/"$thr" tmp_seq.list tmp_id_chunks/seq_ids_
find tmp_id_chunks -type f -size 0 -print -delete
# fire off one LCA.sh per chunk, passing in the ORF and diamond files
parallel --no-notice -j "$thr" ./LCA_loop.sh {} "$orfs" "$diamond" ::: tmp_id_chunks/seq_ids_*

# pause here until all those jobs exit
wait
rm -r tmp_id_chunks
# now concatenate all the per-contig files
if ls *_zclassified.tsv 1> /dev/null 2>&1; then
  echo "Found classified files"
  cat *_zclassified.tsv > "$name"_contigs_classified.tsv
  rm *_zclassified.tsv
else
  echo "No classified files found"
fi
mv "tmp_"$name".diamond.aln" ./"$name".diamond.aln
rm tmp_*
#For z in $(cat tmp_seq.list)
#do	orf_count=$(grep -c '^>"$z"_' "$x"_orfs.faa)
#	awk -F'\t' -v c="$z" '$1 ~ ("^"c"_")' tmp_"$x".diamond.aln > tmp_"${z}_hits.tsv"
#	# if no ORFs OR the hits file is empty ? skip
# 	if [ "$orf_count" -eq 0 ] || [ ! -s tmp_"${z}_hits.tsv" ]; then
#    		echo "No ORFs or hits for $z, moving on"
#    		continue
# 	fi
#	hit_count=$(wc -l < "tmp_${z}_hits.tsv")
#	# Extract sseqids, get taxonomies from IMG/VR metadata
#	awk -F "\t" '{print $5}' "${z}_hits.tsv" | awk -F "|" '{print $1}' > tmp_"$z"_hits.list
#	awk 'FNR==NR {ids[$1]; next} $1 in ids {print $16}' tmp_hits.list "$dbmeta" > tmp_"$z"_hits2.list
#	# Merge taxonomies back with alignment
#	paste tmp_"$x".diamond.aln tmp_hits2.list > tmp_"$z"_diamond_tax.tsv
#	# Get top hit per ORF
#	sort -k1,1 -k6,6nr tmp_"$z"_diamond_tax.tsv | awk '!seen[$1]++' > tmp_"$z"_top_hits.tsv
#	# Run contig-level LCA
#	LCA=$(python $classify tmp_"$z"_top_hits.tsv $f) 
#	# Print final summary for this contig
#	echo -e "$x\t$orf_count\t$hit_count\t$LCA" > "$z"_classified.tsv
#done

