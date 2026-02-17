#!/bin/bash

# run by calling script followed by the name of a fasta file
# e.g., bash run_lca_pipeline.sh input.fasta

#####################################################
# ===== USER CONFIGURATION =====
### Classification FILE PATHS ###
fasta=$1  ##input fasta file
ct1="/work/alex/florida/RNA/assemblies/final_set/cat_325/FL_DT_rnaseq_virus_db_all_seqs_CAT_output_withnames2.txt"  ##CAT contig2classification file with names
checkv="/work/alex/florida/RNA/assemblies/final_set/checkv3_25/quality_summary.tsv"  ##checkv quality summary file
vicat="/work/alex/florida/RNA/assemblies/final_set/vat/FL_DT_rnaseq_virus_db_all_seqs_3_15_25_contigs_classified.tsv"  ## final classification file
vs2="/work/alex/florida/RNA/assemblies/final_set/vs2_325/final-viral-score.tsv" ## virsorter2 final-viral-score.tsv
gnmd1="/work/alex/florida/RNA/assemblies/final_set/genomad3_25/FL_DT_rnaseq_virus_db_all_seqs_3_15_25_summary/FL_DT_rnaseq_virus_db_all_seqs_3_15_25_virus_summary.tsv"
ct31="/work/alex/florida/RNA/assemblies/final_set/ct3_325/ct3_325_virus_summary.tsv" ##ct3 summary file
vitap1="/work/alex/florida/RNA/assemblies/final_set/vitap_325/best_determined_lineages.tsv" ##vitap best_determined_lineages.tsv
### Functional annotation file paths ###
genomadannot="/work/alex/florida/RNA/assemblies/final_set/genomad3_25/FL_DT_rnaseq_virus_db_all_seqs_3_15_25_summary/FL_DT_rnaseq_virus_db_all_seqs_3_15_25_virus_genes.tsv"
### OTHER THINGS ###
ICTV="/home/alex/programs/vSUM/bin/ICTV_masterspecieslist_5_11_25.csv"
thr="10"  ## number of threads
vsumlca="/home/alex/programs/vSUM//bin/vsum_lcaloopv0.3.sh"
#####################################################
### LCA options
THRESHOLD="0.5"  ##the amount of suport a call at t rank position would have to be above
#Program weighting for rank calls within script
VITAP="1.5"
genomad="1.3"
viCAT="1.25"
cenotetaker3="1.2"
CAT="1.0"
#####################################################
agsc="/home/alex/programs/vSUM/bin/lca_aggregatev0.1.py"
fixit="/home/alex/programs/vSUM/bin"
fallback="taxopy" #taxopy, taxopy needs path to taxdump or taxonkit
taxdump="/work/databases/IMG_VR/taxdump/"
#####################################################

# Create the list of sequence headers (scaffold names)
grep '^>' "$fasta" | sed 's/^>//' > tmp_seq.list
#confirm the fasta file 
echo "working with $fasta"
name=$(echo $fasta | awk -F ".fasta" '{print $1}')

#prepare outputs to match vSUM format
echo "Prepping classification files..."
echo "Starting  with VITAP..."
python ${fixit}/fixvitapv0.2.py --vitap ${vitap1}  --ictv ${ICTV} --out "${name}"_vitap_fixedz.csv --fallback ${fallback} --threads ${thr} --taxdump_dir ${taxdump}
cat "${name}"_vitap_fixedz.csv | sed 's/r__Naldaviricetes,/r__Naldaviricetes_U,/g' |sed 's/k__Naldaviricetes,/k__Naldaviricetes_U,/g' | sed 's/p__Naldaviricetes,/p__Naldaviricetes_U,/g' > "${name}"_vitap_fixed.csv
rm "${name}"_vitap_fixedz.csv
vitap=$(pwd)/""${name}"_vitap_fixed.csv"
echo "Next cenotetaker3..."
python ${fixit}/fixct3v0.1.py --ct3 ${ct31} --ictv ${ICTV} --out "${name}"_ct3_fixed.csv --threads ${thr} --fallback ${fallback} --taxdump_dir ${taxdump}
ct3=$(pwd)/""${name}"_ct3_fixed.csv"
echo "Next geNomad..."
python ${fixit}/fixgenomadv0.1.py --genomad ${gnmd1} --ictv ${ICTV} --out "${name}"_genomad_fixed.csv --threads ${thr} --fallback ${fallback} --taxdump_dir ${taxdump}
gnmd=$(pwd)/""${name}"_genomad_fixed.csv"
echo "Next CAT..."
grep -w "Viruses:" ${ct1} | sed 's/\t/,/g' > tempcat.csv
python ${fixit}/fixcatv0.1.py tempcat.csv ${ICTV} "${name}"_cat_fixed.csv
ct=$(pwd)/""${name}"_cat_fixed.csv"
rm tempcat.csv
#read -p "Pause here to inspect fixed files. Press Enter to continue or Ctrl+C to abort..." _

#Break lists into chunks based on threads variable
echo "splitting list into "$thr" chunks"
mkdir -p tmp_id_chunks
split -n l/"$thr" tmp_seq.list tmp_id_chunks/seq_ids_
find tmp_id_chunks -type f -size 0 -print -delete
ls tmp_id_chunks
# fire off one vsumlca per chunk, passing in the ORF and diamond files
echo "Starting classification loops now at $(date) across ${thr} threads... Could take a little, could not."
#echo " parallel command"
#echo "parallel --no-notice -j "$thr" "$vsumlca" {} "$fasta" "$ct" "$checkv" "$vicat" "$vs2" "$gnmd" "$ct3" "$vitap" "$genomadannot" "$ICTV" "$THRESHOLD" "$VITAP" "$genomad" "$viCAT" "$cenotetaker3" "$CAT" "$agsc" ::: ./tmp_id_chunks/seq_ids_*"

DEBUG=true parallel --no-notice -j "$thr" "$vsumlca" {} "$fasta" "$ct" "$checkv" "$vicat" "$vs2" "$gnmd" "$ct3" "$vitap" "$genomadannot" "$ICTV" "$THRESHOLD" "$VITAP" "$genomad" "$viCAT" "$cenotetaker3" "$CAT" "$agsc" ::: tmp_id_chunks/seq_ids_*

# pause here until all those jobs exit
wait

rm -r tmp_id_chunks
# now concatenate all the per-contig files
if ls *_zclassified.csv 1> /dev/null 2>&1; then
  echo "Found classified files"
  cat *_zclassified.csv >> "$name"_classified_vSUM.csv
  rm *_zclassified.csv
else
  echo "No classified files found"
fi
rm tmp_seq.list
echo "Script finished at $(date)"
