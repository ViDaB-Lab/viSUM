#!/bin/bash

### Classification FILE PATHS ###
chunk=$1
fasta=$2  ##input fasta file
ct=$3  ##CAT contig2classification file with names
checkv=$4
vicat=$5  ## final classification file
vs2=$6 ## virsorter2 final-viral-score.tsv
gnmd=$7
ct3=$8 ##ct3 summary file
vitap=$9 ##vitap best_determined_lineages.tsv
### Functional annotation file paths ###
genomadannot="${10}"
### OTHER THINGS ###
ICTV="${11}"
THRESHOLD="${12}"  ##the amount of suport a call at t rank position would have to be above
#Program weighting for rank calls within script
VITAP="${13}"
genomad="${14}"
viCAT="${15}"
cenotetaker3="${16}"
CAT="${17}"
#####################################################
agsc="${18}"
#####################################################


#echo "---------------------------------------------------------------------------------"

#echo "VARIABLES FROM MAIN SCRIPT"
#echo ""$chunk" "$fasta" "$ct" "$checkv" "$vicat" "$vs2" "$gnmd" "$ct3" "$vitap" "$genomadannot" "$ICTV" "$THRESHOLD" "$VITAP" "$genomad" "$viCAT" "$cenotetaker3" "$CAT" "$agsc" ::: ./tmp_id_chunks/seq_ids_*"

for x in $(cat $chunk)
do	#echo "Contig ID = $x ; chunk = $chunk"
	base=$(basename $chunk)
	#echo "starting with genomad results"
	#assigning genom score--G--, virus=1, nothing=0; determining if classification useful 
	if grep -q -w "$x" $gnmd; 
	then	G="1" #assign 1 as genomad score
		#check if we can use classification by first checking if Unclassified
		geclass=$(grep -w "$x" $gnmd | awk -F "," '{print $2","$3","$4","$5","$6","$7","$8","$9","$10}')
      	
	else	#alright, not virus according to genomad, checking if plasmid. Will check diamond output for classification
		#echo "alright, not virus according to genomad, checking if plasmid. Will check diamond output for classification"
		G=0
		geclass=0
	fi
	#echo "starting with CAT results"
	#assigning CAT score--C--, virus=1, cellular=0; extracting classification if useful
	if grep -q -w "$x" $ct; then
    		line=$(grep -w "$x" $ct)
		C=1
    		cclass=$(echo $line | awk -F "," '{print $2","$3","$4","$5","$6","$7","$8","$9","$10}')
    		#echo "$cclass"
	else
    		C=0
		cclass="0"
	fi
	#echo "starting with vs2 results"
	#assigning vs2 score--V--, virus=1, not detected=0; not going to worry about what it was classified as
	if grep -q -w "$x" $vs2; then #check if seq present in virus file from vs2
		V=1
		vs="Detected"
	else
		V=0
		vs="Not detected"
	fi
	#echo "starting with CT3 results"
	#assigning CT# score--CT--, virus=1, not detected=0; will extract taxonomy for comparison
	if [[ $(grep -wc "$x" $ct3) -eq 1 ]]; then #check if seq present in virus summary file from ct3
		#echo "CT hit"
		CT=1
		ct3class=$(grep -w "$x" "$ct3" | awk -F "," '{print $2","$3","$4","$5","$6","$7","$8","$9","$10}')
	else
		CT=0
		ct3class=0
	fi
	#echo "starting with vitap results"
	#assigning vitap score--VT--, virus=1, not detected=0; not going to worry about what it was classified as
	if [[ $(grep -wc "$x" $vitap) -eq 1 ]]; then #check if seq present in virus file from vitap
		#echo "VITAP hit"
		VT=1
		vtclass=$(grep -w "$x" "$vitap" | awk -F "," '{print $2","$3","$4","$5","$6","$7","$8","$9","$10}')
		#echo $vtclass
	else
		VT=0
		vtclass=0
	fi
	#echo "starting with viCAT results"
	#assigning vitap score--VC--, virus=1, not detected=0; not going to worry about what it was classified as
	if [[ $(grep -wc "$x" $vicat) -eq 1 ]]; then 
		VC=1
		vcclass=$(grep -w "$x" $vicat | awk -F "\t" '{print "d__Viruses;"$4}' | sed 's/;/,/g')
		#echo $vcclass
	else
		VC=0
		vcclass=0
	fi
	#checkv info
	line2=$(grep -w "$x" "$checkv")
	chev=$(awk -F '\t' '{print $8}' <<<"$line2")
	len=$(awk -F '\t' '{print $2}' <<<"$line2")
	if [[ $chev == "Not-determined" ]]
	then CV=0
	elif [[ $chev == "Medium-quality" ]]
	then CV=0.75
	elif [[ $chev == "Low-quality" ]]
	then CV=0.5
	elif [[ $chev == "High-quality" || $chev == "Complete" ]]
	then CV=1.0
	fi
	#caclulate support score
	max=7
	score=$(echo "scale=2; ($G + $VT + $CT + $V + $VC + $C + $CV) / $max" | bc)
	rawscore=$(echo "scale=4; $G + $VT + $CT + $V + $VC + $C + $CV" | bc)
	if [[ $(echo "$rawscore > 0" | bc -l) -gt 0 ]]; then
		#Generate final taxonomy call
		echo "genomad,${geclass}" > tmp-"$x".z
		echo "cenotetaker3,${ct3class}" >> tmp-"$x".z
		echo "VITAP,${vtclass}" >> tmp-"$x".z
		echo "viCAT,${vcclass}" >> tmp-"$x".z
		echo "CAT,${cclass}" >> tmp-"$x".z
		grep -v ",0$" tmp-"$x".z > tmp-"$x".2 
		# right after you've built tmp.2…
		lines=$(wc -l < tmp-"$x".2)
		if (( lines > 1 )); then
			#echo "MULTIPLE"
			#echo "		mapfile -t out < <(
  			#python3 $agsc \
    			#--input     tmp-"$x".2 \
    			#--ictv      "$ICTV" \
    			#--threshold "$THRESHOLD" \
    			#--weights   VITAP="${VITAP}" viCAT="${viCAT}" genomad="${genomad}" cenotetaker3="${cenotetaker3}" CAT="${CAT}"
			#)"
			mapfile -t out < <(
  			python3 $agsc \
    			--input     tmp-"$x".2 \
    			--ictv      "$ICTV" \
    			--threshold "$THRESHOLD" \
    			--weights   VITAP="${VITAP}" viCAT="${viCAT}" genomad="${genomad}" cenotetaker3="${cenotetaker3}" CAT="${CAT}"
			)
			ftax=${out[0]}
			frank=${out[1]}
		elif (( lines < 1 )); then
			# No classifications to work with at all
    			ftax="d__Viruses,r__unclassified,k__unclassified,p__unclassified,c__unclassified,o__unclassified,f__unclassified,g__unclassified,s__unclassified"
    			frank="d__0.00;;r__0.00;;k__0.00;;p__0.00;;c__0.00;;o__0.00;;f__0.00;;g__0.00;;s__0.00"
		else
			# single-caller fallback
			ftax=$(awk -F"," '{print $2","$3","$4","$5","$6","$7","$8",g__unclassified,s__unclassified"}' tmp-"$x".2)
			frank="d__0.00;;r__0.00;;k__0.00;;p__0.00;;c__0.00;;o__0.00;;f__0.00;;g__0.00;;s__0.00"
		fi
		rm tmp-"$x".*
	else
		ftax="d__Viruses,r__unclassified,k__unclassified,p__unclassified,c__unclassified,o__unclassified,f__unclassified,g__unclassified,s__unclassified"
		frank="NA"
	fi
	#checking for funcation gene annotation from genomad or ct3
	if grep -qE "^${x}_[0-9]+" "$genomadannot"; then
		genes=$(awk -v C="$x" -F'\t' '$1 ~ "^"C"_[0-9]+$"{ gsub(/\t/,"::"); seq = (seq==""?$0:seq ";;"$0) } END{ print seq }' $genomadannot)
	else
		genes="NA"
	fi
	#echo "Summary for $x:"
	#echo "genomad("$G"),"$geclass""
	#echo "cenotetaker3 ("$CT"),"$ct3class""
	#echo "vitap("$VT"),"$vtclass""
	#echo "viCAT("$VC"),"$vcclass""
	#echo "CAT("$C"),"$cclass""
	#echo "virsorter("$V"), "$vs""
	#echo "CheckV("$CV") = "$chev""
	#echo "Score = $score"
	#echo "Raw score = $rawscore"
	#echo "Final tax = "$ftax""
	#echo "Final rank = "$frank""
	#echo "Genes = "$genes""
	#rm tmp-"$x".*
	if (( $(echo "$rawscore >= 1.0" | bc -l) )); then
		echo ""$x","$len","$score","$rawscore","$chev","$ftax","$frank","$genes"" >> tmp_"$base"_zclassified.csv 
	fi
	unset geclass ct3class vtclass vcclass cclass
	unset ftax frank domain realm kingdom phylum class order family good_genus species
	unset score rawscore G D V DI CV line2
	unset cnt_d cnt_r cnt_k cnt_p cnt_c cnt_o cnt_f
	unset genus_list total
	unset Ceve
done
