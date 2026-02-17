#!/bin/bash


##This is a script to assign a score to all virus sequence classficiations and generate a final taxonomy sheet boiiiii
genomadvir="/scratch1/FINAL_QC_9_28_22/results/taxonomy/polyA/genom/flsctld_polya_genomad_output/all_polya_virus_transcripts.okay_above500_summary/all_polya_virus_transcripts.okay_above500_virus_summary.tsv"
genomadpro="/scratch1/FINAL_QC_9_28_22/results/taxonomy/polyA/genom/flsctld_polya_genomad_output/all_polya_virus_transcripts.okay_above500_summary/all_polya_virus_transcripts.okay_above500_plasmid_summary.tsv"
diamondout="/scratch1/FINAL_QC_9_28_22/results/taxonomy/polyA/imrvr/flscld_polya_imrvr4_dmd.summary"
checkv="/scratch1/FINAL_QC_9_28_22/results/taxonomy/polyA/checkv/flsctld_polyA_checkv/quality_summary.tsv"
virso="/scratch1/FINAL_QC_9_28_22/results/taxonomy/polyA/virfin/flsctld_polyA_vs2.out/final-viral-boundary.tsv"
d6="/scratch1/FINAL_QC_9_28_22/results/taxonomy/polyA/d6/flsctld_polya_vtransresults.csv"
info="/scratch1/FINAL_QC_9_28_22/results/taxonomy/polyA/transcript.info" #just a two column csv with two columns 
tax="/home/ajv5/sctld/finaltaxsheet/tax.csv" ## keep this for now
output="flsctld_polya_final_taxsheet_unchecked.csv"
echo "transcriptid,newid,length,deep6score,genomadscore,virsorter2score,diamondscore,checkVscore,AbsScore,ConfidenceScore,classificationsource,classification,Realm,Kingdom,Phylum,Class,Order,Family,Genus" > $output

for x in $( cat $1)
do	echo "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
	#grab og transcript name
	realid=$(grep -w "$x" $info | awk -F "," '{print $1}')
	echo ""$realid" "$x" at $(date)"
	#assigning genom score--G--, virus=1, plasmid=0.5, nothing=0; determining if classification useful 
	grep -w "$x" $genomadvir > tmp.1 #check if seq present in virus file from genomad
	echo "wc -l tmp.1 = "$(wc -l tmp.1)""
	if [[ $(wc -l tmp.1 | awk '{print $1}') -eq 1 ]] #if present, wc -l tmp.1 should = 1
	then	G="1" #assign 1 as genomad score
		#check if we can use classification by first checking if Unclassified
		gecl=$(cat tmp.1 | awk -F "\t" '{print $11}')
		if [[ $(echo "$gecl" | grep -wc "Unclassified") -eq 1 ]]
		then	#ok its Unclassified, will keep geclass as Unclassified but signal we need to check diamond results for classification by asigning the value 0 to gc
			geclass="Unclassified"
			gc="0"
		else	#ok this is likely a good classification, now lets see if it is complete until genus and if not, we will fill in blanks
			#echo "ok this is likely a good classification, now lets see if it is complete until genus and if not, we will fill in blanks"
			check=$(echo $gecl | cut -d ';' -f 1- --output-delimiter=$'\n'| wc -l) 
			echo "$check is check"
			if [[ $check -eq 8 ]]
			then	#Ok it looks complete, will use this classification in final line
				echo "Ok it looks complete, will use this classification in final line"
				geclass=$(echo $gecl | sed 's/;/,/g')
				gc="1"
				#echo "$gecl"
				#echo "$geclass"
			else	#ok its not complete, lets see how many blanks we need to fill.
				need=$(echo "8-"$check"" | bc)
				#echo "need $need"
				i="1"
				while [ $i -le $need ]
				do	echo "unassigned" > tmp"$i".zz
					i=$(($i+1))
				done
				#ls | grep "tmp"
				paste -d',' tmp*.zz >tmp.zz
				needed=$(cat tmp.zz)
				gecla=$(echo $gecl | sed 's/;/,/g')
				geclass=""$gecla","$needed""
				rm *.zz
				gc="1"
			fi
		fi
      	else	#alright, not virus according to genomad, checking if plasmid. Will check diamond output for classification
		#echo "alright, not virus according to genomad, checking if plasmid. Will check diamond output for classification"
		grep -w "$x" $genomadpro > tmp.1
		if [[ $(wc -l tmp.1 | awk '{print $1}') -eq 1 ]]
		then	G="0.5"
			gc=2
		else	#alright, not virus nor plasmid according to genomad
			#echo "alright, not virus nor plasmid according to genomad"
			G=0
			gc=0
		fi
	fi
	rm tmp.1
	#assign deep6 score and prep annotation
	grep -w "$x" $d6 > tmp.2
	len=$( cat tmp.2 | awk -F "," '{print $2}')
	if [[ $len -ge 500 ]]
	then	D="1"
	else	D="0.5"
	fi
	#prepare deep6 classification
	re=$(cat tmp.2 | awk -F "," '{print $4}')
	realm=$(grep -w "$re" $tax | awk -F "," '{print $2}')
	if [[ $gc -eq 2 ]] 
	then	d6class="Plasmid,$realm,unassigned,unassigned,unassigned,unassigned,unassigned,unassigned"
	else	d6class="Virus,$realm,unassigned,unassigned,unassigned,unassigned,unassigned,unassigned"	
	fi
	rm tmp.2
	#assign virsorter score
	if [[ $(grep -wc "$x" $virso) -eq 1 ]]
	then	V=1
	else	V=0
	fi
	#assign diamond score and prepare classification
	grep -w "$x" $diamondout > tmp.3
	if [[ $(wc -l tmp.3 | awk '{print $1}') -eq 1 ]]
	then	DI=1
		#checking if we need diamond classification
		if [[ $gc != 1 ]]
		then	dc=$(cat tmp.3 | grep -c ",;;;;;;;,")
			if [[ $dc -eq 0 ]]
			then	dcla=$(cat tmp.3 | awk -F "," '{print $9}')
				#check for realm
				if [[ $(echo "$dcla" | awk -F ";" '{print $1}' | grep -c "r__") -eq 1 ]]
				then	drealm=$(echo "$dcla" | awk -F ";" '{print $1}' | sed 's/r__//g')
				else	#clasification is the dsDNA classes unassigned to realm in ICTV
					drealm="dsDNA"
				fi
				#check for kingdom
				if [[ $(echo "$dcla" | awk -F ";" '{print $2}' | grep -c "k__") -eq 1 ]]
				then	dking=$(echo "$dcla" | awk -F ";" '{print $2}' | sed 's/k__//g')
				else	dking="unassigned"
				fi
				#check for phylum
				if [[ $(echo "$dcla" | awk -F ";" '{print $3}' | grep -c "p__") -eq 1 ]]
				then	dphy=$(echo "$dcla" | awk -F ";" '{print $3}' | sed 's/p__//g')
				else	dphy="unassigned"
				fi
				#check for class
				if [[ $(echo "$dcla" | awk -F ";" '{print $4}' | grep -c "c__") -eq 1 ]]
				then	dclx=$(echo "$dcla" | awk -F ";" '{print $4}' | sed 's/c__//g')
				else	dclx="unassigned"
				fi
				#check for order
				if [[ $(echo "$dcla" | awk -F ";" '{print $5}' | grep -c "o__") -eq 1 ]]
				then	dord=$(echo "$dcla" | awk -F ";" '{print $5}' | sed 's/o__//g')
				else	dord="unassigned"
				fi
				#check for family
				if [[ $(echo "$dcla" | awk -F ";" '{print $6}' | grep -c "f__") -eq 1 ]]
				then	dfam=$(echo "$dcla" | awk -F ";" '{print $6}' | sed 's/f__//g')
				else	dfam="unassigned"
				fi
				#check for genus
				if [[ $(echo "$dcla" | awk -F ";" '{print $7}' | grep -c "g__") -eq 1 ]]
				then	dgen=$(echo "$dcla" | awk -F ";" '{print $7}' | sed 's/g__//g')
				else	dgen="unassigned"
				fi
				dclass=""$drealm","$dking","$dphy","$dclx","$dord","$dfam","$dgen""
			fi
		fi
	else	#ok no diamond alignment
		dc=1
		DI=0
	fi
	rm tmp.3
	#assigning checkv score
	chev=$(grep -w "$x" $checkv | awk -F "\t" '{print $8}')
	if [[ $chev == "Complete" ]]
	then CV=1.0
	elif [[ $chev == "High-quality" ]]
	then CV=0.75
	elif [[ $chev == "Medium-quality" ]]
	then CV=0.50
	elif [[ $chev == "Low-quality" ]]
	then CV=0.25
	elif [[ $chev == "Not-determined" ]]
	then CV=0
	fi
	#Selecting classification to use in final line
	if [[ $gc -eq 1 ]]
	then	CLASS="$geclass"
		source="genomad-HC"
	elif [[ $dc -eq 0 ]]
	then	if [[ $gc -eq 2 ]]
		then	CLASS="Plasmid,"$dclass""
			source="diamond-MC"
		else	CLASS="Virus,"$dclass""
			source="diamond-MC"
		fi
	else 	#only info is from deep6
		CLASS="$d6class"
		source="deep6-LC"
	fi
	#calculate confidence score
	max="5"
	score=$(echo "scale=2; ($G + $D + $V + $DI + $CV) / $max" | bc)
	rawscore=$(echo "scale=4; $G + $D + $V + $DI + $CV" | bc)
	#put it all together now
	line=""$realid","$x","$len","$D","$G","$V","$DI","$CV","$rawscore","$score","$source","$CLASS""
	echo "transcriptid,newid,length,deep6score,genomadscore,virsorter2score,diamondscore,checkVscore,AbsScore,ConfidenceScore,classificationsource,classification,Realm,Kingdom,Phylum,Class,Order,Family,Genus"
	echo "$line"
	echo "$line" >> $output
	echo "$x done at $(date)"
	G=z
	D=z
	V=z
	DI=z
	CV=z
	gc=z
	dc=z
done
