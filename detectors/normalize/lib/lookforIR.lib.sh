#!/bin/bash

source "${THIS_DIR}/../../../scripts/bash_init.sh"

SELECTIR="${PROG_DIR}/../../normalize/lib/selectIR.py"

function lookForIR {

	local QUERY="$1"
	local MATCHES=$(basename ${QUERY})
	      MATCHES="${MATCHES/.*/}.matches"
	
	local REPEATS="${MATCHES/.*/}.repseek"
	
	# Blast columns:
	# 	query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
	# We keep blast matches if : 
	#	The match is longer than 1000
	#   The identity is higher than 80%
	#
	# The match file has the following format:
	#	LSC/SSC  begin end  same_strand=1/diff_strand=0
	
	loginfo "Locating SSC and LSC by similarity..."
		blastn -db ${SCDB} \
		       -query ${QUERY} \
		       -outfmt 6 \
		       -max_target_seqs 10000  | \
		  awk -v id_match=80 -v lmin=100 \
		     '($4 > lmin) && (($3+0)>id_match) { 
		             SAME=(($7 < $8) && ($9 < $10)) || (($7 > $8) && ($9 > $10)); 
			 		 if ($7 < $8) 
			 			{print substr($2,1,3),$7,$8,SAME}  
			 		 else 
			 			{print substr($2,1,3),$8,$7,SAME}
			  }' | \
		  sort -nk 2 > ${MATCHES}
	loginfo "Done $(wc -l ${MATCHES} | awk '{print $1}') matches identified"

	loginfo "Looking for long inverted repeats..."
		repseek -c -p 0.001 -i ${QUERY} 2>> /dev/null > ${REPEATS}
		nrepeat="$(wc -l ${REPEATS} | awk '{print $1}')"
	loginfo "Done"

	if (( nrepeat == 0 )) ; then
		logwarning "No inverted repeat identified"
		return 1
	fi

	loginfo " --> ${nrepeat} repeats identified"

	loginfo "Marking and selecting the best inverted repeat..."
		local IR=( $(${SELECTIR} ${MATCHES} ${REPEATS}) )
	loginfo "Done"
	
	loginfo " --> IR size : IRa = ${IR[5]} /  IRb = ${IR[7]}"
	loginfo " --> IR Score: ${IR[8]}"
	
	deltaIR=$((IR[5] - IR[7]))
	
	if (( deltaIR < -10 )) ||  (( deltaIR > 10 )); then
		logwarning "Differences between IR lengths ($deltaIR) is greater than 10"
	fi
	
	
	echo "${IR[@]}"
}

SCDB="${IR_DATA_DIR}/SC_RefDB"

if [[ ! "$1" =~ ^/ ]]; then
	QUERY="${CALL_DIR}/$1"
else
	QUERY="$1"
fi
