source "${SCRIPT_DIR}/../../../scripts/bash_init.sh"

SELECTIR="${PROG_DIR}/../../normalize/lib/selectIR.py"

function lookForIR {

	local QUERY="$1"
	local MATCHES=$(basename ${QUERY})
	      MATCHES="${MATCHES/.*/}.matches"
	
	local REPEATS="${MATCHES/.*/}.repseek"
	
	loginfo "Locating SSC and LSC by similarity..."
		blastn -db ${SCDB} \
		       -query ${QUERY} \
		       -outfmt 6 \
		       -max_target_seqs 10000 | \
		  awk '($4 > 1000) && ($3>80) { \
		             SAME=(($7 < $8) && ($9 < $10)) || (($7 > $8) && ($9 > $10)); \
			 		 if ($7 < $8) \
			 			{print substr($2,1,3),$7,$8,SAME}  \
			 		 else \
			 			{print substr($2,1,3),$8,$7,SAME}}' | \
		  sort -nk 2 > ${MATCHES}
	loginfo "Done"
	  
	loginfo "Looking for long inverted repeats..."
		repseek -c -p 0.001 -i ${QUERY} 2>> /dev/null > ${REPEATS}
		loginfo " --> $(wc -l ${REPEATS} | awk '{print $1}') repeats identified"
	loginfo "Done"
	
	loginfo "Marking and selecting the best inverted repeat..."
		local IR=( $(${SELECTIR} ${MATCHES} ${REPEATS}) )
	loginfo "Done"
	
	loginfo " --> IR size : IRa = ${IR[5]} /  IRb = ${IR[7]}"
	loginfo " --> IR Score: ${IR[8]}"
	
	let "deltaIR=${IR[5]} -  ${IR[7]}"
	
	if (( $deltaIR < -10 )) ||  (( $deltaIR > 10 )); then
		logwarning "Differences between IR lengths ($deltaIR) is greater than 10"
	fi
	
	
	echo "${IR[@]}"
}

SCDB="${IR_DATA_DIR}/SC_RefDB"
QUERY="${CALL_DIR}/$1"


openLogFile "${QUERY/.*/}.log"
