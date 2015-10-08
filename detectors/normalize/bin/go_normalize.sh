#!/bin/bash
#
#                           NORMALISATION D'UN PLASTIDE
#
#========================================================================================
# Ce programme dispose de 4 fonctions pour traiter les donnees fasta issues de genbank
# - seqlength : compte le nombre de paire de base du fichier
# ex : seqlength $1
#
# - cutseq : permet de couper un morceau de la sequence
# cutseq [x] [y]
# [x] : coordonne du debut de la sequence a couper
# [y] : coordonne de la fin de la sequence a couper
# ex : cutseq $1 10 100
#
#
# - revcomp : donne le brin reverse
# ex : $1 | revcomp
#
# - formatfasta : permet de coller a la suite plusieurs morceaux de sequence au moment de 
#   la reecriture
# - joinfasta : enleve les titres au moment de la reecriture du fichier et renvoie les 
#   informations dans la fonction formatfasta
# ex : joinfasta $1
#
#========================================================================================
# Pour lancer le programme, utiliser les commandes : 
# chmod +x normalize_plastid.sh
#./normalize_plastid.sh [fichier].fasta
#
# ex : seqlength $1
#
# cutseq $1 [x] [y]
# [x]:coordonne du debut [y]:coordonne de la fin de la sequence a couper
# ex : cutseq $1 10 100
#
# ex : $1 | revcomp
#
# ex : joinfasta $1
#========================================================================================

# -- CAUTION -- Works as long than the script 
#               is not called through a symlink
SCRIPT_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${SCRIPT_DIR}/../../../scripts/bash_init.sh"

function lookForIR {

	local QUERY="$1"
	local MATCHES=$(basename ${QUERY})
	      MATCHES="${MATCHES/.*/.matches}"
	
	local REPEATS="${MATCHES/.*/.repseek}"
	
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
		local IR=( $(${PROG_DIR}/selectIR.py ${MATCHES} ${REPEATS}) )
	loginfo "Done"
	
	loginfo " --> IR size : IRa = ${IR[5]} /  IRb = ${IR[7]}"
	loginfo " --> IR Score: ${IR[8]}"
	
	let "deltaIR=${IR[5]} -  ${IR[7]}"
	
	if (( $deltaIR < -10 )) ||  (( $deltaIR > 10 )); then
		logwarning "Differences between IR lengths ($deltaIR) is greater than 10"
	fi
	
	
	echo "${IR[@]}"
}

pushTmpDir ORG.normalize

	SCDB="${IR_DATA_DIR}/SC_RefDB"
	QUERY="${CALL_DIR}/$1"
	MATCHES="${1/.*/.matches}"
	REPEATS="${1/.*/.repseek}"

	tmpfasta1="tmp_$$_1.fasta"
	tmpfasta2="tmp_$$_2.fasta"


	openLogFile "${QUERY/.*/.log}"

	loginfo "Computing the genome size..."
		genome_length=$(seqlength $QUERY)
		loginfo " --> $genome_length bp"
	loginfo "Done"
	
	IR=( $(lookForIR ${QUERY}) )
	
	posIR1=${IR[4]}
	posIR2=${IR[6]}
	
	let "lenIR= ( ${IR[5]} +  ${IR[7]} ) / 2 " 

	let "endIR2=$posIR2 + $lenIR - 1"
	let "endIR1=$posIR1 + $lenIR - 1"
	
	if (( "$endIR2" >= "$genome_length" )) ; then	
		loginfo "IRB is at the end of the original sequence"
		
		#
		# We just move the IRB at the begining of the sequence
		#
		
		# Extract the IRB sequence
		let "posCut=($endIR1+$posIR2)/2"
		cutseq ${QUERY} ${posCut} ${genome_length} > ${tmpfasta1}

		# Append the remaining part of the genome		
		let "posCut=$posCut-1"
		cutseq ${QUERY} 1 ${posCut} >> ${tmpfasta1}
		
		# merges both the parts
		joinfasta ${tmpfasta1} > ${tmpfasta2}
		rm -f ${tmpfasta1}
		QUERY=${tmpfasta2}

		loginfo "Recompute location of the IR..."
			declare -a IR=( $(lookForIR ${QUERY}) )
		loginfo "Done"
		
		posIR1="${IR[4]}"
		posIR2="${IR[6]}"
		
		let "lenIR=(${IR[5]} +  ${IR[7]}) / 2 " 
	
		let "endIR2=$posIR2 + $lenIR - 1"
		let "endIR1=$posIR1 + $lenIR - 1"
		
	fi		
	
	tmpIR1="tmp_$$_IR1.fasta"		
	tmpIR2="tmp_$$_IR2.fasta"		
	
	#enregistre les deux fragments IRa et IRb complet
	cutseq ${QUERY} ${posIR1} ${endIR1} > ${tmpIR1}
	cutseq ${QUERY} ${posIR2} ${endIR2} > ${tmpIR2}
	
	let "lenSC1=$posIR1 -1 + ($genome_length - endIR2)"
	let "lenSC2=$posIR2 - $endIR1"
	
	center="${IR[0]}"
		
	tmpLSC="tmp_$$_LSC.fasta"		
	tmpSSC="tmp_$$_SSC.fasta"		
	
	# Extract the first SC present in between the two IRs
	# considering it as LSC

	let "beginLSC=$endIR1+1"
	let "endLSC=$posIR2-1"
	cutseq ${QUERY} ${beginLSC} ${endLSC} > ${tmpLSC}

	strandLSC="${IR[1]}"


	# Extract the second SC present in two parts
	# Considering it as SSC
	
	let "beginSSC=$endIR2+1"
	cutseq ${QUERY} ${beginSSC} ${genome_length} > ${tmpSSC}

	let "endSSC=$posIR1-1"
	cutseq ${QUERY} 1 ${endSSC} >> ${tmpSSC}

	joinfasta ${tmpSSC} > ${tmpfasta1}
	mv ${tmpfasta1} ${tmpSSC}
	
	strandSSC="${IR[3]}"
	
	

	if [[ "$center" == "SSC" ]]; then
	
		# Actually this is the oposite LSC is SSC and SSC is LSC

		# Exchange the SSC and LSC sequences
		mv ${tmpSSC}    ${tmpfasta1}
		mv ${tmpLSC}    ${tmpSSC}
		mv ${tmpfasta1} ${tmpLSC}
		
		# Exchange the IRa and IRb sequences
		mv ${tmpIR1}    ${tmpfasta1}
		mv ${tmpIR2}    ${tmpIR1}
		mv ${tmpfasta1} ${tmpIR2}
		
		tmp=${strandSSC}
		strandSSC=${strandLSC}
		strandLSC=${tmp}
		
	fi
	
	# Reverse complement the SSC if needed
	if [[ "${strandSSC}" == "-" ]]; then
		fastarevcomp -f ${tmpSSC} > ${tmpfasta1}
		mv ${tmpfasta1} ${tmpSSC}
	fi
	
	# Reverse complement the LSC if needed
	if [[ "${strandLSC}" == "-" ]]; then
		fastarevcomp -f ${tmpLSC} > ${tmpfasta1}
		mv ${tmpfasta1} ${tmpLSC}
	fi
	
	# Merges the four parts of the genome.
	cat ${tmpSSC} ${tmpIR1} ${tmpLSC} ${tmpIR2} | joinfasta

	
	
popTmpDir

exit 0

