#!/bin/bash
#
#                           NORMALISATION D'UN PLASTIDE
#
#========================================================================================
#
# Normalize the way the chloroplaste genome sequence is linearized in the fasta file
# The normalized sequence is: 
#
#               LSC + IRB + SSC + IRA
#
# The SSC and LSC are approximatively mapped by similarity with a reference database
# Inverted repeats (IRs) are identified for maximizing the segregation between 
# LSC and SSC match 
# 
#
#  go_normalize.sh <FASTAFILE>
#
#		- <FASTAFILE> : The fasta file containing the genome to normalize
#
#  Results are printed to the standart output
#
#========================================================================================

# -- CAUTION -- Works as long than the script 
#               is not called through a symlink
THIS_DIR="$(dirname ${BASH_SOURCE[0]})"
source ${THIS_DIR}/../lib/lookforIR.lib.sh

ORG_DEBUG=1

pushTmpDir ORG.normalize

	tmpfasta1="tmp_$$_1.fasta"
	tmpfasta2="tmp_$$_2.fasta"

	logdebug "Running on : $QUERY"

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

		loginfo "Recomputing location of the IR..."
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
	cat ${tmpLSC} ${tmpIR2} ${tmpSSC} ${tmpIR1} | joinfasta

	
	
popTmpDir

exit 0

