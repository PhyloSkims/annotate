#!/bin/bash
#
#                      Annotate the Inverted Repeats of a plastide genome
#
#========================================================================================
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
SCRIPT_DIR="$(dirname ${BASH_SOURCE[0]})"
source ${SCRIPT_DIR}/../../normalize/lib/lookforIR.lib.sh

pushTmpDir ORG.ir

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

	beginLSC=1
	let "endLSC=$posIR1-1"


	let "beginSSC=$endIR1+1"
	let "endSSC=$posIR2-1"
	
	
	echo "FT   misc_feature    ${beginLSC}..${endLSC}"
	echo "FT                   /note=\"large single copy region (LSC)\""
	echo "FT   repeat_region   ${posIR1}..${endIR1}"
	echo "FT                   /rpt_type=INVERTED"
	echo "FT                   /note=\"left inverted repeat B; IRB\""
	echo "FT   misc_feature    ${beginSSC}..${endSSC}"
	echo "FT                   /note=\"small single copy region (SSC)\""
	echo "FT   repeat_region   ${posIR2}..${endIR2}"
	echo "FT                   /rpt_type=INVERTED"
	echo "FT                   /note=\"left inverted repeat A; IRA\""
	

popTmpDir

exit 0