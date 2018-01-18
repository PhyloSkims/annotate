#!/bin/bash
#
#           Annotate the Intergenic Spacer (ITS) of nuclear rDNA cluster
#
#========================================================================================
#
# This script is based on ITSx
# 
#
#  go_its.sh <FASTAFILE>
#
#		- <FASTAFILE> : The fasta file containing the cluster to annotate
#
#  Results are printed to the standart output
#
#========================================================================================

# -- CAUTION -- Works as long than the script 
#               is not called through a symlink
THIS_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${THIS_DIR}/../../../scripts/bash_init.sh"

pushTmpDir ORG.its

	loginfo "Normalizing nuclear rDNA cistron..."
	RRNADB="${NUCRRNA_DATA_DIR}/plants/nuc_RRNA.hmm"
	
	if [[ ! "$1" =~ ^/ ]]; then
		QUERY="${CALL_DIR}/$1"
	else
		QUERY="$1"
	fi

	loginfo "Sequence length $(seqlength ${QUERY})"

	strand=( $(hmmsearch --max ${RRNADB} ${QUERY} | \
				$AwkCmd '/Query: / { \
				                profil=$2; \
				                match($3,"[0-9][0-9]*");\
				                lprof=substr($3,RSTART,RLENGTH)} \
				     / [0-9][0-9]* ! / { \
				                print profil,lprof,$7,$8,$10,$11}' | \
				$AwkCmd '($3 <=5) && (($2-$4) <=5) { \
				                full=1;$5=$5-$3+1;$6=$6+($2-$4)}  \
				               {loc="Forward"} \
				     ($1 ~ /_RC$/) { \
				                loc="Reverse"} \
				     (full==1) {match($1,"_..*S");\
				                rrna=substr($1,RSTART+1,RLENGTH-1);\
				                print loc;\
				                full=0
				                }' | sort | uniq) )
				               
	if [[ "${#strand[@]}" == 1 ]] ; then
		if [[ "${strand[0]}" == "Forward" ]] ; then
			cat ${QUERY}
		else
			loginfo "Revert complement rDNA cluster cistron"
			fastarevcomp  ${QUERY} 
		fi
	else
		logerror "Cannot determine the Cistron orientation"
		exit 1
	fi
	 
    loginfo "Done."

popTmpDir

exit 0

