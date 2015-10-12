#!/bin/bash
#
#                           Annotate rRNA 
#
#========================================================================================
#
#  Annotate rRNA genes based on HMMER3 profils.
#
#  go_rrna.sh <FASTAFILE>
#
#		- <FASTAFILE> : The fasta file containing the genome to annotate
#
#  Results are printed to the standart output
#
#========================================================================================

# -- CAUTION -- Works as long than the script 
#               is not called through a symlink
SCRIPT_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${SCRIPT_DIR}/../../../scripts/bash_init.sh"

pushTmpDir ORG.rrna

	RRNADB="${RRNA_DATA_DIR}/plast_RRNA.hmm"
	export CAUTRNADB
	
	if [[ ! "$1" =~ ^/ ]]; then
		QUERY="${CALL_DIR}/$1"
	else
		QUERY="$1"
	fi

	RRNA=$(basename ${QUERY})
	
	hmmsearch --max ${RRNADB} ${QUERY} | \
		awk '/Query: / { \
		                profil=$2; \
		                match($3,"[0-9][0-9]*");\
		                lprof=substr($3,RSTART,RLENGTH)} \
		     / [0-9][0-9]* ! / { \
		                print profil,lprof,$7,$8,$10,$11}' | \
		awk '($3 <=5) && (($2-$4) <=5) { \
		                full=1;$5=$5-$3+1;$6=$6+($2-$4)}  
		               {loc=$5".."$6} \
		     ($1 ~ /_RC$/) { \
		                loc="complement("loc")"} \
		               {match($1,"_..*S");\
		                rrna=substr($1,RSTART+1,RLENGTH-1);\
		                print "FT   rRNA            " loc; \
		                print "FT                   /gene=\"rrn"rrna"\""
		                print "FT                   /product=\""rrna" ribosomal RNA\""\
		                }'	

popTmpDir

exit 0
