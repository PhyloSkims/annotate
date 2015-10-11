#!/bin/bash
#
#                           Annotate tRNA 
#
#========================================================================================
#
#  Annotate tRNA based on the Aragorn software predictions.

#  go_trna.sh <FASTAFILE>
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

pushTmpDir ORG.trna

	CAUTRNADB="${TRNA_DATA_DIR}/CAU_tRNA_DB.fasta"
	export CAUTRNADB
	
	if [[ ! "$1" =~ ^/ ]]; then
		QUERY="${CALL_DIR}/$1"
	else
		QUERY="$1"
	fi

	TRNA=$(basename ${QUERY})
	
	aragorn -i -w -seq ${QUERY} | \
		${PROG_DIR}/../lib/aragorn_wrapper.awk
	

popTmpDir

exit 0

