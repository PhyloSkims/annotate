#!/bin/bash
#
#                           Annotate tRNA 
#
#========================================================================================
#
#  Annotate tRNA based on the Aragorn software predictions.
#
#  go_trna.sh <FASTAFILE>
#
#		- <FASTAFILE> : The fasta file containing the genome to annotate
#
#  Results are printed to the standart output
#
#========================================================================================

# -- CAUTION -- Works as long as the script 
#               is not called through a symlink
THIS_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${THIS_DIR}/../../../scripts/bash_init.sh"

pushTmpDir ORG.trna

	CAUTRNADB="${TRNA_DATA_DIR}/CAU_tRNA_DB.fasta"
	export CAUTRNADB
	
	if [[ ! "$1" =~ ^/ ]]; then
		QUERY="${CALL_DIR}/$1"
	else
		QUERY="$1"
	fi

	TRNA=$(basename ${QUERY})
	
	aragorn -i -w -seq -gcbact ${QUERY} | \
		${AwkCmd} -f ${PROG_DIR}/../lib/aragorn_wrapper.awk
	

popTmpDir

exit 0

