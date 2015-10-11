#!/bin/bash
#
#

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
source "${SCRIPT_DIR}/scripts/bash_init.sh"

pushTmpDir ORG.organnot

	if [[ ! "$1" =~ ^/ ]]; then
		QUERY="${CALL_DIR}/$1"
	else
		QUERY="$1"
	fi

	RESULTS=$(basename ${QUERY/.*/})
	LOG="${CALL_DIR}/${RESULTS}.log"
	
	rm -f ${LOG}
	openLogFile ${LOG}

	${PROG_DIR}/detectors/normalize/bin/go_normalize.sh ${QUERY} > "${RESULTS}.norm.fasta"
	
	${PROG_DIR}/detectors/ir/bin/go_ir.sh ${QUERY} > "${RESULTS}.annot"		
	${PROG_DIR}/detectors/trna/bin/go_trna.sh ${QUERY} >> "${RESULTS}.annot"
	
	cat "${RESULTS}.annot"
	
popTmpDir

