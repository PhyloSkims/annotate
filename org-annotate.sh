#!/bin/bash
#
#
#
#                           Annotate Organelle 
#
#========================================================================================
#
#
#========================================================================================

# -- CAUTION -- Works as long as the script 
#               is not called through a symlink
THIS_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${THIS_DIR}/scripts/bash_init.sh"

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

	loginfo "Normalizing the structure of the Chloroplast sequence..."
		loginfo "   LSC + IRB + SSC + IRA"
		${PROG_DIR}/detectors/normalize/bin/go_normalize.sh ${QUERY} > "${RESULTS}.norm.fasta"
	loginfo "Done."
	
	loginfo "Annotating the Inverted repeats and Single copies (LSC and SSC)..."
		${PROG_DIR}/detectors/ir/bin/go_ir.sh ${QUERY} > "${RESULTS}.annot"		
	loginfo "Done."
	
	loginfo "Annotating the tRNA..."
		${PROG_DIR}/detectors/trna/bin/go_trna.sh ${QUERY} >> "${RESULTS}.annot"
	loginfo "Done."
	
	loginfo "Annotating the rRNA genes..."
		${PROG_DIR}/detectors/rrna/bin/go_rrna.sh ${QUERY} >> "${RESULTS}.annot"
	loginfo "Done."

	loginfo "Annotating the CDS..."
		${PROG_DIR}/detectors/cds/bin/go_cds.sh ${QUERY} >> "${RESULTS}.annot"
	loginfo "Done."
	
	loginfo "Printing annotations header..."
		echo "XX"
    	echo "FH   Key             Location/Qualifiers"
	loginfo "Done."
	
	loginfo "Ordering annotations..."
		awk '/^.....(misc|repeat|rRNA|tRNA|gene)/ { \
		        match($3,"[0-9][0-9]*"); \
		        pos=substr($3,RSTART,RLENGTH)*1000 + 1; \
		        print pos,$0;    \
		        next} \
		      { pos++; \
		        print pos,$0}' "${RESULTS}.annot" | \
		sort -nk1 |\
		awk '{ \
		        match($0,"^[0-9]* ");\
		        line=substr($0,RLENGTH+1);\
		        print line}' 
	loginfo "Done."
	
	loginfo "Closing annotations table..."
		echo "XX"
	loginfo "Done."
	
	loginfo "Computing statistics on nucleotide usage..."
		awk '! /^>/ { \
			    seq=toupper($0); \
				gsub(" ","",seq); \
			    lseq=length(seq); \
				for (i=0; i < lseq; i++) { \
					freq[substr(seq,i,1)]++}\
					} \
			 END { \
			 	other=0; \
			 	for (i in freq) { \
			 		if (i!="A" && i!="C" && i!="G" && i!="T") {\
			 			other+=freq[i] \
			 			} \
			 		}; \
			 		print "SQ   Sequence "\
			 		      (freq["A"]+freq["C"]+freq["G"]+freq["T"]+other) \
			 		      " BP; "\
			 		      freq["A"]" A; "\
			 		      freq["C"]" C; "\
			 		      freq["G"]" G; "\
			 		      freq["T"]" T; "\
			 		      other" other;" \
			 }' ${QUERY}
	loginfo "Done."
	
	loginfo "Reformating sequences..."
		lines=$(wc -l ${QUERY} | awk '{print $1}')
		awk -v lines=$lines ' \
			! /^>/ { \
					seq=tolower($0); \
					gsub(" ","",seq); \
					printf("     ") ;\
					for (i=0; i < 6; i++) { \
						f=substr(seq,i * 10, 10); \
						pos+=length(f); \
						f = f  substr("          ",1,10-length(f)); \
						printf("%s ",f) \
					}; \
					if (NR==lines) \
					  {pos-=1}; \
					printf("   %6d\n",pos) \
			   }' ${QUERY}
	loginfo "Done."
	
	loginfo "Closing sequence part..."
		echo "//"
	loginfo "Done."

popTmpDir

