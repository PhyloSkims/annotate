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

#
# Management of options
#

taxid="no"
normalization="yes"
irdetection="yes"

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt -o t:ih -l ncbi-taxid:,no-ir-detection,help -- "$@")
then
    # something went wrong, getopt will put out an error message for us
    exit 1
fi

eval set -- "$options"

while [ $# -gt 0 ]
do
    case $1 in
    -t|--ncbi-taxid) taxid="$2" ; shift;;
    -i|--no-ir-detection)  irdetection="no" ;;
    -h|--help)  echo "Usage:" ;  
    			echo "    $0 "'[-t|--ncbi-taxid ###] [-n|--no-normalization] \' 
    			echo "       [-i|--no-ir-detection] [-h|--help] <FASTAFILE>"
    			echo
    			echo "Options:"
    			echo '  -t ### | --ncbi-taxid ###'
    			echo '      ### represents the ncbi taxid associated to the sequence'
    			echo
    			echo '  -i     | --no-ir-detection'
    			echo '      Does not look for inverted repeats in the plastid genome'
               	exit 0;;
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    (*) break;;
    esac
    shift
done

#############################

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
		
	if [ "$irdetection"=="yes" ]; then

		loginfo "Normalizing the structure of the Chloroplast sequence..."
			loginfo "   LSC + IRB + SSC + IRA"
			${PROG_DIR}/detectors/normalize/bin/go_normalize.sh ${QUERY} > "${RESULTS}.norm.fasta"
		loginfo "Done."
		
		loginfo "Annotating the Inverted repeats and Single copies (LSC and SSC)..."
			${PROG_DIR}/detectors/ir/bin/go_ir.sh "${RESULTS}.norm.fasta" > "${RESULTS}.annot"		
		loginfo "Done."
		
	fi
	
	loginfo "Annotating the tRNA..."
		${PROG_DIR}/detectors/trna/bin/go_trna.sh "${RESULTS}.norm.fasta" >> "${RESULTS}.annot"
	loginfo "Done."
	
	loginfo "Annotating the rRNA genes..."
		${PROG_DIR}/detectors/rrna/bin/go_rrna.sh "${RESULTS}.norm.fasta" >> "${RESULTS}.annot"
	loginfo "Done."

	loginfo "Annotating the CDS..."
		${PROG_DIR}/detectors/cds/bin/go_cds.sh "${RESULTS}.norm.fasta" >> "${RESULTS}.annot"
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
			 }' "${RESULTS}.norm.fasta"
	loginfo "Done."
	
	loginfo "Reformating sequences..."
		lines=$(wc -l "${RESULTS}.norm.fasta" | awk '{print $1}')
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
			   }' "${RESULTS}.norm.fasta"
	loginfo "Done."
	
	loginfo "Closing sequence part..."
		echo "//"
	loginfo "Done."

popTmpDir

