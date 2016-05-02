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
organism="no"
types="chloro"

function usage {
	echo "Usage:" ;  
	echo "    $1 "'[-t|--ncbi-taxid ###] [-n|--no-normalization] \' 
	echo '       [-i|--no-ir-detection] [-h|--help] \ '
	echo '       [-o|--organism <organism_name>]  \ '
	echo '       [-c|--chloroplast|-r|--nuclear-rdna|-m|--mitochondrion] <FASTAFILE>' 
	echo
	echo "Options:"
	echo '  -t ### | --ncbi-taxid ###'
	echo '      ### represents the ncbi taxid associated to the sequence'
	echo
	echo '  -i     | --no-ir-detection'
	echo '      Does not look for inverted repeats in the plastid genome'
	echo
	echo '  -o     | --organism <organism_name>'
	echo '      Allows for specifiying the organism name in the embl generated file'
	echo '      Spaces have to be substituted by underscore ex : Abies_alba'
	echo
	echo '  -c     | --chloroplast'
	echo '      Selects for the annotation of a chloroplast genome'
	echo '      This is the default mode'
	echo
	echo '  -r     | --nuclear-rdna'
	echo '      Selects for the annotation of the rDNA nuclear cistron'
	echo
	echo '  -m     | --mitochondrion'
	echo '      Selects for the annotation of an animal mitochondrion genome'
   	exit $2
}


# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt -o t:o:icrmh -l ncbi-taxid:,organism,no-ir-detection,chloroplast,nuclear-rdna,mitochondrion,help -- "$@")
then
    # something went wrong, getopt will put out an error message for us
    exit 1
fi

eval set -- "$options"

while [ $# -gt 0 ]
do
    case $1 in
    -t|--ncbi-taxid) taxid="$2" ; shift ;;
    -i|--no-ir-detection)  irdetection="no" ;;
    -o|--organism) organism="$2" ; shift ;;
    -c|--chloroplast) types="chloro" ;;
    -r|--nuclear-rdna) types="nucrdna" ;;
    -m|--mitochondrion) types="mito" ;;
    -h|--help)  usage $0 0;;
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    (*) break;;
    esac
    shift
done

echo $type
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

	case "$types" in 
		chloro) 
			loginfo "Annotating a plant chloroplast genome..."
	
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
				tcsh -f ${PROG_DIR}/detectors/cds/bin/go_cds.sh "${RESULTS}.norm.fasta" >> "${RESULTS}.annot"
			loginfo "Done."
			
			topology="circular"
			defline="plastid, complete genome"
			;;
		nucrdna) 
			loginfo "Annotating a plant rDNA cistron..."
			
			loginfo "Normalizing the structure of the cistron sequence..."
				${PROG_DIR}/detectors/normalizerdna/bin/go_normalizerdna.sh ${QUERY} > "${RESULTS}.norm.fasta"
			loginfo "Done."
			
			loginfo "Annotating the rRNA genes..."
				${PROG_DIR}/detectors/nucrrna/bin/go_nucrrna.sh "${RESULTS}.norm.fasta" > "${RESULTS}.annot"
			loginfo "Done."

			topology="linear"
			defline="18S rRNA gene, ITS1, 5.8S rRNA gene, ITS2 and 28S rRNA gene"
			;;
		mito) 
			loginfo "Annotating an animal mitochondrial genome..."
			logerror "Not yet implemented"

			topology="circular"
			defline="mitochondrion, complete genome"

			exit 1
			;;
		*) 
			echo usage $0 1;;
	esac
					
	if [[ "${organism}" == "no" ]]; then
		organism="{organism}"
	else
		organism="$(echo ${organism} | tr '_' ' ')"
	fi
	
	loginfo "Printing minimal header..."		
		echo "ID   XXX; XXX; ${topology}; genomic DNA; XXX; XXX; $(seqlength ${RESULTS}.norm.fasta) BP."
		echo "XX"
		echo "AC   XXX;"
		echo "DE   ${organism} ${defline}."
		echo "XX"
	loginfo "Done."

	loginfo "Printing annotations header..."
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

