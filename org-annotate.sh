#!/bin/bash
#
#
#
#                           Annotate Organelle 
#
# The org-annotate pipeline aims to annotate fasta files produced by assembling 
# genome skimming. It has been developped in the context of the PhyloAlps 
# (http://phyloalps.org) and of the PhyloNorway (http://phylonorway.no) projects.
#
# Today it is able to produce EMBL flat files suitable for submission to ENA/EBI
# It provides annotation procedure for :
#
#     - Plant chloroplast genomes.
#     - Nuclear rDNA Region.
# 
#
#========================================================================================
#
# The template used for generating the EMBL files follows the recommendation presented
# at ENA documentation website (at the date of 2021-11-04).
# 
# https://ena-docs.readthedocs.io/en/latest/submit/fileprep/sequence-flatfile.html
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
cdsdetection="yes"
cdsdetection_pass1="yes"
cdsdetection_pass2="yes"
trnadetection="yes"
rrnadetection="yes"
idprefix="no"
organism="no"
country="no"
specimen="no"
project="no"
resetorganism="yes"
types="chloro"
partial=0
minlength=0

function usage {
	echo "Usage:" ;  
	echo "    $1 "'[-t|--ncbi-taxid ###] [-n|--no-normalization] \' 
	echo '       [-i|--no-ir-detection] [-h|--help] \ '
	echo '       [-o|--organism <organism_name>]  \ '
	echo '       [-c|--chloroplast|-r|--nuclear-rdna|-m|--mitochondrion] <FASTAFILE>' 
	echo
	echo "Options:"
	echo 
	echo ' Defining the sequence category'
	echo '  -c     | --chloroplast'
	echo '      Selects for the annotation of a chloroplast genome'
	echo '      This is the default mode'
	echo
	echo '  -r     | --nuclear-rdna'
	echo '      Selects for the annotation of the rDNA nuclear cistron'
	echo
#	echo '  -m     | --mitochondrion'
#	echo '      Selects for the annotation of an animal mitochondrion genome'
	echo
	echo ' Providing information about the sequence'
	echo '  -s     | --specimen ###'
	echo '      Represents the specimen voucher identifier '
	echo '      (e.g. for herbarium sample TROM_V_991090'
	echo '      will be added to the /specimen_voucher qualifier'
	echo
	echo '  -t     | --ncbi-taxid ###'
	echo '      Represents the ncbi taxid associated to the sequence'
	echo
	echo '  -o     | --organism <organism_name>'
	echo '      Allows for specifiying the organism name in the embl generated file'
	echo '      Spaces have to be substituted by underscore ex : Abies_alba'
	echo
	echo '  -b     | --country "<country_value>[:<region>][, <locality>]"'
	echo '      Location (at least country) of collection'
	echo
	echo '  -f     | --not-force-ncbi'
	echo '      Do not force the name of the organism to match the NCBI scientific name.'
	echo '      if the provided NCB taxid is publically defined'
	echo
	echo ' Information related to an ENA project'
	echo '  -P     | --project <ENA Project>'
	echo '      The accession number of the related project in ENA'
	echo
	echo '  -i     | --id-prefix <prefix>'
	echo '      prefix used to build the sequence identifier'
	echo '      the number of the contig is append to the prefix'
	echo '      to build the complete id. This id is used only if'
	echo '      an ENA project is specified.'
 	echo
	echo ' Annotation of partial sequences'
	echo '  -p     | --partial'
	echo '      Indicates that the genome sequence is partial and therefore in several contigs'
	echo
	echo '  -l     | --min-length'
	echo '      Indicates for partial mode the minimum length of contig to annotate'
	echo
	echo ' Setting up the sequence annotation'
	echo '  -N     | --no-normalization'
	echo '      Does not normalize the sequence befire annotation'
	echo
	echo '  -I     | --no-ir-detection'
	echo '      Does not look for inverted repeats in the plastid genome'
	echo
	echo '  -C     | --no-cds'
	echo '      Do not annotate CDS'
	echo
	echo '  -D     | --no-cds-pass1'
	echo '      Do not annotate core CDS'
	echo
	echo '  -E     | --no-cds-pass2'
	echo '      Do not annotate rps12 CDS'
	echo
	echo '  -T     | --no-trna'
	echo '      Do not look for transfert RNA'
	echo
	echo '  -R     | --no-rrna'
	echo '      Do not look for ribosomal RNA'
   	exit $2
}


function ncbiscientificname {
	local efetch='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
	local params='db=taxonomy&id='$1'&mode=text&report=xml'
	local url="${efetch}?${params}" 
	curl $url \
	| grep '<ScientificName>' \
	| sed 's/<ScientificName>//' \
	| sed 's@</ScientificName>@@' \
	| sed 's@^ *@@' \
	| sed 's@ *$@@' \
	| head -1
}

function split80 {
	local text=$1
	local pattern
	local header

	if (( $# >= 2 )) ; then
		header=$2
	else
		header=""
	fi

	if (( $# >= 3 )) ; then
		pattern=$3
	else
		pattern=" "
	fi

	echo $text \
	| $AwkCmd -v pattern="$pattern" -v header="$header" '
		BEGIN { line = header }
		{
			n = split($0,parts,pattern)
			j = 1
			for (i = 1; i <= n; i++) {
				if (length(line) + length(parts[i]) > 79) {
					print line
					line = header
					j = i
				} 
				if (i > j) line = line pattern
				line = line parts[i]
			}
			$0 = line    
        }		
		END {
			print $0
		}
	  '
}

function fastaIterator() {
	$AwkCmd '/^>/ {if (seq) printf("%s\f",seq); seq=""} \
	          {if (seq) seq=seq"\n"; seq=seq $1} \
	          END {print seq}' "$1"
}

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt -o s:t:o:b:P:i:fcrmhpl:NICDETR -l specimen:,ncbi-taxid:,organism:,country:,project:,id-prefix:,not-force-ncbi,chloroplast,nuclear-rdna,mitochondrion,partial,min-length:,help,no-normalization,no-ir-detection,no-cds,no-cds-pass1,no-cds-pass2,no-trna,no-rrna -- "$@")
then
    # something went wrong, getopt will put out an error message for us
    usage $0 1
fi

eval set -- "$options"

while [ $# -gt 0 ]
do
    case $1 in
    -c|--chloroplast) types="chloro" ;;
    -r|--nuclear-rdna) types="nucrdna" ;;
#    -m|--mitochondrion) types="mito" ;;
    -t|--ncbi-taxid) taxid="$2" ; shift ;;
	-s|--specimen) specimen="$2" ; shift ;;
	-b|--country) country="$2" ; shift ;;
    -o|--organism) organism="$2" ; shift ;;
    -P|--project) project="$2" ; shift ;;
    -i|--id-prefix) idprefix="$2" ; shift ;;
    -f|--not-force-ncbi)  resetorganism="no" ;;
    -p|--partial) partial="1" ;;
    -l|--min-length) minlength="$2" ; shift ;;
    -h|--help)  usage $0 0;;
    -N|--no-normalization)  normalization="no" ;;
    -I|--no-ir-detection)  irdetection="no" ;;
	-C|--no-cds) cdsdetection="no";;
	-D|--no-cds-pass1) cdsdetection_pass1="no";;
	-E|--no-cds-pass2) cdsdetection_pass2="no";;
	-T|--no-trna) trnadetection="no";;
	-R|--no-rrna) rrnadetection="no";;
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    (*) break;;
    esac
    shift
done

loginfo "NCBI taxid provided......: $taxid"

if [[ "$taxid" != "no" ]] ; then
	scientificname=$(ncbiscientificname $taxid)
	loginfo "NCBI scientific name.....: $scientificname"	
	if [[ -z "$scientificname" ]] ; then
		loginfo "    Unknown taxid."
	else
		loginfo "Organism name from taxid.: $resetorganism"
		if [[ "$resetorganism" == "yes" ]] ; then
			organism=$(echo $scientificname | tr ' ' '_')
		fi
	fi
fi

loginfo "Annotating mode..........: $types"
loginfo "Sequence normalization...: $normalization"
loginfo "IR detection mode........: $irdetection"
loginfo "CDS detection mode.......: $cdsdetection"
loginfo "           pass 1........: $cdsdetection_pass1"
loginfo "           pass 2........: $cdsdetection_pass2"
loginfo "tRNA detection mode......: $trnadetection"
loginfo "rRNA detection mode......: $rrnadetection"
loginfo "Organism.................: $organism"
loginfo "Country..................: $country"
loginfo "Partial mode.............: $partial"
loginfo "Minimum length...........: $minlength"


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

	IFS=$'\f'

	sequence_number=0
	
	for sequence in $(fastaIterator "${QUERY}") ; do
	    let sequence_number=sequence_number+1
		unset IFS
		if [[ ! -z "${sequence}" ]] ; then
			echo "${sequence}" > toannotate.fasta
			sl=$(seqlength "toannotate.fasta")
			
			if (( sl >= minlength )) ; then
			
				seqid=$($AwkCmd '(NR==1) {print substr($1,2,1000)}' toannotate.fasta)
				
				case "$types" in 
					chloro) 
						loginfo "Annotating a plant chloroplast genome..."
						
						if [[ "$irdetection" == "yes" ]] && (( partial == 0 )) ; then
					
							if [[ "$normalization" == "yes" ]] ; then
								loginfo "Normalizing the structure of the Chloroplast sequence..."
								loginfo "   LSC + IRB + SSC + IRA"
								${PROG_DIR}/detectors/normalize/bin/go_normalize.sh toannotate.fasta > "${RESULTS}.norm.fasta"
								loginfo "Done."
							else
								loginfo "No normalization the structure of the Chloroplast sequence..."
								cat toannotate.fasta > "${RESULTS}.norm.fasta"
							fi
							
							loginfo "Annotating the Inverted repeats and Single copies (LSC and SSC)..."
								${PROG_DIR}/detectors/ir/bin/go_ir.sh "${RESULTS}.norm.fasta" > "${RESULTS}.annot"		
							loginfo "Done."
							
						else
							loginfo "No normalization of the structure of the Chloroplast sequence..."
							cat toannotate.fasta > "${RESULTS}.norm.fasta"
							rm -f "${RESULTS}.annot"
							touch "${RESULTS}.annot"
						fi
						
						if [[ "$trnadetection" == "yes" ]] ; then
							loginfo "Annotating the tRNA..."
							${PROG_DIR}/detectors/trna/bin/go_trna.sh "${RESULTS}.norm.fasta" >> "${RESULTS}.annot"
							loginfo "Done."
						fi
						
						if [[ "$rrnadetection" == "yes" ]] ; then
							loginfo "Annotating the rRNA genes..."
							${PROG_DIR}/detectors/rrna/bin/go_rrna.sh "${RESULTS}.norm.fasta" >> "${RESULTS}.annot"
							loginfo "Done."
						fi
					
						if [[ "$cdsdetection" == "yes" ]] ; then
							loginfo "Annotating the CDS..."
							cdsdetection_pass1=$cdsdetection_pass1 \
								cdsdetection_pass2=$cdsdetection_pass2 \
								${PROG_DIR}/detectors/cds/bin/go_cds.sh "${RESULTS}.norm.fasta" >> "${RESULTS}.annot"
							loginfo "Done."
						fi
						
						if (( partial == 0 )) ; then 
							topology="circular"
							defline="plastid, complete genome"
						else
							topology="linear"
							defline="plastid, partial sequence"
						fi
						;;
						
					nucrdna) 
						loginfo "Annotating a plant rDNA cistron..."
						
						if [[ "$normalization" == "yes" ]] ; then
							loginfo "Normalizing the structure of the cistron sequence..."
							${PROG_DIR}/detectors/normalizerdna/bin/go_normalizerdna.sh toannotate.fasta > "${RESULTS}.norm.fasta"								
							loginfo "Done."
						else
							loginfo "No normalization of the structure of the cistron sequence..."
							cat toannotate.fasta > "${RESULTS}.norm.fasta"
						fi
						
						loginfo "Annotating the rRNA genes..."
							${PROG_DIR}/detectors/nucrrna/bin/go_nucrrna.sh "${RESULTS}.norm.fasta" > "${RESULTS}.annot"
						loginfo "Done."
			
						topology="linear"
						defline=$(cat "${RESULTS}.annot" \
						        | grep "/gene=" \
								| sed -E 's@^FT */gene="([^"]*)".*$@\1@' \
								| sort -u \
								| awk '{printf("%s;",$0)}' \
								| awk 'BEGIN        {i=1} 
								       /18S rRNA;/  {gene[i]="18S rRNA"; i++} 
									   /ITS1;/      {gene[i]="ITS1"; i++} 
									   /5.8S rRNA;/ {gene[i]="5.8S rRNA"; i++} 
									   /ITS2;/      {gene[i]="ITS2"; i++} 
									   /28S rRNA;/  {gene[i]="28S rRNA"; ii++} 
									   END          {
										for (i=1;i <= length(gene); i++) {
											separator =""
											if (i  < (length(gene)-1)) separator=", " 
											if (i == (length(gene)-1)) separator=" and " 
											printf("%s gene%s",gene[i],separator)
											}
										}') 
						# defline="18S rRNA gene, ITS1, 5.8S rRNA gene, ITS2 and 28S rRNA gene"
						;;
						
					mito) 
						loginfo "Annotating an animal mitochondrial genome..."
						logerror "Not yet implemented"
			
						if (( partial == 0 )) ; then 
							topology="circular"
							defline="mitochondrion, complete genome"
						else
							topology="linear"
							defline="mitochondrion, partial sequence"
						fi
						
						exit 1
						;;
						
					*) 
						usage $0 1;;
				esac
									
				if [[ "${organism}" == "no" ]]; then
					organism="{organism}"
				else
					organism="$(echo ${organism} | tr '_' ' ')"
				fi

				if [[ "$specimen" != "no" ]] ; then
					defline="${defline}, voucher ${specimen}"
				fi
				
				sl=$(seqlength "${RESULTS}.norm.fasta")
				
				loginfo "Printing minimal header..."		
					echo "ID   XXX; XXX; ${topology}; genomic DNA; XXX; XXX; ${sl} BP."
					echo "XX"
					echo "AC   XXX;"
					echo "XX"
					if [[ "${project}" != "no" ]] ; then 
					sequence_number
						if [[ "$idprefix" != "no" ]] ; then
							seqid="${idprefix}${sequence_number}"
						fi
						echo "AC * _${seqid};"
						echo "XX"
						echo "PR   Project:${project};"              
						echo "XX"
						echo "DE   XXX"
					else
						split80 "${organism} ${defline}." "DE   "
					fi	
					
					echo "XX"
				loginfo "Done."
			
				loginfo "Printing annotations header..."
			    	echo "FH   Key             Location/Qualifiers"
				loginfo "Done."
	
				loginfo "Printing the source feature"
						echo "FT   source          1..${sl}"                               
	
					if [[ "${organism}" != "no" ]] ; then 
						echo "FT                   /organism=\"${organism}\""              
					fi	
					
					case "${types}" in 
						chloro)  
							echo "FT                   /organelle=\"plastid:chloroplast\"" 
						;;
						mito)    
							echo "FT                   /organelle=\"mitochondrion\""       
						;;
						*) 
							loginfo "Nuclear sequence"  
						;;
					esac
					
						echo "FT                   /mol_type=\"genomic DNA\""              
					
					if [[ "${specimen}" != "no" ]] ; then 
						echo "FT                   /specimen_voucher=\"${specimen}\""            
					fi
					
					if [[ "${taxid}" != "no" ]] ; then 
						echo "FT                   /db_xref=\"taxon:${taxid}\""            
					fi
					
					if [[ "${country}" != "no" ]] ; then 
						echo "FT                   /country=\"${country}\""            
					fi
					
				loginfo "Done."
				
				loginfo "Ordering annotations..."
					$AwkCmd '(entry && /^.....(misc|repeat|rRNA|tRNA|CDS|source)/) { \
	                           print pos,entry } \
						 /^.....(misc|repeat|rRNA|tRNA|CDS|source)/ { \
					        match($3,"[0-9][0-9]*"); \
					        pos=substr($3,RSTART,RLENGTH)*1000 + 1; \
					        entry=$0;    \
					        next} \
					      { entry=entry "@" $0} \
	 					END {print pos,entry}' "${RESULTS}.annot" | \
					sort -nk1 |\
					$AwkCmd '{ \
					        match($0,"^[0-9]* ");\
					        line=substr($0,RLENGTH+1);\
							gsub("@","\n",line); \
					        print line}' > "${RESULTS}.sorted.annot"
				loginfo "Done."
				
				if [[ "$idprefix" != "no" ]] ; then
				    loginfo "Adding locus tags..."
					cat "${RESULTS}.sorted.annot" \
					| $AwkCmd -v idprefix="$idprefix" '
					    BEGIN {n=1}
						/^FT +\/locus_tag=""/ {
							sub(/locus_tag=""/,"locus_tag=\""idprefix"_"n"\"",$0);
							n++;
						}
						{
							print $0
						}
					'
				    loginfo "Locus tags done."
				else
				    loginfo "Clearing locus tags done."
					egrep -v '^FT +\/locus_tag=""' \
					        "${RESULTS}.sorted.annot"
				    loginfo "Clearing of tags done."
				fi
				
				loginfo "Closing annotations table..."
					echo "XX"
				loginfo "Done."
				
				loginfo "Computing statistics on nucleotide usage..."
					$AwkCmd '! /^>/ { \
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
					lines=$(wc -l "${RESULTS}.norm.fasta" | $AwkCmd '{print $1}')
					
					loginfo "Sequence length $(seqlength ${RESULTS}.norm.fasta)"
					loginfo "lines $lines"
					formatfasta "${RESULTS}.norm.fasta" | \
					$AwkCmd -v lines=$lines ' \
						! /^>/ { \
								seq=tolower($0); \
								gsub(" ","",seq); \
								printf("     ") ;\
								for (i=0; i < 6; i++) { \
									f=substr(seq,i * 10 + 1, 10); \
									pos+=length(f); \
									f = f  substr("          ",1,10-length(f)); \
									printf("%s ",f) \
								}; \
								printf("   %6d\n",pos) \
						   }'
				loginfo "Done."
				
				loginfo "Closing sequence part..."
					echo "//"
				loginfo "Done."
			fi # End of the minimum length condition
		fi  # End of not empty sequence condition
		
		IFS=$'\f'
	done # End of the loop over the sequences
popTmpDir
loginfo "Annotation done."

