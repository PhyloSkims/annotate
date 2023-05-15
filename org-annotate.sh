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
cdsdetection_pass3="yes"
trnadetection="yes"
rrnadetection="yes"
idprefix="no"
tagprefix="no"
locusshift=1
organism="no"
country="no"
specimen="no"
project="no"
resetorganism="yes"
types="chloro"
partial=0
minlength=0
listfile="no"

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
	echo '  -L     | --locus-prefix <prefix>'
	echo '      Prefix used to build the locus tag of every annotated genes'
	echo '      generated locus tags follow the pattern : prefix_###,'
	echo '      where ### is a number following the order of gene in the embl file'
	echo '      starting at locus tag shift (default 1).'
	echo
	echo '  -S     | --locus-shift <###>'
	echo '      Start number for building locus tags'
	echo
	echo '  --list-file <FILENAME>'
	echo '      The chomosome list file name to file for ENA submmission'
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
	echo '      Do not annotate CDSs'
	echo
	echo '  -D     | --no-cds-pass1'
	echo '      Do not annotate core CDSs'
	echo
	echo '  -E     | --no-cds-pass2'
	echo '      Do not annotate rps12 CDS'
	echo
	echo '  -F     | --no-cds-pass3'
	echo '      Do not annotate shell and dust CDSs'
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

function over_junction() {
	local genome_length=$1

	$AwkCmd -v genome_length=$genome_length '
		function split_location(location,      p, status, newloc) {
			delete(state)
			delete(newloc)
			level = 0
			state[level] = ""
			subloc = location
			newloc[level] = ""
			locpart = ""
			start = 0
			pending = 0

			while(length(subloc)>0) {
				first_letter = substr(subloc,1,1)
				
				switch (first_letter) {
					case "j" :
						print "entereing in join" > /dev/stderr
						subloc = substr(subloc,6)
						level++
						state[level] = "join"
						break
					case "c" :
						print "entereing in complement" > /dev/stderr
						subloc = substr(subloc,12)
						level++
						state[level] = "complement"
						break
					case "," :
						print "next exon" > /dev/stderr
						subloc = substr(subloc,2)
						newloc[level] = newloc[level] ","
						break
					case ")" :
						print "Ending "state[level] > /dev/stderr
						subloc = substr(subloc,2)
						level--
						newloc[level] = newloc[level] state[level+1] "(" newloc[level+1] ")"
						if (pending == 1) {
							level--
							newloc[level] = newloc[level] state[level+1] "(" newloc[level+1] ")"
						}
						break
					default :
						print "Matching location: " subloc > /dev/stderr
						p = match(subloc,/[0-9]+\.\.[0-9]+/)
						piece = substr(subloc,RSTART,RLENGTH)
						subloc= substr(subloc,RLENGTH+1)
						p = match(piece,/[0-9]+/)                   
						from = substr(piece,RSTART,RLENGTH)
						piece = substr(piece,RLENGTH+1)
						p = match(piece,/[0-9]+/)                   
						to = substr(piece,RSTART,RLENGTH)
						if ((genome_length+0 >= from+0) &&
							(genome_length+0 <= to+0)) {
								status = "overlap"
								level++
								state[level] = "join"
								pending = 1
								piece = from ".." genome_length ",1.."(to - genome_length)
							} else {
								status = "out"
							}
						newloc[level] = newloc[level] piece
						print p, RSTART, from, to, status, piece > /dev/stderr
						break
				}
			}
			return newloc[0]
		}

		function cut_location(location) {
			pos = match(location,/,[^,]*$/)
			right = ""

			if (pos > 80) {
				while (pos > 80) {
					right = substr(location,pos+1) right
					location = substr(location,1,pos)
					pos = match(location,/,[^,]*,$/)
				}

				if (pos > 0) {
					right = substr(location,pos+1) right
					location = substr(location,1,pos)
				}
			}
			
			if (right !="") {
				right = "\nFT                   " right
				if (length(right) > 80) {
					right = cut_location(right)
				}
				location = location right
			} 
			return location
		}
		
		/^FT   [^ ]/ && (feature != "") {
			if (status != "remove") {
				if (status == "overlap") {
					print cut_location(fttype location)
					print "FT                   /note=\"" fttype " crosses the IR boundary\""
				} else {
					print cut_location(fttype location)
				}
				print feature
			}

			feature = ""
			fttype  = ""
			position= ""
		} 

		(feature != "") {
			feature = feature "\n" $0
		}

		#
		# Begining of a feature
		#
		/^FT   [^ ]/ {
			fttype=$2
			plocation=$3
			fttype=substr($0,1,21)
		} 

		#
		# Following of a location
		#
		/^FT    +[^ \/]/ && (plocation != "") {
			plocation = plocation $2
		}

		#
		# End of the location
		# And begining of the feature description
		# 
		/^FT +\// && (plocation != "") {
			feature = $0
			location = plocation
			split(location,parts,/\.\./)
			delete p
			pmin=1000000000
			pmax=0

			for (i in parts) {
				pos = match(parts[i],/[0-9]+/)
				if (pos+0 > 0) {
					j++
					p[i] = substr(parts[i], RSTART , RLENGTH)
					if (p[i]+0 > pmax+0) { pmax = p[i]}
					if (p[i]+0 < pmin+0) { pmin = p[i]}
				}
			}

			status = "ok"
			if (pmin+0 > genome_length+0) {
				status = "remove"
			} else {
				if (pmax+0 > genome_length+0) {
				status = "overlap"
				location = split_location(location)
				}
			}
			plocation=""
		}

		END {
			if (status != "remove") {
				if (status == "overlap") {
					print fttype location
					print "FT                   /note=\"CDS crosses the IR boundary\""
				} else {
					print fttype location
				}
				print feature
			}

		}
		' $2
}

function fastaIterator() {
	$AwkCmd '/^>/ {if (seq) printf("%s\f",seq); seq=""} \
	          {if (seq) seq=seq"\n"; seq=seq $1} \
	          END {print seq}' "$1"
}

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt -o s:t:o:b:P:i:fcrmhpl:NICDEFTRL:S: -l specimen:,ncbi-taxid:,organism:,country:,project:,id-prefix:,not-force-ncbi,chloroplast,nuclear-rdna,mitochondrion,partial,min-length:,help,no-normalization,no-ir-detection,no-cds,no-cds-pass1,no-cds-pass2,no-cds-pass3,no-trna,no-rrna,locus-prefix:,locus-shift:,list-file: -- "$@")
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
    -L|--locus-prefix) tagprefix="$2" ; shift ;;
    -S|--locus-shift) locusshift="$2" ; shift ;;
    -f|--not-force-ncbi)  resetorganism="no" ;;
    -p|--partial) partial="1" ;;
    -l|--min-length) minlength="$2" ; shift ;;
    -h|--help)  usage $0 0;;
    -N|--no-normalization)  normalization="no" ;;
    -I|--no-ir-detection)  irdetection="no" ;;
	-C|--no-cds) cdsdetection="no";;
	-D|--no-cds-pass1) cdsdetection_pass1="no";;
	-E|--no-cds-pass2) cdsdetection_pass2="no";;
	-F|--no-cds-pass3) cdsdetection_pass3="no";;
	-T|--no-trna) trnadetection="no";;
	-R|--no-rrna) rrnadetection="no";;
	--list-file) listfile="$2" ; shift;;
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    (*) break;;
    esac
    shift
done

loginfo "Locus tag prefix provided: $tagprefix"
loginfo "Locus tag numbered from..: $locusshift"
loginfo "NCBI taxid provided......: $taxid"

if [[ "$listfile" != "no" ]] ; then
	loginfo "contig info saved in.....: $listfile"
fi

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
loginfo "           pass 3........: $cdsdetection_pass3"
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
				# seqid=$(tr "." "_" <<< ${seqid})

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

						#
						# We are annotating a circular sequence
						#
						# 2kb of the beginnig of the sequence is added to its end
						# to allow for overlaping feature
						#
						add_circular_extension=0
						if (( partial == 0 )) ; then
							loginfo "Extends the sequence to allows for features overlaping juction"
							cp "${RESULTS}.norm.fasta" "${RESULTS}.norm.orig.fasta"
							cp "${RESULTS}.norm.fasta" "${RESULTS}.norm.frgs.fasta"
							cutseq "${RESULTS}.norm.orig.fasta" 1 2000 >> "${RESULTS}.norm.frgs.fasta"
							joinfasta "${RESULTS}.norm.frgs.fasta" > "${RESULTS}.norm.fasta"
							add_circular_extension=1
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
								cdsdetection_pass3=$cdsdetection_pass3 \
								${PROG_DIR}/detectors/cds/bin/go_cds.sh "${RESULTS}.norm.fasta" >> "${RESULTS}.annot"
							loginfo "Done."
						fi

						if (( add_circular_extension == 1 )) ; then
							mv "${RESULTS}.norm.orig.fasta" "${RESULTS}.norm.fasta"
							cp "${RESULTS}.annot" "${CALL_DIR}/${RESULTS}.annot.circular"
							over_junction $sl "${RESULTS}.annot" > "${RESULTS}.circular.annot"
							mv "${RESULTS}.circular.annot" "${RESULTS}.annot"
							add_circular_extension=0
						fi
						
						if (( partial == 0 )) ; then 
							topology="circular"
							defline="plastid, complete genome"
						else
							topology="linear"
							defline="plastid, partial sequence"
						fi

						if [[ "$listfile" != "no" ]] ; then
							if (( partial == 0 )) ; then
								echo "${seqid}	CHL	Circular-Chromosome	Chloroplast" >>  ${CALL_DIR}/$listfile
							else
								echo "${seqid}	CHL" >>  ${CALL_DIR}/$listfile
							fi
						fi

						notAnnoted  "${RESULTS}.annot" "${RESULTS}.norm.fasta" 100 > ${CALL_DIR}/not_annotated.fasta
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
								| $AwkCmd '{printf("%s;",$0)}' \
								| $AwkCmd 'BEGIN        {i=1} 
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
				
				loginfo "Unifying gene names"
					$AwkCmd '
					(FNR==NR) && /^FT                   \/gene="/ {
						gene = substr($0,29,100)
						gene = substr(gene,0,length(gene)-1)
						occurrence[gene]++
					}

					(FNR==1) && (FNR!=NR) {
						for(gene in occurrence){
						if (occurrence[gene]==1) {
							delete occurrence[gene]
						} else {
							occurrence[gene] = 1
						}
						}
					}

					(FNR!=NR) && /^FT                   \/gene="/ {
						gene = substr($0,29,100)
						gene = substr(gene,0,length(gene)-1)
						n = occurrence[gene]
						if (n > 0) {
							$0="FT                   /gene=\""gene"_"n"\""
							occurrence[gene]++
						}
					}

					(FNR!=NR) {
						print $0
					}
					' "${RESULTS}.sorted.annot" "${RESULTS}.sorted.annot" \
					  > "${RESULTS}.uniq_gene.annot"
				loginfo "Done."

				if [[ "$tagprefix" != "no" ]] ; then
				    loginfo "Adding locus tags from number: $locusshift..."
					cat "${RESULTS}.uniq_gene.annot" \
					| $AwkCmd -v tagprefix="$tagprefix" \
					          -v locusshift="$locusshift" '
						/^FT +\/locus_tag=""/ {
							sub(/locus_tag=""/,"locus_tag=\""tagprefix"_"locusshift"\"",$0);
							locusshift++;
						}
						{
							print $0
						}
					'
				    loginfo "Locus tags done."
				else
				    loginfo "Clearing locus tags done."
					egrep -v '^FT +\/locus_tag=""' \
					        "${RESULTS}.uniq_gene.annot"
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

