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

	loginfo "Annotating ITS and TSU..."

	RRNADB="${NUCRRNA_DATA_DIR}/plants/nuc_RRNA.hmm"
	
	if [[ ! "$1" =~ ^/ ]]; then
		QUERY="${CALL_DIR}/$1"
	else
		QUERY="$1"
	fi


    ITSx -p "${ITS_DATA_DIR}/ITSx_db/HMMs" -i "${QUERY}" -o "output.itsx"
    
    ITS1=( $(sed -E 's/.*ITS1: *([0-9]+)-([0-9]+).*/\1 \2/' "output.itsx.positions.txt" ) )
    ITS2=( $(sed -E 's/.*ITS2: *([0-9]+)-([0-9]+).*/\1 \2/' "output.itsx.positions.txt" ) )
    TSU=( $(sed -E 's/.*5\.8S: *([0-9]+)-([0-9]+).*/\1 \2/' "output.itsx.positions.txt" ) )


    
    if [[ ${#ITS1[@]}=="2" ]]  ; then
    	echo "FT   misc_RNA        ${ITS1[0]}..${ITS1[1]}"
		echo 'FT                   /gene="ITS1"'
		echo 'FT                   /note="internal transcribed spacer 1, ITS1"'
	fi
	
    if [[ ${#TSU[@]}=="2" ]]  ; then
		echo "FT   rRNA            ${TSU[0]}..${TSU[1]}"
		echo 'FT                   /gene="5.8S rRNA"'
		echo 'FT                   /product="5.8S ribosomal RNA"'
	fi
	
    if [[ ${#ITS2[@]}=="2" ]]  ; then
		echo "FT   misc_RNA        ${ITS2[0]}..${ITS2[1]}"
		echo 'FT                   /gene="ITS2"'
		echo 'FT                   /note="internal transcribed spacer 2, ITS2"'
	fi
    

	hmmsearch --max ${RRNADB} ${QUERY} | \
		$AwkCmd '/Query: / { 
		                profil=$2
		                match($3,"[0-9][0-9]*");
		                lprof=substr($3,RSTART,RLENGTH)} 
		     / [0-9][0-9]* ! / { 
		                print profil,lprof,$7,$8,$10,$11}' | \
		$AwkCmd '($3 <=5) && (($2-$4) <=5) { 
		                  full=1;$5=$5-$3+1;$6=$6+($2-$4)
						}  
		                { loc=$5".."$6 } 
		         ($1 ~ /_RC$/) { 
		                  loc="complement("loc")" 
						} 
		         (full==1) { 
					      match($1,"_..*S")
		                  rrna=substr($1,RSTART+1,RLENGTH-1)
		                  print "FT   rRNA            " loc
		                  print "FT                   /gene=\""rrna" rRNA\""
		                  print "FT                   /product=\""rrna" ribosomal RNA\""
		                  full=0
		                }'	

    loginfo "Done."

popTmpDir

exit 0

