#!/bin/bash
#
#                           BUILD REFERENCE THE CAU TRNA LIBRARy
#
#========================================================================================

# -- CAUTION -- Works as long than the script 
#               is not called through a symlink
SCRIPT_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${SCRIPT_DIR}/../../../scripts/bash_init.sh"

function fasta1li {

    awk '/^>/ {if (sequence) \
                  {print sequence}; \
               print $0; \
               sequence=""} \
         !/^>/ {sequence = sequence $0} \
         END {print sequence}' $1
}

function dereplicate {
	DATA=$1
	sumaclust -t 1 $DATA | \
		fasta1li | \
		grep -A 1 '^>' | \
		grep -A1 'cluster_center=True;' | \
		grep -v -- -- | \
		sed -E "s/count=[0-9]+; //" | \
		sed 's/cluster_weight/count/' | \
		awk ' /^>/ {SEQ++;$1=$1"_"SEQ;print $0} \
			 !/^>/ {print $0}'
}

function extractSeqs {
	rm -f $$.index
	fastaindex -f $1 \
	           -i $$.index
	
	for id in `cat $2`; do
		fastafetch -f $1 \
				   -i $$.index \
				   -q $id 
	done
	
	rm -f $$.index
}

function goodtrna {
    local QUERY=$1
    local REF=$2
	sumatra -t 0.90 -x $QUERY $REF | \
		sed -E 's/.(trn.M?)[_A-Z0-9]+/ \1 /' | \
		sort -k 1,2 | \
		awk '(OLD) && ($1!=OLD) {print OLD,c["trnM"],c["trnfM"],c["trnI"]} \
		     (OLD !=$1)         {c["trnM"]=0;c["trnfM"]=0;c["trnI"]=0;OLD=$1} \
		                        {c[$2]+=$5}' | awk '{p=0;} \
		     ($2 > $3) && ($2 > $4) { print $0,"trnM";p=1 } \
		     ($3 > $2) && ($3 > $4) {print $0,"trnfM";p=1} \
		     ($4 > $2) && ($4 > $3) {print $0,"trnI";p=1} \
		     (p==0) {print $0,"----"}' |      sed 's/_/ /' | \
		awk '{print $1"_"$2,$3,$4,$5,$1,$6}' | \
		awk '(($2+$3+$4) > 1) && ($5==$6) {print $1}'
}

pushTmpDir ORG.buildSCDB

	CAUFILE=CAU.fasta	

	openLogFile "${TRNA_DATA_DIR}/CAU_tRNA_DB.log"
	
	loginfo "Selecting Viridiplantae genebank entries..."
		VIRIDIPLANTAE=$(${PROG_DIR}/../../normalize/tools/selectViridiplantae.sh $*) 
		loginfo " --> $(echo ${VIRIDIPLANTAE} | wc -w) entries selected"
	loginfo "Done"
	
	loginfo "Extracting the CAU tRNA from the plants entries..."
		${PROG_DIR}/extract_refCAUtRNA.sh ${VIRIDIPLANTAE} | \
			fasta1li | \
			egrep -A 1 '^>trn(I|M|fM)' | \
			grep -v -- -- > ${CAUFILE}
	loginfo "Done"

	loginfo "Sorting the CAU tRNA..."
	
		loginfo "Extract and dereplicate trnI..."
		egrep -A 1 '^>trnI_'  ${CAUFILE} | grep -v -- -- > trnI.fasta
		dereplicate trnI.fasta > trnCAU.fasta
		
		loginfo "Extract and dereplicate trnM..."
		egrep -A 1 '^>trnM_'  ${CAUFILE} | grep -v -- -- > trnM.fasta
		dereplicate trnM.fasta >> trnCAU.fasta
		
		loginfo "Extract and dereplicate trnfM..."
		egrep -A 1 '^>trnfM_' ${CAUFILE} | grep -v -- -- > trnfM.fasta
		dereplicate trnfM.fasta >> trnCAU.fasta
		
	loginfo "Done"
	
	loginfo "Cleaning the CAU tRNA..."

		loginfo "First pass..."
		
		goodtrna 	trnCAU.fasta trnCAU.fasta 	 > trnCAU.good.ids
		extractSeqs trnCAU.fasta trnCAU.good.ids > trnCAU.good.fasta

		loginfo " --> $(wc -l trnCAU.good.ids) sequences kept"
		
		loginfo "Second pass..."
		
		goodtrna 	trnCAU.fasta trnCAU.good.fasta > trnCAU.good.ids
		extractSeqs trnCAU.fasta trnCAU.good.ids   > trnCAU.good.fasta	
		
		loginfo " --> $(wc -l trnCAU.good.ids) sequences kept"

	loginfo "Done"
	
	loginfo "Installing the CAU tRNA database..."
		
		cp trnCAU.good.fasta "${TRNA_DATA_DIR}/CAU_tRNA_DB.fasta"

	loginfo "Done"

popTmpDir