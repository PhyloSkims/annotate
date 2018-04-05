#!/bin/bash
#
#                           BUILD REFERENCE THE CAU TRNA LIBRARy
#
#========================================================================================

# -- CAUTION -- Works as long than the script 
#               is not called through a symlink
THIS_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${THIS_DIR}/../../../scripts/bash_init.sh"

function fasta1li {

    $AwkCmd '/^>/ {if (sequence) \
                  {print sequence}; \
                   print $0; \
                   sequence=""} \
            !/^>/ {sequence = sequence $0} \
            END   {print sequence}' $1
}


function filtertrna() {

	$AwkCmd -F '_' 'BEGIN       {RS=">"}                                      \
            (! /^$/)    {trna=$1;                                     \
                         ac=$2"_"$3;}                                 \
            (ac!=oldac &&                                             \
             trnas["trnfM"]==1 &&                                     \
             trnas["trnM"]==1 &&                                      \
             trnas["trnI"]==1                                         \
            )           {print seqs}                                  \
            (ac!=oldac) {trnas["trnfM"]=0;                            \
                         trnas["trnM"]=0;                             \
                         trnas["trnI"]=0;                             \
                         seqs="";                                     \
                         oldac=ac                                     \
                         }                                            \
            (! /^$/)    {seqs=seqs"\n>"$0;                            \
            			 trnas[trna]=1;}                              \
            END         {if (trnas["trnfM"]==1 &&                     \
                             trnas["trnM"]==1 &&                      \
                             trnas["trnI"]==1)                        \
                             print seqs}' $1 |                        \
    egrep -v "^ *$"
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
			grep -v -- --  | \
			filtertrna > ${CAUFILE}
	loginfo "Done"
		
	loginfo "Installing the CAU tRNA database..."
		
		cp ${CAUFILE} "${TRNA_DATA_DIR}/CAU_tRNA_DB.fasta"

	loginfo "Done"

popTmpDir