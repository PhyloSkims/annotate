#!/bin/bash
#
#                           BUILD REFERENCE SINGLE COPY LIBRARIES
#
#========================================================================================

# -- CAUTION -- Works as long than the script 
#               is not called through a symlink
THIS_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${THIS_DIR}/../../../scripts/bash_init.sh"

pushTmpDir ORG.buildSCDB
	
	openLogFile "${IR_DATA_DIR}/SingleCopyDB.log"
	
	loginfo "Selecting Viridiplantae genebank entries..."
		VIRIDIPLANTAE=$(${PROG_DIR}/selectViridiplantae.sh $*) 
		loginfo " --> $(echo ${VIRIDIPLANTAE} | wc -w) entries selected"
	loginfo "Done"
	
	
	#
	# Deals with Long Single Copies (LSC)
	#
	
	loginfo "Extracting Long Single Copies (LSC)..."
		${PROG_DIR}/extract_refLSC.sh ${VIRIDIPLANTAE} > LSC.fasta
		loginfo " --> $(fastaCount LSC.fasta) retreived sequences"
	loginfo "Done"
	
	
	
	loginfo "Building LSC coorientation graph..."
		${PROG_DIR}/coorienteSC.sh LSC.fasta 20000 ${ORG_LOGFILE} > LSC.tgf
		${PROG_DIR}/cc.py LSC.tgf > LSC.cc
		loginfo " --> $(awk '{print $1}' LSC.cc | uniq | wc -l) connected componants"
	loginfo "Done"


	loginfo "Indexing LCS..."
		fastaindex -f LSC.fasta -i LSC.index
	loginfo "Done"
	
	
	
	loginfo "Extracting main connected components for LCS..."
		rm -f LSC.direct.fasta
		touch LSC.direct.fasta
		for id in `awk '($1==0) {print $2}' LSC.cc`; do
			fastafetch -f LSC.fasta -i LSC.index -q "${id}" >> LSC.direct.fasta
		done
		loginfo " --> $(fastaCount LSC.direct.fasta) sequences"
	loginfo "Done"
	
	
	
	loginfo "Extracting second connected components for LCS..."
		rm -f LSC.reverse.fasta
		touch LSC.reverse.fasta
		for id in `awk '($1==1) {print $2}' LSC.cc`; do
			fastafetch -f LSC.fasta -i LSC.index -q "${id}" >> LSC.reverse.fasta
		done		
		loginfo " --> $(fastaCount LSC.reverse.fasta) sequences"
	loginfo "Done"
	
	
	
	loginfo "merging both connected components for LCS..."
		fastarevcomp LSC.reverse.fasta >> LSC.direct.fasta
		loginfo " --> $(fastaCount LSC.direct.fasta) sequences in total"
	loginfo "Done"



	loginfo "Checking LCS homogeneity..."
		${PROG_DIR}/coorienteSC.sh LSC.direct.fasta 20000 ${ORG_LOGFILE} > LSC_RefDB.tgf		
		${PROG_DIR}/cc.py LSC_RefDB.tgf > LSC_RefDB.cc
		NCC=$(awk '{print $1}' LSC_RefDB.cc | uniq | wc -l)
		if (( $NCC == 1 )); then
			loginfo " --> $NCC connected componants"
		else
			logwarning " --> $NCC connected componants"
		fi
	loginfo "Done"
	
	
	
	loginfo "Installing LCS reference databases..."
		cp LSC.direct.fasta "${IR_DATA_DIR}/LSC_RefDB.fasta"
		cp LSC_RefDB.tgf    "${IR_DATA_DIR}/LSC_RefDB.tgf"
		
	loginfo "Done"
		
	#
	# Deals with Short Single Copies (SSC)
	#
	
	loginfo "Extracting Short Single Copies (SSC)..."
		${PROG_DIR}/extract_refSSC.sh ${VIRIDIPLANTAE} > SSC.fasta
		loginfo " --> $(fastaCount SSC.fasta) retreived sequences"
	loginfo "Done"


	
	loginfo "Building SSC coorientation graph..."
		${PROG_DIR}/coorienteSC.sh SSC.fasta 5000 ${ORG_LOGFILE} > SSC.tgf
		${PROG_DIR}/cc.py SSC.tgf > SSC.cc
		loginfo " --> $(awk '{print $1}' SSC.cc | uniq | wc -l) connected componants"
	loginfo "Done"



	loginfo "Indexing SSC..."
		fastaindex -f SSC.fasta -i SSC.index
	loginfo "Done"
	
	
	
	loginfo "Extracting main connected components for SSC..."
		rm -f SSC.direct.fasta
		touch SSC.direct.fasta
		for id in `awk '($1==0) {print $2}' SSC.cc`; do
			fastafetch -f SSC.fasta -i SSC.index -q "${id}" >> SSC.direct.fasta
		done
		loginfo " --> $(fastaCount SSC.direct.fasta) sequences"
	loginfo "Done"
	
	
	
	loginfo "Extracting second connected components for SSC..."
		rm -f SSC.reverse.fasta
		touch SSC.reverse.fasta
		for id in `awk '($1==1) {print $2}' SSC.cc`; do
			fastafetch -f SSC.fasta -i SSC.index -q "${id}" >> SSC.reverse.fasta
		done		
		loginfo " --> $(fastaCount SSC.reverse.fasta) sequences"
	loginfo "Done"
	
	
	
	loginfo "merging both connected components for SSC..."
		fastarevcomp SSC.reverse.fasta >> SSC.direct.fasta
		loginfo " --> $(fastaCount SSC.direct.fasta) sequences in total"
	loginfo "Done"



	loginfo "Checking SSC homogeneity..."
		${PROG_DIR}/coorienteSC.sh SSC.direct.fasta 5000 ${ORG_LOGFILE} > SSC_RefDB.tgf		
		${PROG_DIR}/cc.py SSC_RefDB.tgf > SSC_RefDB.cc
		NCC=$(awk '{print $1}' SSC_RefDB.cc | uniq | wc -l)
		if (( $NCC == 1 )); then
			loginfo " --> $NCC connected componants"
		else
			logwarning " --> $NCC connected componants"
		fi
	loginfo "Done"
	
	
	
	loginfo "Installing SSC reference databases..."
		cp SSC.direct.fasta "${IR_DATA_DIR}/SSC_RefDB.fasta"
		cp SSC_RefDB.tgf    "${IR_DATA_DIR}/SSC_RefDB.tgf"		
	loginfo "Done"
		
	loginfo "Installing blast version of the SC_RefDB reference databases..."
		cat "${IR_DATA_DIR}/LSC_RefDB.fasta" \
		    "${IR_DATA_DIR}/SSC_RefDB.fasta" > SC_RefDB.fasta
		    
		makeblastdb -in SC_RefDB.fasta \
		            -dbtype nucl \
		            -out "${IR_DATA_DIR}/SC_RefDB" >& /dev/null
	loginfo "Done"

popTmpDir

