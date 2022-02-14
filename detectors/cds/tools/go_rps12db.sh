#!/bin/bash
#
#                           BUILD REFERENCE : THE RPS12 LIBRARY
#
#========================================================================================

# -- CAUTION -- Works as long than the script 
#               is not called through a symlink

THIS_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${THIS_DIR}/../../../scripts/bash_init.sh"
source "${THIS_DIR}/lib/clusterize_prot.sh"

function extract_rps12() {
$AwkCmd '                                          \
      /^LOCUS/    {LOCUS=$2;}                      \
      /^     [^ ]/ { if (CDS) {                    \
                        print LOCUS "/" feature;   \
                        print "#################"  \
                       }                           \
                    CDS=0;                         \
                   }                               \
      /^     CDS / {CDS=1;                         \
                    $1="";                         \
                    feature=""}                    \
      (CDS)        { sub(/^ */,"",$0);             \
                     sub(/ *$/,"",$0);             \
                     feature=feature $0}           \
    ' \
    | egrep -i '"rps12"'  \
    | $AwkCmd -F"/" '                              \
       function printfasta(seq) {                  \
             seqlen=length(seq);                   \
             for (i=1; i <= seqlen; i+=60)         \
                 print substr(seq,i,60);           \
             }                                     \
                                                   \
                ($1 != current) {current=$1;       \
                                  n=1              \
                                 }                 \
                 {$1=$1 "_rps12_" n;               \
                  n++;                             \
                  delete keys;                     \
                  for (i=3; i<=NF; i++) {          \
                      split($i,key,"=");           \
                      keys[key[1]]=key[2]          \
                  }                                \
                  prot = keys["translation"];      \
                  gsub(/"/,"",prot);               \
                  print ">" $1,"location=" $2 ";"; \
                  printfasta(prot)                 \
                 }                                 \
    ' 
}


pushTmpDir ORG.buildRPS12DB

	RPS12FILE=RPS12_prot.fst	

	openLogFile "${CDS_DATA_DIR}/chlorodb/RPS12_DB.log"

	loginfo "Selecting Viridiplantae genbank entries..."
		VIRIDIPLANTAE=$(${PROG_DIR}/../../normalize/tools/selectViridiplantae.sh $*) 
		loginfo " --> $(echo ${VIRIDIPLANTAE} | wc -w) entries selected"
	loginfo "Done"

    loginfo "Extracting the RPS12 protein sequences from the plants entries..."
        ( for gbk in ${VIRIDIPLANTAE} ; do
            gzcat $gbk | \
                extract_rps12 
          done ) > ${RPS12FILE}
	loginfo "Done"

	loginfo "Installing the RPS12 protein sequence database..."
		
		cp ${RPS12FILE} "${CDS_DATA_DIR}/chlorodb/RPS12_DB.fst"

	loginfo "Done"

popTmpDir

pushd "${CDS_DATA_DIR}/chlorodb"

	loginfo "Clusterizing the RPS12 protein sequence database..."
    rm -rf RPS12_DB.clean.fst
    clusterize RPS12_DB
	loginfo "Done"

    loginfo "  formatting Blast RPS12 DB"
    timeoutcmd 300 makeblastdb -dbtype prot -in RPS12_DB.clean.fst >& /dev/null
	loginfo "Done"



popd

#
# format blast protein database
#




loginfo "Done"

