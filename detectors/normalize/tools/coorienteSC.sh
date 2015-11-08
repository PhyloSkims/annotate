#!/bin/bash
#
#                           BUILD THE COORIENTATION GRAPH
#
# From a set of sequences, this commande build a graph where:
#             - vertices are sequence id
#             - relation is : s oriented in the same way
#
# The same ortientation is defined from the similiarity measured by blast
# The result graph is formated following the tgf format (trivial graph format)
#                 NODE1 NODE2 EDGE_LABEL
#
# Run the command
#
#   coorienteSC.sh <FASTA_FILE> <MINLENGTH> [LOGFILE]
#
#   <FASTA_FILE> : the fasta file containing sequences to analyse
#   <MINLENGTH>  : the cumulative minimum length strand have to share to conclude
#   <LOGFILE>    : an optional logfile name
#
#========================================================================================


# -- CAUTION -- Works as long than the script 
#               is not called through a symlink
THIS_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${THIS_DIR}/../../../scripts/bash_init.sh"


pushTmpDir  ORG.coorienteSC

	if [[ ! -z $3 ]]; then
		openLogFile $3
	fi

	DATA="${CALL_DIR}/${1}"
	MINLENGTH=$2
	BLASTDB=${1/.fasta$/}

	loginfo "Build temporary blast DB..."
	makeblastdb -in "${DATA}" -dbtype nucl -out "${BLASTDB}" >& /dev/null
	loginfo "Done"
	
	loginfo "Running Blast..."
		blastn -db "${BLASTDB}" -query "${DATA}" -outfmt 6 | \
	 	$AwkCmd ' \
	 	      ($4 > 1000) && ($3 > 70) \
	 	      ($1==QUERY) && \
	 	      ($2==SUBJECT) && \
	 	      (($7<$8) && ($9<$10)) || (($7>$8) && ($9>$10)) {LSAME+= ($3/100.*$4)} \
	 	      ($4 > 1000) &&  ($3 > 70) \
	 	      ($1==QUERY) && \
	 	      ($2==SUBJECT) && \
	 	      (($7<$8) && ($9>$10)) || (($7>$8) && ($9<$10)) {LDIFF+= ($3/100.*$4)} \
	 	      (QUERY!="") && \
	 	      (($1!=QUERY) || ($2!=SUBJECT)) {print QUERY,SUBJECT,LSAME,LDIFF,(LSAME>LDIFF)} \
	 	      (($1!=QUERY) || ($2!=SUBJECT)) {QUERY=$1;SUBJECT=$2; \
	 	                                      LSAME=0; \
	 	                                      LDIFF=0; \
	 	                                      if (($4 > 1000) && ($3 > 70)) { \
	 	                                        if ((($7<$8) && ($9<$10)) || \
	 	                                            (($7>$8) && ($9>$10))) {\
	 	                                         LSAME= ($3/100.*$4) } \
	 	                                      else { \
	 	                                         LDIFF= ($3/100.*$4) }} \
	 	                                     } \
	          END {print QUERY,SUBJECT,LSAME,LDIFF,(LSAME>LDIFF)}' | \
     	$AwkCmd -v minlength="${MINLENGTH}" \
     		  ' (($3>minlength) || \
       			($4 > minlength)) && \
       			($3/($4+1) > 2) && \
       			($1!=$2) { if ($1 > $2) \
       						 { print $2,$1,$5 } \
                           else \
                             {print $1,$2,$5}}' | \
        sort | \
        uniq -c | \
        $AwkCmd '($1==2) {print $2,$3,$4}'
	loginfo "Done"
	          
popTmpDir

 