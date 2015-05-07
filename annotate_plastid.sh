#!/bin/bash
#
#

export ORGANNOT_HOME=`dirname $0`

REPSEEK=${ORGANNOT_HOME}/repseek
SUMATRA=${ORGANNOT_HOME}/sumatra
ARAGORN=${ORGANNOT_HOME}/aragorn
WRAPARAGORN=${ORGANNOT_HOME}/aragorn_wrapper.awk
ECOFIND=${ORGANNOT_HOME}/ecofind

function annotateCAU {
	QUERY="$$.query.fasta"
	echo $1 | sed 's/&/ /' | tr '@' '\n' > ${QUERY}
	${SUMATRA} -d -n ${QUERY} $2  2> /dev/null | \
	awk '    {n[$2]+=1;d[$2]+=$3} \
	     END {for (i in n) \
	            print i, n[i],d[i], d[i]/n[i]\
	         }' | \
	sort -rnk4 | \
	egrep '^trn(I|M|fM)' | \
	tail -1 | \
	awk '{print $1,$NF}'
	rm -rf ${QUERY}
}

function gffTRNA {

	${ARAGORN} -w -io -seq $3 | awk -v gid=${1} -f ${WRAPARAGORN}

}

# s'alimente avec un fichier.fasta
# $3 : nb de caractere du fichier, t : nb de caractere du titre,
# $1+1 : nb de retour chariot du fichier
function seqlength {
  cat $1 | \
  wc |\
  awk -v t="`head -1 $1 | wc -c`" '{print $3 - t - $1 + 1}'
}


# recupere les informations issues du programme repseek avec l'origine des deux
# IR et leur taille
function lookforIR {
	${REPSEEK} -c -p 0.001 $1 |                                     \
     grep 'Distant.inv'       |                                     \
     sort -n -k4              |                                     \
     tail -1                  |                                     \
     awk '{print $7}'         |                                     \
     sed 's/-/ /g'
}

# recupere le nom de la sequence analyse
function seqName {
		 head -n1 $1|                                                          \
	     awk '{print $1}' |                                                    \
	     sed 's/^>//' |                                                        \
	     sed -E 's/.*\|([^|]+)\|/\1/'
}

# cree un resume du fichier analyse au format gff
# ex : GFF (NC_***	Repseek	IR1 start	end .	+	.	)
function gffIR {
	lseq=$2
	nom=$1
    lookforIR $3 |                                                             \
    awk -v nom="$nom" -v  lseq="$lseq"                                         \
    	'BEGIN {OFS="\t"}                                                      \
    	{  startIR1=$1;                  \
    	   startIR2=$2; \
    	   endIR1=startIR1 + $3 -1; \
    	   endIR2=startIR2 + $3 -1; \
    	   startSSC=1; \
    	   endSSC=startIR1-1; \
    	   startLSC=endIR1+1; \
    	   endLSC=startIR2-1; \
    	                \
    	   print nom,"RepSeek","misc_feature",startSSC,endSSC,"\.","+","\.","ID=SSC;note=small single copy region";\
    	   print nom,"RepSeek","repeat_region",startIR1,endIR1,"\.","+","\.","ID=IRA;note=inverted repeat A";\
    	   print nom,"RepSeek","misc_feature",startLSC,endLSC,"\.","+","\.","ID=LSC;note=large single copy region";\
    	   print nom,"RepSeek","repeat_region",startIR2,endIR2,"\.","-","\.","ID=IRB;note=inverted repeat B";\
    	}'
}


echo "##gff-version 3"


genome=$1
genome_name=`seqName $1`
genome_length=`seqlength $1`

gffIR   ${genome_name} ${genome_length} ${genome}| grep -v '^ *$'
gffTRNA ${genome_name} ${genome_length} ${genome}| grep -v '^ *$'
