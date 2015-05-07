#!/bin/bash

SUMATRA=`dirname $0`/sumatra

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


for seq in `awk '(on==1) && /^>/ {print "";on=0} /^>/ {printf("%s@",$0)} ! /^>/ {on=1;printf($0)}' $1 | tr ' ' '&'`; do
    new=(`annotateCAU ${seq} $1`)
	echo $seq | sed 's/&/ /g' | sed -E 's/>([^ ]+) />'${new[0]}' /' | tr '@' '\n'
done