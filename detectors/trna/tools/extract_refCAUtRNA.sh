#!/bin/bash
#
#                           BUILD REFERENCE THE CAU TRNA LIBRARy
#
#========================================================================================

# -- CAUTION -- Works as long than the script 
#               is not called through a symlink
THIS_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${THIS_DIR}/../../../scripts/bash_init.sh"

function taxid {
	egrep '/db_xref="taxon:[0-9]+"' $1 | \
	sed -E 's@ +/db_xref="taxon:([0-9]+)"@\1@'
}

function ac {
	head -1 $1 | $AwkCmd '{print $2}'
}

function definition {
	$AwkCmd '/^DEFINITION/      {on=1}                         \
	     (on==1)            {printf("%s ",$0)}             \
	     (/\.$/ && (on==1)) {on=0;print ""}' $1 |          \
	sed 's/^DEFINITION *//' | \
	sed 's/ *$//'
}

function gb2fasta {
	AC=`ac $1`
	TAXID=`taxid $1`
	DEFINITION=`definition $1`

    echo ">${AC} taxid=${TAXID}; ${DEFINITION}"
    
    $AwkCmd '/^\/\//    {on=0}                                             \
         (on==1)    {print $0}                                         \
         /^ORIGIN / {on=1}' $1 |                                       \
    sed -E 's/^ *[0-9]+ +//'   |                                       \
    sed 's/ //g'
}

function findCAUtrna {

	FASTATMP="$$.genome.fasta"

	gb2fasta $1 > ${FASTATMP}
	aragorn -i -w -seq ${FASTATMP} |                               \
	$AwkCmd '(on==1) && /^ *[0-9]+/ {on=0;print ""}                       \
	     (on==1)                {printf($0)}                          \
	     /\(cat\)$/             {on=1; printf("%s ",$0)}              \
	     END {print ""}' | \
	$AwkCmd '{print $3,$6}'  | \
	sed -E 's/c?\[([0-9]+),([0-9]+)\]/\1 \2/' | \
	sed 's/ /:/g'
	
	rm ${FASTATMP}
}

function trnaAnnotations {
    $AwkCmd '/^ORIGIN/    {on=0} \
         (on==1)      {print $0} \
         /^FEATURE/   {on=1}' $1 | \
    $AwkCmd '/^     [^ ]/ {print ""} \
                      {printf("%s ",$0)} \
         END          {print ""}' | \
    sed 's/^ *//' | \
    sed -E 's/ +/ /g' | \
    grep '^tRNA' | grep '/gene="' | \
    sed -E 's/([0-9]+)\.\.([0-9]+)/\1 \2/g' | \
    sed -E 's/ [0-9]+,[0-9]+ / /g' | \
    grep -v '>' | \
    grep -v '<' | \
    sed -E 's/join\(([0-9]+ [0-9]+)\)/\1/' | \
    sed -E 's/complement\(([0-9]+ [0-9]+)\)/\1/' | \
    sed -E 's/join\(([0-9]+ [0-9]+)\)/\1/' | \
    sed 's/^tRNA *//' |  \
    sed -E 's@([0-9]+) +([0-9]+).*/gene="([^"]+)"@\1 \2 \3@' | \
    $AwkCmd '{print $1,$2,$3}' 
}

function annotateCAU {
    DISTTMP="$$.trna.dist"
	trna=(`echo $1 | sed 's/:/ /g'`) 
	$AwkCmd -v b=${trna[0]} -v e=${trna[1]} \
	    '{printf("sqrt((%d - %d)^2 + (%d - %d)^2)\n",$1,b,$2,e)}' $2 | \
	bc -l | \
	sed 's/\..*$//' > ${DISTTMP}
	paste ${DISTTMP} $2 | sort -nk 1 | head -1 | $AwkCmd '{print $1,$4}'
	rm -f ${DISTTMP}
}

function writeTRNA {
	TAXID=`taxid $1`	
	AC=`ac $1`
	DEFINITION=`definition $1`

	TRNATMP="$$.trna.txt"

	trnaAnnotations $1 > ${TRNATMP}
	ntrna=`wc -l ${TRNATMP} | $AwkCmd '{print $1}'`

	if (( ntrna > 0 )); then 
		trnacau=`findCAUtrna $1`

		for t in $trnacau; do
			AA=(`annotateCAU $t ${TRNATMP}`)
			distance=${AA[0]}
			aa=`echo ${AA[1]} | sed -E 's/(t(rn|RNA)-?)?(I|M|fM).*$/trn\3/'`
		
			if (( distance <= 10 )); then
				echo ">${aa}_${AC} gbac=${AC}; trna=${aa}; taxid=${TAXID}; distance=${distance}; ${DEFINITION}"
				echo "$t" | $AwkCmd -F ':' '{print $3}' 
			fi
		done
	fi

	rm -f ${TRNATMP}
}

pushTmpDir ORG.buildCAUtRNA
	
	for gb in $*; do
		writeTRNA $gb
	done

popTmpDir


