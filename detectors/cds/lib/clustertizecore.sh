#!/bin/bash

THIS_DIR="$(dirname ${BASH_SOURCE[0]})"

source "${THIS_DIR}/../../../scripts/bash_init.sh"

CORELIB="${CDS_DATA_DIR}/chlorodb/core"

CDHIT_ID=0.7
CDHIT_DELTA=0.8


function clusterize() {
 
 	local prot="${1}"
 	local fastain="${prot}.fst"
 	local cdhitout="${prot}.cdhit"
 	local fasta1="${prot}.1l.fst"
 	
 	rm -rf "${prot}"
 	mkdir -p "${prot}"
 	
 	cd-hit 	-i "${fastain}" \
 	    	-o "${prot}/${cdhitout}" \
 	    	-c ${CDHIT_DELTA} \
 	    	-G 1 \
			-g 1 \
			-aL 0.95 \
			-s ${CDHIT_ID} \
 	    	-b 350 -p 1 \
 	    	-d 0 -n 3
 	       
 	fasta1line "${fastain}" > "${prot}/${fasta1}"
 	
 	pushd "$prot"
 	
 	rm -rf "*.cluster.*.ids"
 	$AwkCmd -v prot="$prot" \
 	        '/^>/   {cluster=$2} \
 	         ! /^>/ {sub("\\.\\.*$","",$3); \
 	                 filename=prot".cluster."cluster".ids"; \
 	                 print $3 >> filename ; \
 	                 close(filename) }' "${cdhitout}.clstr"
 	                
 	rm -f "../$prot.clean.fst"
 	                 
 	for ids in *.cluster.*.ids ; do
 		cluster=$(printf "%03d" $(echo "${ids}" | $AwkCmd -F'.' '{print $3}'))
 		size=$(wc -l "$ids" | $AwkCmd  '{print $1}')
 		fsize=$(printf '%05d' $size)

 		alignment="${prot}.cluster.${fsize}.${cluster}.fst"
 	    
 	    if (( size > 1 )) ; then
	 	    egrep -f "$ids" -A 1 "${fasta1}" | \
	 	    egrep -v -- '^--$' | \
	 	    clustalo -i -  > "$alignment"
	 	else
	 	    egrep -f "$ids" -A 1 "${fasta1}" | \
	 	    egrep -v -- '^--$' > "$alignment"
	    fi
	 		
 	    if (( size >= 5 )) ; then
 	    	egrep -f "$ids" -A 1 "${fasta1}" | \
	 	    egrep -v -- '^--$' | \
	 	    formatfasta >> "../$prot.clean.fst"
 	    fi
 	    
 	    rm -f "$ids"
 	    
 	done 
 	
 	rm -f "${fasta1}"
 	rm -f "${cdhitout}"
 	
 	popd
 	       
}
 


pushd $CORELIB

rm -rf *.clean.fst

for prot in *.fst ; do
	prot="${prot/.fst/}"
	clusterize "$prot"
done

popd