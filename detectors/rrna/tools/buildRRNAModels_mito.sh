#!/bin/bash
#
#                           BUILD RRNA models
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
		$AwkCmd ' /^>/ {SEQ++;\
		            match($0,"count=[0-9][0-9]*;");\
		            count=substr($0,RSTART,RLENGTH);\
		            $1=$1"_"SEQ;\
		            print $1,count} \
			 !/^>/ {print $0}'
}


function clustering {
	DATA=$1
	rm -rf $DATA
	mkdir $DATA
	sumaclust -t 0.9 $DATA.fasta | \
		fasta1li > $DATA.clust.fasta
	cluster=$(grep '^>' $DATA.clust.fasta | \
	            sed -E 's/.*cluster=([^;]+);.*$/\1/' | \
	            sort -u)
	for c in $cluster; do
		w=$(grep "$c" "${DATA}.clust.fasta" | \
			head -1 | \
			sed -E 's/.*cluster_weight=([^;]+);.*$/\1/')
	    out=$(printf "${DATA}/%05d_%s" $w $c)
        grep -A1  "$c" "${DATA}.clust.fasta" | \
           grep -v -- -- > "$out.fasta"
        muscle -in "$out.fasta" -out "$out.align.fasta"
	done
}

function revcomp {
    $AwkCmd 'function printfasta(seq) {                                  \
            seqlen=length(seq);                                       \
            for (i=1; i <= seqlen; i+=60)                              \
              print substr(seq,i,60);                                 \
         }                                                          \
        function comp(seq) {                                           \
            "echo "seq" | tr acgtACGT tgcaTGCA " | getline res;     \
            close("echo "seq" | tr acgtACGT tgcaTGCA "); \
            return res;                                                \
        }                                                              \
        function rev(seq) {                                            \
            "echo "seq" | rev " | getline res;                         \
            close("echo "seq" | rev ");                                \
            return res;                                                \
        }                                                              \
        function revcomp(seq) {                                        \
            res=rev(comp(seq));                                        \
            return res;                                                \
        }                                                              \
                                                                       \
        (seq) && /^>/ {print head;                                     \
                       printfasta(revcomp(seq));                       \
                       seq=""}                                         \
        /^>/   {head=$0}                                               \
        ! /^>/ {seq=seq$0}                                             \
        END { print head;                                     \
              printfasta(revcomp(seq));                       \
            }' $1
}


pushTmpDir ORG.buildRRNAMito
	loginfo "Tempdir: $(pwd)"

	openLogFile "${RRNA_DATA_DIR}/rRNA_mito_models.log"
	
	loginfo "Selecting Metazoa genebank entries..."
		METAZOA=$(${PROG_DIR}/../../normalize/tools/selectMetazoa.sh $*) 
		loginfo " --> $(echo ${METAZOA} | wc -w) entries selected"
	loginfo "Done"
	
	loginfo "Extracting 12S rRNA sequences..."
		rm -f raw_12S.fasta
		for f in ${METAZOA}; do
			loginfo "Extracting 12S rRNA sequences from ${f}..."
			${PROG_DIR}/extract_ref12S.sh ${f} >> raw_12S.fasta
		done
		loginfo " --> $(fastaCount raw_12S.fasta) retreived sequences"
		dereplicate raw_12S.fasta >> 12S.fasta
		loginfo " --> $(fastaCount 12S.fasta) distinct sequences"
	loginfo "Done"

	loginfo "Clustering 12S rRNA sequences..."
		clustering 12S
	loginfo "Done"
	
	loginfo "Installing 12S rRNA sequences..."
		cp -r 12S 	"${RRNA_DATA_DIR}/RRNA_12S_mito"
	loginfo "Done"


	loginfo "Extracting 16S rRNA sequences..."
		rm -f raw_16S.fasta
		for f in ${METAZOA}; do
			${PROG_DIR}/extract_ref16S.sh ${f} >> raw_16S.fasta
		done
		loginfo " --> $(fastaCount raw_16S.fasta) retreived sequences"
		dereplicate raw_16S.fasta > 16S.fasta
		loginfo " --> $(fastaCount 16S.fasta) distinct sequences"
	loginfo "Done"

	loginfo "Clustering 16S rRNA sequences..."
		clustering 16S
	loginfo "Done"
	
	loginfo "Installing 16S rRNA sequences..."
		cp -r 16S 	"${RRNA_DATA_DIR}/RRNA_16S_mito"
	loginfo "Done"


popTmpDir
