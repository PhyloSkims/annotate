#!/bin/bash
#
#                           BUILD REFERENCE : From SwissProt
#
#========================================================================================

# -- CAUTION -- Works as long than the script 
#               is not called through a symlink

THIS_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${THIS_DIR}/../../../scripts/bash_init.sh"

SPDIR="$CDS_DATA_DIR/sp_chlorodb"
SPGENESDIR="$SPDIR/genes"

CORE_GENES="accD"
CORE_GENES="$CORE_GENES atpA atpB atpE atpF atpH atpI"
CORE_GENES="$CORE_GENES ccsA"
CORE_GENES="$CORE_GENES cemA"
CORE_GENES="$CORE_GENES clpP"
CORE_GENES="$CORE_GENES infA"
CORE_GENES="$CORE_GENES matK"
CORE_GENES="$CORE_GENES ndhA ndhB ndhC ndhD ndhE ndhF ndhG ndhH ndhI ndhJ ndhK"
CORE_GENES="$CORE_GENES petA petB petD petG petL petN"
CORE_GENES="$CORE_GENES psaA psaB psaC psaI psaJ psbA"
CORE_GENES="$CORE_GENES psbB psbC psbD psbE psbF psbH psbI psbJ psbK psbL psbM psbN psbT psbZ"
CORE_GENES="$CORE_GENES rbcL"
CORE_GENES="$CORE_GENES rpl14 rpl16 rpl2 rpl20 rpl22 rpl23 rpl32 rpl33 rpl36"
CORE_GENES="$CORE_GENES rpoA rpoB rpoC1 rpoC2"
CORE_GENES="$CORE_GENES rps2 rps3 rps4 rps7 rps8 rps11 rps14 rps15 rps16 rps18 rps19"
CORE_GENES="$CORE_GENES ycf1 ycf2 ycf3 ycf4"

RPS12_GENE="rps12"

function download_swissprot() {
    URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz"

    curl $URL
}

function extract_chloro_entries() {
    awk -F'\n' '
        BEGIN {RS="//\n"; ORS=RS; OFS="\n"} 
        /OC   Eukaryota; Viridiplantae; Streptophyta;/ && /OG   Plastid; Chloroplast./ && !/DE   Flags: Fragment;/ {print $0}
    ' $1
}

function extract_chloro_gene_entries() {

    awk -v gene=$1 -F'\n' '
        BEGIN {RS="//\n"; ORS=RS; OFS="\n"} 
        ($0 ~ gene"_") && /OC   Eukaryota; Viridiplantae; Streptophyta;/ && /OG   Plastid; Chloroplast./ {print $0}
    ' $2
}

function extract_chloro_gene_frg() {

    awk -F'\n' '
        BEGIN {RS="//\n"; ORS=RS; OFS="\n"} 
        /DE   Flags: Fragment;/ {print $0}
    ' $1
}

function extract_chloro_gene_ac() {

    awk -v ac=$1 -F'\n' '
        BEGIN {RS="//\n"; ORS=RS; OFS="\n"} 
        ($0 ~ "AC   " ac ";") && /OC   Eukaryota; Viridiplantae; Streptophyta;/ && /OG   Plastid; Chloroplast./ {print $0}
    ' $2
}

function extract_fasta_protein() {
    $AwkCmd '
        /^ID/ {ID=$2} 
        /^AC/ {AC=$2; sub(";","",AC)} 
        /^DR   EMBL;/ {
            EMBL=$4; 
            sub(";","",EMBL); 
            sub("-","xxx",EMBL);
            if ( EMBL != "xxx" ) {
                EBI="curl \"https://www.ebi.ac.uk/ena/browser/api/embl/" EMBL "?download=true\" | egrep \"FT +/gene=\" | cut -d \"=\" -f 2 | tr -d \"\\\"\""
                EBI | getline GENE
                close(EBI)
            } else {
                GENE = "xxx"
            }
            } 
        /^   / {gsub(/ /,"",$0); SEQ=SEQ $0}
        /^\/\// {
            if (GENE != "xxx") {
                print ">"ID,"SP_AC="AC";","EMBL_AC="EMBL";","gene="GENE";";
                print SEQ 
            }
            ID="xxx";
            AC="xxx"; 
            EMBL="xxx"; 
            GENE="xxx"
            SEQ=""
        }
        ' $1 \
    | formatfasta
}

function rename_rules() {
    egrep "^>" $1 \
    | sed 's/^>//' \
    | sed 's/_/=/' \
    | tr -d ";" \
    | awk -F"=" '{print $1,$NF}' \
    | sort \
    | uniq -c \
    | awk '{
        x=$1
        $1=$2
        $2=sprintf("%06d",x)
        print $0}' \
    | sort -rn \
    | sed 's/ /=/' \
    | sed 's/ /=/' \
    | awk -F "=" '
        ($1 == current) {
            sub(/^0+/,"",occurrence)
            print "# based on",occurrence,"observations renames",$NF,"in",gene
            print "/^>" current "/ s@gene=" $NF "@gene=" gene "@"
        }
        ($1 != current) {
            current= $1
            gene=$NF 
            occurrence=$2
        }
        '
}

function rename_genes() {
    local n=1
    local input=$1
    cat $input > __tmp__$$__fasta
    rules=${input/fst/rules}
    #echo $rule 1>&2

    rename_rules __tmp__$$__fasta > $rules.$n
    while [[ -s $rules.$n ]] ; do
        sed -f $rules.$n __tmp__$$__fasta > __tmp2__$$__fasta
        mv __tmp2__$$__fasta __tmp__$$__fasta
        ((n++))
        rename_rules __tmp__$$__fasta > $rules.$n
    done

    cat __tmp__$$__fasta
    rm __tmp__$$__fasta
}


function clean_strange_gene_name() {
    local input=$1
    local output=${input/.fst/_sp_ebi.genes}
    local to_keep=${input/.fst/_to_keep.lst}
    local to_be_removed=${input/.fst/_to_be_removed.lst}

    grep '^>' $input \
    | sed 's/^>//' \
    | sed 's/_/=/' \
    | sed 's/;$//' \
    | awk -F'=' '{print $NF,$1}' \
    | sort \
    | uniq -c \
    | $AwkCmd '{
        x=$1
        $1=$2
        $2=sprintf("%06d",x)
        print $0}' \
    | sort -r > $output

    $AwkCmd '
        ($1!=current) {current=$1; print $3}
        ' $output > $to_keep

    $AwkCmd '
        ($1==current) {print $3} 
        ($1!=current) {current=$1}
        ' $output > $to_be_removed

    filter_sp_fasta_db $to_be_removed $input > ${input/.fst/_strange_genes.fst}
    filter_sp_fasta_db $to_keep $input
}

function filter_sp_fasta_db() {
    gene_pattern="$(echo $(cat $1) | tr ' ' '|')"

    $AwkCmd -F "_" -v gene_pattern=$gene_pattern '
        /^>/ && (gene ~ gene_pattern) {
            print entry
        }

        /^>/ {
            entry = $0
            gene = $1
            sub(/^>/, "", gene)
        }

        !/^>/ {
            entry = entry "\n" $0
        }

        END {
            if (gene ~ gene_pattern)
                print entry
        }
    ' $2

}

function filter_out_sp_fasta_db() {
    gene_pattern="$(echo $(cat $1) | tr ' ' '|')"

    $AwkCmd -v gene_pattern=$gene_pattern '
        /^>/ && (gene !~ gene_pattern) {
            print entry
        }

        /^>/ {
            entry = $0
            gene = $1
            sub(/^>/, "", gene)
        }

        !/^>/ {
            entry = entry "\n" $0
        }

        END {
            if (gene !~ gene_pattern)
                print entry
        }
    ' $2

}

function split_by_gene() {
    $AwkCmd -F "=" -v SPGENES=$SPGENESDIR '
        function writegene(gene,entry) {
            outputdir = SPGENES "/" gene 
            output = outputdir "/" gene ".fst"
            system("mkdir -p " outputdir)
            print entry >> output
            close(output)
        }

        /^>/ && gene!="" {
            writegene(gene,entry)
        }

        /^>/ {
            entry = $0
            gene = $NF
            sub(/;$/, "", gene)
        }

        !/^>/ {
            entry = entry "\n" $0
        }

        END {
            writegene(gene,entry)
        }
    ' $1
}

function dereplicate() {
    local CDHIT_ID=0.95
    local CDHIT_DELTA=0.95

 	local gene="${1}"
 	local fastain="${gene}/${gene}.fst"
 	local cdhitout="${gene}/${gene}.cdhit.fst"

     	cd-hit 	-i "${fastain}" \
 	    	-o "${cdhitout}" \
 	    	-c ${CDHIT_DELTA} \
 	    	-G 1 \
			-g 1 \
			-aL 0.95 \
			-s ${CDHIT_ID} \
 	    	-b 350 -p 1 \
 	    	-d 0 -n 3

 	local fasta1="${gene}/${gene}.1l.fst"
}

function dereplicate_genes() {
    pushd $SPGENESDIR
    for g in * ; do 
        dereplicate $g ; 
    done
    popd
}

function buildGeneBlastDB() {
 	local gene="${1}"
 	local fastain="${gene}/${gene}.cdhit.fst"

    loginfo "  formatting Blast $gene DB"
    timeoutcmd 300 makeblastdb -dbtype prot -in ${fastain} >& /dev/null
	loginfo "Done"
}

function buildBlastDBs() {
    pushd $SPGENESDIR
    for g in * ; do 
        buildGeneBlastDB $g ; 
    done
    popd
}

function list_shell_genes() {
    pushd $SPDIR

    ls genes \
    | grep -v '\.' \
    | egrep -iv $(tr " " "|" <<< $CORE_GENES) \
    | grep -iv '^orf' \
    | grep -iv "$RPS12_GENE"

    popd
}

function list_dust_genes() {
    pushd $SPDIR 1>&2

    ls genes \
    | grep -v '\.' \
    | egrep -iv $(tr " " "|" <<< $CORE_GENES) \
    | grep -i '^orf' 

    popd 1>&2
}

function build_core_libraries() {
    pushd $SPDIR 1>&2

    rm -rf core
    mkdir -p core

    for gene in $CORE_GENES ; do 
        cp genes/$gene/$gene.cdhit.fst core/$gene.fst
        cp genes/$gene/$gene.cdhit.fst.phr core/$gene.fst.phr
        cp genes/$gene/$gene.cdhit.fst.pin core/$gene.fst.pin
        cp genes/$gene/$gene.cdhit.fst.psq core/$gene.fst.psq
    done

    popd 1>&2
}

function build_rps12_library() {
    pushd $SPDIR 1>&2

    local gene=$RPS12_GENE

    rm -rf RPS12
    mkdir -p RPS12
    
    cp genes/$gene/$gene.cdhit.fst RPS12/$gene.fst
    cp genes/$gene/$gene.cdhit.fst.phr RPS12/$gene.fst.phr
    cp genes/$gene/$gene.cdhit.fst.pin RPS12/$gene.fst.pin
    cp genes/$gene/$gene.cdhit.fst.psq RPS12/$gene.fst.psq

    popd 1>&2
}

function build_shell_libraries() {
    pushd $SPDIR 1>&2

    rm -rf shell
    mkdir -p shell

    for gene in $(list_shell_genes) ; do 
        cp genes/$gene/$gene.cdhit.fst shell/$gene.fst
        cp genes/$gene/$gene.cdhit.fst.phr shell/$gene.fst.phr
        cp genes/$gene/$gene.cdhit.fst.pin shell/$gene.fst.pin
        cp genes/$gene/$gene.cdhit.fst.psq shell/$gene.fst.psq
    done

    popd 1>&2
}

function build_dust_libraries() {
    pushd $SPDIR 1>&2

    rm -rf dust
    mkdir -p dust

    for gene in $(list_dust_genes) ; do 
        cp genes/$gene/$gene.cdhit.fst dust/$gene.fst
        cp genes/$gene/$gene.cdhit.fst.phr dust/$gene.fst.phr
        cp genes/$gene/$gene.cdhit.fst.pin dust/$gene.fst.pin
        cp genes/$gene/$gene.cdhit.fst.psq dust/$gene.fst.psq
    done

    popd 1>&2
}

function get_product_line() {
    pushd $SPGENESDIR 1>&2

    local gene=$1
    local spac=$(head -1 $gene/$gene.cdhit.fst \
                 | awk '
                     { AC=$2; 
                       sub(/^SP_AC=/,"",AC); 
                       sub(/;$/,"",AC); 
                       print AC}')

    popd 1>&2

    pushd $SPDIR 1>&2

    extract_chloro_gene_ac $spac rawdata/SP_Chloro.dat \
    | grep "^DE   " \
    | $AwkCmd -v gene=$gene '
        function remove_tails(line) {
            sub(/ *(\{[^}]+\});/,"",line)
            st = index(line,"=")
            return substr(line,st+1)
        } 

        /DE +RecName:/ {
            full = remove_tails($0)
        }
        /DE +Short=/ {
            ns = remove_tails($0)
            if (length(ns) > length(short)) {
                short = ns
            }
        }
        /DE +EC=/ {
            ec = remove_tails($0)
        }

        END {
            if (length(short) > 10) {
                product = short
            } else {
                product = full
            }

            if (ec != "") {
                product = product " (EC:" ec ")"
            }

            if (product == "") {
                product = "Hypothetical protein of unknown function"
            } 

            gsub(/ /,"_",product)

            print gene,gene,"--","--","--",product
        }
      '
    popd 1>&2
}

function build_annotate_lst() {
    pushd $SPGENESDIR 1>&2

    for gene in * ; do 
        get_product_line $gene
    done

    popd 1>&2
}

pushd $SPDIR

####
#
# Download and prepare raw library from Swissprot FTP site
#
###

rm -rf rawdata
mkdir -p rawdata

pushd rawdata

download_swissprot | extract_chloro_entries > SP_Chloro.dat

extract_fasta_protein SP_Chloro.dat > SP_Chloro_gene_db.fst

popd

####
#
# Clean swiss-prot fasta file for gene name annotation
#
###

pushd rawdata

rename_genes SP_Chloro_gene_db.fst > SP_Chloro_gene_db.clean_name.fst

clean_strange_gene_name SP_Chloro_gene_db.clean_name.fst \
                      > SP_Chloro_gene_db.good_gene.fst

popd

####
#
# Prepare the database for all genes
#
###

rm -rf genes
split_by_gene rawdata/SP_Chloro_gene_db.good_gene.fst

dereplicate_genes

buildBlastDBs

####
#
# Prepare the differente gene databases for CDS annotation
#
###

build_core_libraries
build_shell_libraries
build_dust_libraries
build_rps12_library

####
#
# Build Annotation file
#
###

build_annotate_lst > Annot.lst

# ls ../../chlorodb/core/ | grep -v '\.' | sort > core_genes.lst
# ls | grep -v '\.' | awk '{print tolower($1),$1}' | sort > sp_genes.lst
# join -a 1 -e xxxx core_genes.lst sp_genes.lst  > join_core.lst

popd