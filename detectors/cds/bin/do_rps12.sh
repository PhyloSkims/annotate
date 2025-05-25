#!/bin/bash
#
#                      Annotate the RPS12 gene of a plastide genome
#
#========================================================================================
#
# The RPS12 gene is one of the CDS coding for a riboosomal protein
# Depending on the species, the gene is constituted of oone to three exons.
# The exon one is not located close to the others and a trans-splicing is needed
# to reconstruct the spliced mRNA. The exons 2 and eventuually 3 can be located in the
# inverted repeats (IRs) and therfore they can exist in two copies. This can lead to two
# ways to annotate RPS12
# 
#
#  go_rps12.sh <FASTAFILE>
#
#		- <FASTAFILE> : The fasta file containing the normalized genome to annotate 
#
#  Results are printed to the standart output
#
#========================================================================================

# -- CAUTION -- Works as long than the script 
#               is not called through a symlink
THIS_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${THIS_DIR}/../../../scripts/bash_init.sh"

if [[ ! "$1" =~ ^/ ]]; then
	QUERY="${CALL_DIR}/$1"
else
	QUERY="$1"
fi

shift

GENOME_LENGTH="$1"

if (( $# > 1 )) ; then
    TEMP=$2
else
    TEMP=""
fi

DBROOT="$CDS_DATA_DIR/sp_chlorodb/RPS12"
RPS12DB="${DBROOT}/rps12.fst"
DELTA=50

AnnotFile="$CDS_DATA_DIR/sp_chlorodb/Annot.lst"
ModelsDir="$CDS_DATA_DIR/sp_chlorodb/models"

SEQLEN=$(seqlength "${QUERY}")

pushTmpDir ORG.RPS12

# localize the gene on the chloroplast genome using blast
loginfo "Locating RPS12 gene by similarity..."
loginfo "  Considered genome length : $GENOME_LENGTH"

blastx \
    -query ${QUERY} \
    -db ${RPS12DB}  \
    -query_gencode 11 \
    -outfmt 7 \
    | $AwkCmd ' # Blast HSPs are filtered to keep only
                # at maximum the 20 first ones having an 
                # e-value below 1e-20
                BEGIN {BEST_EVAL = 1e-40; 
                       OUT = 0} 
                /^#/ {next} 
                ($2 == PREV_CDS) && (($11 + 0.0) < (1e-5 + 0.0)) { 
                    HSPs = HSPs "\n" $0;
                    } 
                
                (OUT < 20) && ($2 != PREV_CDS) && (BEST_EVAL < (1e-20 + 0.0)) { 
                    if (PREV_CDS) print HSPs; 
                    HSPs = $0;  
                    BEST_EVAL = 1; 
                    PREV_CDS = $2; 
                    OUT++ 
                } 
                {PREV_CDS = $2;}
                
                (BEST_EVAL > ($11 + 0.0)) {BEST_EVAL = ($11 + 0.0)} 
          ' \
    | $AwkCmd -v glength="${GENOME_LENGTH}"  \
              '!(($7 + 0) > (glength +  0) && $8  + (0 > glength +  0))' \
    | $AwkCmd ' BEGIN { FS = OFS = "\t" }

                {
                    subject = $2
                    bitscore = $12 + 0
                    s_start = $9 + 0
                    s_end = $10 + 0

                    # Ajuster l ordre des coordonnées
                    actual_start = (s_start < s_end) ? s_start : s_end
                    actual_end = (s_start < s_end) ? s_end : s_start

                    # Stocker les données
                    count[subject]++
                    idx = count[subject]
                    starts[subject, idx] = actual_start
                    ends[subject, idx] = actual_end
                    scores[subject, idx] = bitscore
                    lines[subject, idx] = $0
                }

                END {
                    for (subject in count) {
                        n = count[subject]
                        
                        # Tri par score (décroissant) puis position (croissante)
                        for (i = 1; i <= n; i++) {
                            for (j = i + 1; j <= n; j++) {
                                if (scores[subject, i] < scores[subject, j] || 
                                (scores[subject, i] == scores[subject, j] && 
                                    starts[subject, i] > starts[subject, j])) {
                                    swap(starts, subject, i, j)
                                    swap(ends, subject, i, j)
                                    swap(scores, subject, i, j)
                                    swap(lines, subject, i, j)
                                }
                            }
                        }

                        # Sélection des HSPs optimaux
                        selected_count = 0
                        delete selected_starts
                        delete selected_ends
                        delete selected_scores
                        delete selected_lines

                        for (i = 1; i <= n; i++) {
                            current_start = starts[subject, i]
                            current_end = ends[subject, i]
                            current_score = scores[subject, i]
                            current_line = lines[subject, i]
                            
                            # Vérifier les chevauchements avec les HSPs sélectionnés
                            overlap_max = 0
                            overlap_indices = ""
                            for (j = 1; j <= selected_count; j++) {
                                os = (current_start > selected_starts[j]) ? current_start : selected_starts[j]
                                oe = (current_end < selected_ends[j]) ? current_end : selected_ends[j]
                                if (os <= oe && (oe - os + 1) >= 15) {
                                    if (selected_scores[j] > overlap_max) overlap_max = selected_scores[j]
                                    overlap_indices = overlap_indices j " "
                                }
                            }

                            # Appliquer les règles de sélection
                            if (overlap_indices == "") {
                                # Aucun chevauchement : sélectionner
                                add_to_selected(current_start, current_end, current_score, current_line)
                            }
                            else if (current_score >= overlap_max) {
                                # Supprimer les HSPs moins bons et conserver les ex æquo
                                new_selected_count = 0
                                delete temp_selected

                                # Garder les non-chevauchants et les ex æquo
                                for (j = 1; j <= selected_count; j++) {
                                    if (!index(overlap_indices, j " ") || selected_scores[j] == current_score) {
                                        temp_selected[++new_selected_count] = j
                                    }
                                }
                                
                                # Mettre à jour la liste sélectionnée
                                selected_count = 0
                                for (j = 1; j <= new_selected_count; j++) {
                                    idx = temp_selected[j]
                                    selected_count++
                                    selected_starts[selected_count] = selected_starts[idx]
                                    selected_ends[selected_count] = selected_ends[idx]
                                    selected_scores[selected_count] = selected_scores[idx]
                                    selected_lines[selected_count] = selected_lines[idx]
                                }
                                add_to_selected(current_start, current_end, current_score, current_line)
                            }
                        }

                        # Affichage des résultats
                        for (j = 1; j <= selected_count; j++) {
                            print selected_lines[j]
                        }
                    }
                }

                function add_to_selected(start, end, score, line) {
                    selected_count++
                    selected_starts[selected_count] = start
                    selected_ends[selected_count] = end
                    selected_scores[selected_count] = score
                    selected_lines[selected_count] = line
                }

                function swap(arr, subject, i, j,    tmp) {
                    tmp = arr[subject, i]
                    arr[subject, i] = arr[subject, j]
                    arr[subject, j] = tmp
                }
                ' \
    > "rps12_locate.hsps"

    #
    # Extracting protein ids from selected blast HSPs
    #

    $AwkCmd '{print $2}' "rps12_locate.hsps" \
        | sort \
        | uniq > "dbsel.txt"

    #
    # Extract corresponding protein sequences
    # from the RPS12 database.
    #

    mkdir -p RPS12
    $AwkCmd -v FILE="dbsel.txt" \
            -f $LIB_DIR/subdb.awk ${RPS12DB} \
            > "RPS12/rps12.fasta"

    cat "rps12_locate.hsps" \
        | $AwkCmd '# Normalizes the writing of the forward and reverse strand matches
                   ($7 <= $8) {print $7,$8,$9,$10,"F"} 
                   ($7 >  $8) {print $8,$7,$9,$10,"R"}' \
        | sort -n  \
        | uniq \
        | $AwkCmd  -f $LIB_DIR/rps12_filter_1.awk  \
        | sort -nk 3  \
        | $AwkCmd '($3 != old3 || $4 != old4) { 
                          i++
                          old3=$3 
                          old4=$4 
                        } 
                    {print $0,i} 
                  ' \
        | sort -nk 6  \
        | $AwkCmd -f $LIB_DIR/rps12_filter_2.awk \
        | $AwkCmd -v delta="$DELTA" \
                  -v seqlen="$SEQLEN" \
                  -v chloro="${QUERY}" \
                  -f $LIB_DIR/rps12_filter_3.awk

    
    nrps12=$(ls -1 rps12_fragments_*.fasta | wc -l)

    if (( nrps12 > 1 )) ; then
        message="$nrps12 versions"
    else
        message="$nrps12 version"
    fi

    loginfo "$message of the gene rps12 detected."

    #
    # Run exonarate on every fragment of chloroplast
    #
    # It should be one or two fragments
    #
        export PASS1_SPEEDUP=0
        nbseq=0
        for fasta in rps12_fragments_*.fasta ; do
            tcsh -f ${PROG_DIR}/do_exonerate.csh \
                Pass2 \
                $fasta \
                "RPS12/rps12.fasta" \
                $AnnotFile \
                $ModelsDir $(pwd)
            ((nbseq=nbseq+1))
        done

    #
    # Rewrite the coordinates of the genes on the extracted
    # fragment to the chloroplaste genome coordinates 
    #

        n=0
        for f in *.res ; do
            loginfo "processing $f"
            xxx=$(cat $f)
            echo -e "\n==============\n$xxx\n==============\n" 1>&2
            
            ((n=n+1))
            mv $f $f.ori
            if [[ -z "$TEMP" ]] ; then
                dest="/dev/stdout"
            else
                dest="$TEMP/$f"
            fi
            loginfo "Destination file $dest"
            header=$(head -1 ${f/.rps12.res/.fasta})
            loginfo "Header: $header"
            L2=$(sed -E 's/^.*limit=([0-9]+);.*$/\1/' <<< $header)
            S1=$(sed -E 's/^.*strand1=(R|F);.*$/\1/' <<< $header)
            S2=$(sed -E 's/^.*strand2=(R|F);.*$/\1/' <<< $header)
            F1=$(sed -E 's/^.*from1=([0-9]+);.*$/\1/' <<< $header)
            F2=$(sed -E 's/^.*from2=([0-9]+);.*$/\1/' <<< $header)
            T1=$(sed -E 's/^.*to1=([0-9]+);.*$/\1/' <<< $header)
            T2=$(sed -E 's/^.*to2=([0-9]+);.*$/\1/' <<< $header)
            cat $f.ori \
            | $AwkCmd -v S1="$S1" -v F1="$F1" -v T1="$T1" \
                      -v S2="$S2" -v F2="$F2" -v T2="$T2" -v L2="$L2" \
                      -f $LIB_DIR/rps12_filter_4.awk \
            | $AwkCmd '
              #
              # Normalize join(complement(A),complement(B),complement(C)) locations
              #      into complement(join(C,B,A)) 
              #
              /join\((complement\([0-9]+\.\.[0-9]+\),)+complement\([0-9]+\.\.[0-9]+\)\)/ \
                {
                    sub(/join\(complement/,"complement(join",$0)
                    gsub(/\),complement\(/,",",$0)
                    match($0,/[0-9]+\.\.[0-9]+(,[0-9]+\.\.[0-9]+)*/)
                    positions=substr($0,RSTART,RLENGTH)
                    n = split(positions,exons,",")
                    for (i=1; i<=n; i++) {
                        if (i > 1)
                            rexons = exons[i] "," rexons
                        else
                            rexons = exons[i]
                    }
                    sub(positions,rexons,$0)
                }
                { print $0}
            '  \
        | $AwkCmd '
          /^FT   [^ ]/ && (length($0) > 80) {
                 n = split($0,parts,",")
                 j = 1
                 for (i = 1; i <= n; i++) {
                     if (length(line) + length(parts[i]) > 78) {
                         print line ","
                         line = "FT                   "
                         j = i
                     } 
                     if (i > j) line = line ","
                     line = line parts[i]
                 }
                 $0 = line
              }
          {print $0}
        ' | \
        $AwkCmd -v n=$n -v nbseq=$nbseq '
            /^FT +\/gene="rps12"/ && (nbseq > 1) {
                sub(/rps12/,"rps12_" n,$0)
            }
            {
                print $0
            }       
        ' | \
        $AwkCmd '
            #
            # Adds the trans_splicing qualifier
            #
            /^FT                   \/translation=/ {
                print "FT                   /trans_splicing"
            }
            {
                print $0
            }       
        '  > "$dest" 
        done   

popTmpDir

exit 0

# NC_010654.fst
# location=complement(join(77925..77967,78465..78700,52867..52980));
# location=join(complement(52867..52980),109583..109818,110316..110358);
# 52837 52980 1 48 R 1
# 77928 77981 113 130 R 3
# 78458 78712 39 132 R 2
# 109571 109825 39 132 F 2
# 110302 110355 113 130 F 3
# /translation="MPTNPQLIRDARQQKKKKRGSRGLQRCPQRRGVCARVYNINPKK
# -->           MPTNPQLIRDARQQKKKKRGSRGLQRCPQRRGVCARVSNINPKK
# ==>           MPTNPQLIRDARQQKKKKRGSRGLQRCPQRRGVCARVSNINPKK
#               PNSALRKVARVRLTSGFEITAYIPGIGHNLQEHSVVLVRGGRVKDLPGVKYRIVRGTL
# -->           PNSALRKVARVRLTSGFEITAYIPGIGHNLQEHSVVLVRGGRVKDLPGVKYRIVRGTL
# ==>           PNSALRKVARVRLTSGFEITAYIPGIGHNLQEHSVVLVRGGRVKDLPGVKYRIVRGTL
#               DAVAVKNRQQGRSSAIWSQKAEKKVIHF"
# -->           DAVAVKNRQQGRSSAIWSQKAEKKVIHF
# ==>           DAVAVKNRQQGRSSAIWSQKAEKKVIHF
# ADL.norm.fasta
# 69300 69425 1 42 R 1
# 97365 97670 36 137 R 2
# 130601 130906 36 137 F 2

# NC_008822
# location=90942..91313;
# 90942 91310 1 123 F 1

# location=complement(join(77925..77967,78465..78700,52867..52980));
# >RPS12_1 parts=2 limit=255; from1=52787; to1=53030; strand1=R; from2=77878; to2=78762; strand2=R;
# join(51..159,312..553,1051..1092)
# location=join(complement(52867..52980),109583..109818,110316..110358);
# >RPS12_2 parts=2 limit=254; from1=52787; to1=53030; strand1=R; from2=109521; to2=110405; strand2=F;
# join\((complement\([0-9]+\.\.[0-9]+\),)+complement\([0-9]+\.\.[0-9]+\)\)

# cat NC_010654.annot.embl | awk -v tag="agser" 'BEGIN {n=0} /\locus_tag/ {n++; sub(/""/, sprintf("\"%s%04d\"",tag,n),$0)} {print $0}' | less

# BRR.chloro_1	NC_018117_rps12_2	96.08	102	4	0	99022	98717	36	137	9e-57	195
# BRR.chloro_1	NC_018117_rps12_2	96.08	102	4	0	141187	141492	36	137	9e-57	195
# BRR.chloro_1	NC_018117_rps12_2	94.87	39	2	0	70611	70495	1	39	2e-16	78.6


# /trans_splicing
