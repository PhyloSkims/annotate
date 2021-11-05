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

if (( $# > 1 )) ; then
    TEMP=$2
else
    TEMP=""
fi

DBROOT="$CDS_DATA_DIR/chlorodb/RPS12"
RPS12DB="${DBROOT}/RPS12_DB.clean.fst"
DELTA=50

SEQLEN=$(seqlength "${QUERY}")
SEQUENCE=$(readfirstfastaseq "${QUERY}") 

pushTmpDir ORG.RPS12

# localize the gene on the chloroplast genome using blast
loginfo "Locating RPS12 gene by similarity..."


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
                ($2 == PREV_CDS) { HSPs = HSPs "\n" $0;} 
                
                (OUT < 20) && ($2 != PREV_CDS) && (BEST_EVAL < (1e-20 + 0.0)) { 
                    if (PREV_CDS) print HSPs; 
                    HSPs = $0;  
                    BEST_EVAL = 1; 
                    PREV_CDS = $2; 
                    OUT++ 
                } 
                {PREV_CDS = $2;}
                
                (BEST_EVAL > ($11 + 0.0)) {BEST_EVAL = ($11 + 0.0)} 
          ' > "rps12_locate.hsps"

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
        | $AwkCmd  'function overlap(x1,y1,x2,y2) {  
                      return (((x1+0 <= x2+0) && ((y1+1) >= x2+0)) || 
                              ((x2+0 <= x1+0) && ((y2+1) >= x1+0))) 
                    } 
                    function min(a,b) {return (a <= b) ? a:b } 
                    function max(a,b) {return (a >= b) ? a:b }
                    (NR==1) {i=0
                             frg[i]=$0 
                            } 
                    (x1 && y1) { 
                        if (overlap(x1,y1,$1,$2)) {
                            $1 = min(x1,$1) 
                            $2 = max(y1,$2) 
                            if (overlap(v1,w1,$3,$4)) { 
                                $3 = min(v1,$3) 
                                $4 = max(w1,$4) 
                                } 
                            } 
                            else i++ 
                    } 
                    (x1 && y1) { 
                       frg[i] = $0
                    } 
                    { x1 = $1 
                      y1 = $2 
                      v1 = $3 
                      w1 = $4 
                    } 
                    END {
                        for (j = 0; j <= i; j++) {
                           print frg[j]
                        }
                    }
        ' \
        | sort -nk 3 \
        | $AwkCmd '($3 != old3 || $4 != old4) { 
                          i++
                          old3=$3 
                          old4=$4 
                        } 
                    {print $0,i} 
                  ' \
        | sort -nk 6  \
        | $AwkCmd 'function min(a,b) {return (a <= b) ? a:b } 
                   (old6 == 1) {
                        print old
                        oldprint = 1
                   }
                   ((old6 == 2 && $6==2) ||
                    full == 1) {
                       print old
                       full = 0
                   }
                   (((old6 == 2 && $6==3) ||
                    (old6 == 3 && $6==2)) && full != 1) {
                        $1 = old1
                        $6 = min(old6,$6)
                        full = 1
                   }
                   END {print old}
                   {
                     old = $0
                     old1 = $1
                     old6= $6
                   }'  \
        | $AwkCmd -v delta="$DELTA" \
                  -v seqlen="$SEQLEN" \
                  -v chloro="$SEQUENCE" \
             'function min(a,b) {return (a <= b) ? a:b } 
              function max(a,b) {return (a >= b) ? a:b }
              function rev(s) {
                x = ""
                for (i=length(s);i!=0;i--) 
                    x=x substr(s,i,1)
                return x
              }           
              function swapchar(s,a,b) {
                gsub(a,"@",s)
                gsub(b,a,s)
                gsub(/@/,b,s)
                return s
              }
              function revcomp(s) {
                s = swapchar(s,"A","T")
                s = swapchar(s,"C","G")
                s = swapchar(s,"M","K")
                s = swapchar(s,"R","Y")
                s = swapchar(s,"W","S")
                s = swapchar(s,"B","V")
                s = swapchar(s,"D","H")
                s = swapchar(s,"a","t")
                s = swapchar(s,"c","g")
                s = swapchar(s,"m","k")
                s = swapchar(s,"r","y")
                s = swapchar(s,"w","s")
                s = swapchar(s,"b","v")
                s = swapchar(s,"d","h")
                return rev(s)
              }
              {   from = max(1,$1 - delta)
                  to   = min($2 + delta,seqlen)
                  sequence = substr(chloro,from,to-from+1)
                  if ($5 == "R") sequence = revcomp(sequence)
                  nparts[$6]+=1
                  n = nparts[$6]
                  parts[$6][n][1] = from
                  parts[$6][n][2] = to
                  parts[$6][n][3] = $3
                  parts[$6][n][4] = $4
                  parts[$6][n][5] = $5
                  parts[$6][n][6] = $6
                  parts[$6][n][7] = sequence
              }
              END {
                  l =  length(parts)
                  if (l==1) {
                      n = nparts[1]
                      for (i =1; i <= n; i++) {
                          print ">RPS12_" i,"parts=1; limit=" length(parts[1][i][7]) + 1 \
                                "; from1=" parts[1][i][1] \
                                "; to1=" parts[1][i][2] "; strand1=" parts[1][i][5] \
                                ";" > "rps12_fragments_" i ".fasta"
                          print parts[1][i][7] \
                                > "rps12_fragments_" i ".fasta"
                      }
                    }

                  if (l==2) {
                      n1 = nparts[1]
                      n2 = nparts[2]
                      for (i =1; i <= n1; i++) 
                        for (j =1; j <= n2; j++) {
                          k = (i-1)*n2+j
                          print ">RPS12_" k,"parts=2", \
                                "limit=" (length(parts[1][i][7]) + 10 + 1) \
                                "; from1=" parts[1][i][1] "; to1=" parts[1][i][2] "; strand1=" parts[1][i][5] \
                                "; from2=" parts[2][j][1] "; to2=" parts[2][j][2] "; strand2=" parts[2][j][5] \
                                ";" > "rps12_fragments_" k ".fasta"
                          print parts[1][i][7] "nnnnnnnnnn" parts[2][j][7] \
                                > "rps12_fragments_" k ".fasta"

                      }
                    }
               }
             '

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
        cp $DBROOT/Annot.lst RPS12

        for f in rps12_fragments_*.fasta ; do
            tcsh -f ${PROG_DIR}/do_exonerate.csh \
                $f \
                "RPS12/rps12.fasta" \
                $DBROOT/../models $(pwd)
        done

    #
    # Rewrite the coordinates of the genes on the extracted
    # fragment to the chloroplaste genome coordinates 
    #


        for f in *.res ; do
            mv $f $f.ori
            if [[ -z "$TEMP" ]] ; then
                dest="/dev/stdout"
            else
                dest="$TEMP/$f"
            fi
            header=$(head -1 ${f/.rps12.res/.fasta})
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
                '
                function convert1p(p) {
                  if (p+0 < L2) {
                    I = 1
                    if (S1=="F") {
                        S = 1
                        B = F1
                    } else {
                        S = -1
                        B = T1
                    }
                  } else {
                    I = L2
                    if (S2=="F") {
                        S = 1
                        B = F2
                    } else {
                        S = -1
                        B = T2
                    }
                  }
                  return S*(p - I) + B
                }
                function convert(p1,p2) {
                  p1  = convert1p(p1)
                  p2  = convert1p(p2)
                  if (p1 < p2)
                     res = p1 ".." p2
                  else
                     res = "complement(" p2 ".." p1 ")"
                  return res
                }
                /[0-9]+\.\.[0-9]+/ {
                    s = $0
                    r = $0
                    while (length(s) > 0) {
                        match(s,/[0-9]+\.\.[0-9]+/)
                        range = substr(s,RSTART,RLENGTH)
                        s = substr(s,RSTART+RLENGTH+1)
                        match(range,/^[0-9]+/)
                        from = substr(range,RSTART,RLENGTH)
                        match(range,/[0-9]+$/)
                        to = substr(range,RSTART,RLENGTH)
                        sub(range,convert(from,to),r)
                    }
                    $0=r
                 }
                {print $0}
                ' \
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
            ' \
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
        ' > $dest
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



