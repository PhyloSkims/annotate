function min(a,b) {return (a <= b) ? a:b } 
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

function extractFirstSequence(file) {
    # Variables
    sequence = ""
    inSequence = 0

    # Lecture du fichier
    while ((getline line < file) > 0) {
        # Ignorer les lignes commençant par ">"
        if (substr(line, 1, 1) == ">") {
            if (inSequence == 1)
                break;  # Sortir si nous avons déjà extrait la première séquence
            else
                inSequence = 1;  # Commencer à extraire la séquence
        } else if (inSequence == 1) {
            # Supprimer les espaces et retours chariot de la séquence
            gsub(/[[:space:]]/, "", line);
            sequence = sequence line;
        }
    }

    # Fermer le fichier
    close(file);

    # Retourner la première séquence
    return sequence;
}

BEGIN {
    chloro = extractFirstSequence(chloro)
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