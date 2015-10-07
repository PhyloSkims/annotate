#!/bin/bash
#
#                           NORMALISATION D'UN PLASTIDE
#
#========================================================================================
# Ce programme dispose de 4 fonctions pour traiter les donnees fasta issues de genbank
# - seqlength : compte le nombre de paire de base du fichier
# ex : seqlength $1
#
# - cutseq : permet de couper un morceau de la sequence
# cutseq [x] [y]
# [x] : coordonne du debut de la sequence a couper
# [y] : coordonne de la fin de la sequence a couper
# ex : cutseq $1 10 100
#
#
# - revcomp : donne le brin reverse
# ex : $1 | revcomp
#
# - formatfasta : permet de coller a la suite plusieurs morceaux de sequence au moment de 
#   la reecriture
# - joinfasta : enleve les titres au moment de la reecriture du fichier et renvoie les 
#   informations dans la fonction formatfasta
# ex : joinfasta $1
#
#========================================================================================
# Pour lancer le programme, utiliser les commandes : 
# chmod +x normalize_plastid.sh
#./normalize_plastid.sh [fichier].fasta
#
# ex : seqlength $1
#
# cutseq $1 [x] [y]
# [x]:coordonne du debut [y]:coordonne de la fin de la sequence a couper
# ex : cutseq $1 10 100
#
# ex : $1 | revcomp
#
# ex : joinfasta $1
#========================================================================================

FINDCDS=`dirname $0`/findcds
REPSEEK=`dirname $0`/repseek

# s'alimente avec un fichier.fasta
# $3 : nb de caractere du fichier, t : nb de caractere du titre, 
# $1+1 : nb de retour chariot du fichier
function seqlength {
  cat $1 | \
  wc |\
  awk -v t="`head -1 $1 | wc -c`" '{print $3 - t - $1 + 1}'
}


# selectionne une sequence parmi le fichier
# $2 : debut de la sequence a couper, $3 : fin de la sequence a couper
function cutseq {
	awk -v from=$2 -v end=$3 'function printfasta(seq) {            \
									seqlen=length(seq);             \
									for (i=1; i <= seqlen; i+=60)    \
									   print substr(seq,i,60);      \
								}                                   \
																	\
								/^>/   {print $0}                   \
								! /^>/ {seq=seq$0}                  \
								END {printfasta(substr(seq,from,end-from+1))}' $1
}

# donne le brin reverse de la sequence  : 
# la sous-fonction comp reecrit la sequence a l'envers 
# la sous-fonction rev remplace les bases par leurs bases associe
# la sous-fonction revcomp reprend les deux precedente
function revcomp {
	awk 'function printfasta(seq) {          					    \
			seqlen=length(seq);           						    \
			for (i=1; i <= seqlen; i+=60)      						\
			  print substr(seq,i,60);         						\
		 }                                      					\
		function comp(seq) {                   						\
			"echo "seq" | tr acgtACGT tgcaTGCA " | getline res; 	\
			return res;                        						\
		}                                      						\
		function rev(seq) {                    						\
			"echo "seq" | rev " | getline res; 						\
			return res;                        						\
		}                                      						\
		function revcomp(seq) {                						\
			res=rev(comp(seq));                						\
			return res;                        						\
		}                                      						\
																	\
		/^>/   {print $0}                    					    \
		! /^>/ {seq=seq$0}                     						\
		END {printfasta(revcomp(seq))}' $1
}


function formatfasta {
	awk  'function printfasta(seq) {                                \
									seqlen=length(seq);             \
									for (i=1; i <= seqlen; i+=60)   \
									   print substr(seq,i,60);      \
								   }                                \
								/^>/   { print $0 }                 \
								! /^>/ { seq=seq $0 }               \
								END    { printfasta(seq)}' $1
}



# colle bout a bout deux sequence en mettant le meme nombre de paire de base par ligne
# sur le fichier, ici regle a 60
# enleve les titres intermediaire entre deux sequences recollees si il y en a
function joinfasta {
	awk '(NR==1 && /^>/) {print $0}                                 \
	     ! /^>/          {print $0}' $1 |                           \
		 formatfasta
}

# recupere les informations issues du programme repseek avec l'origine des deux 
# IR et leur taille
function lookforIR {
	repseek -c -p 0.001 $1|                                         \
     grep 'Distant.inv'|                                            \
     sort -n -k4              |                                     \
     tail -1                  |                                     \
     awk '{print $7}'         |                                     \
     sed 's/-/ /g'
}

# determine si le fragment analyse doit etre recolle en forward ou reverse dans le 
# nouveau .fasta
function maxCDS {
	${FINDCDS} -F $1 -c -l 150 | \
	awk '/^[^#]/ && ($2 == "Watson") { Watson+=$5-$4+1} \
	     /^[^#]/ && ($2 == "Crick")  {Crick+=$5-$4+1} \
	     END                         {print Watson - Crick}'
}

#
# Exemple results from repseek
#
#	$1			$2		$3		$4		$5		$6		$7						$8		$9			$10		$11	 $12
# Class			pos_r1	pos_r2	len_r1	len_r2	Delta	Seed					ident	score		Rmean	Rmod frac
# Distant.inv	86608	130934	25319	25319	19007	86608-130934-25319-2.01	100.000	25262.43	2.01	2	 0.99

# Test where the sequence is cut

# Definie les variables utilisé : le debut, fin et taille des deux IR
genome=$1

genome_length=`seqlength $1`
IRS=(`lookforIR ${genome}`)

posIR1=${IRS[0]}
posIR2=${IRS[1]}
lenIR=${IRS[2]}

let "endIR2=$posIR2 + $lenIR - 1"
let "endIR1=$posIR1 + $lenIR - 1"

# Defini la coupe a adopter en fonction de : 
# - la position de la fin de l'IR2 par rapport a la sequence total, pour identifier une 
#   coupe au sein d'une IR 
# Le programme repseek considere toujours que la position maximal de la fin d'IR2 ne peut 
# pas depasser celle de la sequence, que l'IR2 soit coupe ou non. Donc dans tous les cas
# on choisit de recouper la sequence a mi-distance entre la fin de l'IR1 et le debut de 
# l'IR2
if (( endIR2 == genome_length )) ; then	
	tmpfasta1="tmp_$$_1.fasta"
	tmpfasta2="tmp_$$_2.fasta"
# defini la localisation de la coupure entre les deux IR
	let "posCut=($endIR1+$posIR2)/2"
# realise la coupure du fichier d'entre du nucleotide calcule jusqu'a la fin de la sequence
	cutseq ${genome} ${posCut} ${genome_length} > ${tmpfasta1}
	let "posCut=$posCut-1"
# realise la coupure du fichier d'entre du debut de la sequence jusqu'au nucleotide calcule
	cutseq ${genome} 1 ${posCut} >> ${tmpfasta1}
# ces deux fragment sont rassembles dans un fichier temporaire
	joinfasta ${tmpfasta1} > ${tmpfasta2}
	rm -f ${tmpfasta1}
	genome=${tmpfasta2}
# recalcul la nouvelle position des IR
	IRS=(`lookforIR ${genome}`)
	posIR1=${IRS[0]}
	posIR2=${IRS[1]}
	lenIR=${IRS[2]}
	let "endIR2=$posIR2 + $lenIR - 1"
	let "endIR1=$posIR1 + $lenIR - 1"
fi

tmpIR1="tmp_$$_IR1.fasta"		
tmpIR2="tmp_$$_IR2.fasta"		

#enregistre les deux fragments IR1 et IR2 complet
cutseq ${genome} ${posIR1} ${endIR1} > ${tmpIR1}
cutseq ${genome} ${posIR2} ${endIR2} > ${tmpIR2}

let "lenSC1=$posIR1 -1 + ($genome_length - endIR2)"
let "lenSC2=$posIR2 - $endIR1"

tmpLSC="tmp_$$_LSC.fasta"		
tmpSSC="tmp_$$_SSC.fasta"		

# Defini la coupe a adopter en fonction de : 
# - la taille de la SC1 par rapport a la taille de la SC2, pour identifier une 
#   coupe au sein d'une SC. La coupe a lieu au sein de la SC1, le but est d'identifier la 
#   LSC et la SSC parmis les SC1 et la SC2
# si la SC1 est plus grande que la SC2, alors la SC1 est la LSC et la coupe a eu lieu 
# dans la LSC
if (( lenSC1 > lenSC2 )); then
# defini le debut de la LSC
	let "beginLSC=$endIR2+1"
	cutseq ${genome} ${beginLSC} ${genome_length} > ${tmpLSC}
# defini la fin de la LSC
	let "endLSC=$posIR1-1"
	cutseq ${genome} 1 ${endLSC} >> ${tmpLSC}
	tmpfasta1="tmp_$$_1.fasta"
# rejoint les deux morceaux pour former la LSC
	joinfasta ${tmpLSC} > ${tmpfasta1}
	mv ${tmpfasta1} ${tmpLSC}

# donc la SC2 est la SSC, 
# definit l'origine et la fin a couper pour avoir le fragment SSC
	let "beginSSC=$endIR1+1"
	let "endSSC=$posIR2-1"
	cutseq ${genome} ${beginSSC} ${endSSC} > ${tmpSSC}
	
	tmp=${tmpIR1}
	tmpIR1=${tmpIR2}
	tmpIR2=${tmp}
else
# sinon la SC2 est la LSC, et la coupe a eu lieu dans la SSC
# definit l'origine et la fin a couper pour avoir le fragment LSC
	let "beginLSC=$endIR1+1"
	let "endLSC=$posIR2-1"
	cutseq ${genome} ${beginLSC} ${endLSC} > ${tmpLSC}

# definit le debut de la SSC et coupe la premiere partie
	let "beginSSC=$endIR2+1"
	cutseq ${genome} ${beginSSC} ${genome_length} > ${tmpSSC}
# definit la fin de la SSC	et coupe la seconde partie 
	let "endSSC=$posIR1-1"
	cutseq ${genome} 1 ${endSSC} >> ${tmpSSC}
# joint les deux parties afin de reformer la SSC
	joinfasta ${tmpSSC} > ${tmpfasta1}
	mv ${tmpfasta1} ${tmpSSC}
fi




if [[ ! -z $tmpfasta2 ]]; then
	rm -f $tmpfasta2
fi
# determine si les fragments doivent etre recolle en reverse ou forward
maxSSC=`maxCDS ${tmpSSC}`

# si maxSSC est negatif, le rapport Watson - Crick est negatif, le fragment est 
# donc reverse
if (( maxSSC < 0 )); then
	revcomp ${tmpSSC} > ${tmpfasta1}
	mv ${tmpfasta1} ${tmpSSC}
fi

maxLSC=`maxCDS ${tmpLSC}`

# si maxLSC est negatif, le rapport Watson - Crick est negatif, le fragment est 
# donc reverse
if (( maxLSC < 0 )); then
	revcomp ${tmpLSC} > ${tmpfasta1}
	mv ${tmpfasta1} ${tmpLSC}
fi

# Les quatre fragments sont recolle ensemble sans erreur de coupure et dans un ordre connus.
cat ${tmpSSC} ${tmpIR1} ${tmpLSC} ${tmpIR2} | joinfasta

exit 0
