#!/bin/csh -f
#
# usage: go_subdb.sh prot.fst pat.txt [deltalen covmin pmax idmin sizmin]
# usage: prot.fst : proteins fasta file
# usage: pat.txt  : text file containing patterns and names for families to extract
# usage: output directory containig subdbs : basename <pat:r>.db
#

unsetenv ORG_SOURCED

setenv ORG_HOME `dirname $0`/../../../../..
source $ORG_HOME/scripts/csh_init.sh

NeedArg 2

set ProtFile = $Argv[1]; Shift
set PatFile  = $Argv[1]; Shift

NeedFile $ProtFile
NeedFile $PatFile

#
# parameters
#

set Delta  = 0.5
set Covmin = 30
set Pmax   = 1e-6
set Idmin  = 30
set Sizmin = 5

if ($#Argv > 0) then
  set Delta = $Argv[1]; Shift
endif

if ($#Argv > 0) then
  set Covmin = $Argv[1]; Shift
endif

if ($#Argv > 0) then
  set Pmax = $Argv[1]; Shift
endif

if ($#Argv > 0) then
  set Idmin = $Argv[1]; Shift
endif

if ($#Argv > 0) then
  set Sizmin = $Argv[1]; Shift
endif

#
# output directory
#

set OutDir = `basename $PatFile:r`.db

if (-d $OutDir) \rm -r $OutDir
mkdir $OutDir

set OutLog = `basename $PatFile:r`.log

echo -n '' > $OutLog

alias Report 'egrep "^>" \!:1 | wc -l | awk -v P=`basename \!:1` -v H=\!:2 '"'{print H,P,"'$1}'"'"' >> $OutLog'

#
# remove entries with bad symbols
#

Notify "cleanup $ProtFile"

Report $ProtFile "init_size"

$AwkCmd -f $LIB_DIR/db.filter.sym.awk $ProtFile > P_$$

Report $ProtFile "cleanup_size"

#
# select by name pattern
#

Notify "select by patterns"

mkdir D_$$
mkdir E_$$
mkdir F_$$

set noms = `awk '{print $1}' $PatFile`

foreach nom ($noms)
  set pat = `egrep "^$nom " $PatFile | awk '{print $2}'`
  $AwkCmd -f $LIB_DIR/db.filter.pat.awk -v PAT="$pat" P_$$ > D_$$/$nom.fst
  set n = `egrep '^>' D_$$/$nom.fst | wc -l`
  Notify "  pattern : $nom : $n"
  Report D_$$/$nom.fst "pattern_filter"
  if ($n <= $Sizmin) \rm -f D_$$/$nom.fst
end

set ok = `ls D_$$ | wc -l`
if ($ok == 0) goto fin

#
# select by length
#

Notify "select by length"

foreach f (D_$$/*.fst) 
  set nom = `basename $f:r`
  $AwkCmd -f $LIB_DIR/db.getlen.awk $f > L_$$
  $LIB_DIR/db.filter.len.r L_$$ $Delta |\
    $AwkCmd '($NF == "TRUE") {print $2}' > M_$$
  $AwkCmd -v FILE=M_$$ -f $LIB_DIR/db.subdb.awk $f > E_$$/$nom.fst
  Report E_$$/$nom.fst "length_filter"
  set n = `egrep '^>' E_$$/$nom.fst | wc -l`
  Notify "  length filter : $nom : $n"
  if ($n <= $Sizmin) \rm -f E_$$/$nom.fst
end

set ok = `ls E_$$ | wc -l`
if ($ok == 0) goto fin


#
# select by similarity
#

Notify "select by similarity"

foreach f (E_$$/*.fst) 
  set nom = `basename $f:r`

  Notify "  blasting $nom"
  
  makeblastdb -dbtype 'prot' -in $f >>& db.log
  blastp -db $f -query $f -outfmt 7 > $f.blast.out
  \rm -f $f.p??
  
  $AwkCmd -v COVMIN=$Covmin -v PMAX=$Pmax -v IDMIN=$Idmin \
      -f $LIB_DIR/db.blastlink.awk $f.blast.out |\
      $AwkCmd -f $LIB_DIR/db.todl.awk > G_$$
      
  ($LIB_DIR/db.cc.r G_$$ > $f.cc.txt) >>& db.log
  
  awk -v NAME=$nom -f $LIB_DIR/db.reportcc.awk $f.cc.txt >> $OutLog 
  
  $AwkCmd -f $LIB_DIR/db.selcc.awk $f.cc.txt > S_$$
  $AwkCmd -v FILE=S_$$ -f $LIB_DIR/db.subdb.awk $f > F_$$/$nom.fst

  Report F_$$/$nom.fst "similarity_filter"

  set n = `egrep '^>' F_$$/$nom.fst | wc -l`
  Notify "  blast filter : $nom : $n"
  if ($n <= $Sizmin) \rm -f F_$$/$nom.fst
  
end

set ok = `ls D_$$ | wc -l`
if ($ok == 0) goto fin

#
# annotations
#

echo -n "" > J_$$

foreach f (F_$$/*.fst) 
  $AwkCmd -f $LIB_DIR/db.annot.awk $f >> J_$$
end

awk '(NF >= 3) {print $1, $NF}' $PatFile | sort > A_$$
sort J_$$ | egrep -v '^ *$' > B_$$
join A_$$ B_$$ > F_$$/Annot.lst

#
# copy files
#

set n = `ls F_$$/* | wc -l`
Notify "copy $n files to $OutDir"

\mv -f F_$$/* $OutDir

#
# end
#

fin:
Notify "output directory : $OutDir"

\rm -r ?_$$


exit 0
