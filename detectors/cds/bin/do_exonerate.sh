#!/bin/csh -f
#
#                           Annotate CDS - Exonerate
#
#========================================================================================
#
#  Annotate CDS using exonerate
#
#  do_exonerate.sh <FASTAGENOM> <FASTAPROT> [<OUTDIR>]
#
#		- <FASTAGENOM> : The fasta file containing the genome to annotate
#		- <FASTAPROT>  : The fasta file containing the protein family
#
#  Results are in file : `basename <FASTAGENOM>:r`.`basename <FASTAPROT>:r`.res 
#
#========================================================================================
#
# usage: do_exonerate.sh dna.fasta prot.fasta [outdir]
#
unsetenv ORG_SOURCED

setenv ORG_HOME `dirname $0`/../../..
source $ORG_HOME/scripts/csh_init.sh

set PARAMS_DIR = $LIB_DIR/../params
set MODELS_DIR = $LIB_DIR/../models

alias Override 'if (-e \!:2) set \!:1 = \!:2'

NeedArg 2

set GenoFile = $Argv[1]
set GenoName = `basename $GenoFile:r`

set ProtFile = $Argv[2]
set ProtDir  = `dirname $ProtFile`
set ProtName = `basename $ProtFile:r`

NeedFile $GenoFile
NeedFile $ProtFile
NeedFile $ProtDir/Annot.lst

set OutDir = .
if ($#Argv >= 3) set OutDir = $3
if (! -d $OutDir) mkdir $OutDir

#
# general parameters
#

source $PARAMS_DIR/default

#
# family specific parameters
#

if (-e $PARAMS_DIR/$ProtName) then
  source $PARAMS_DIR/$ProtName
endif

#
# start/stop/splice models
#

if ($?STARTMODEL == 0) then
  set STARTMODEL = $MODELS_DIR/start.default.frq
  Override STARTMODEL $MODELS_DIR/start.$ProtName.frq
endif

if ($?STOPMODEL == 0) then
  set STOPMODEL = $MODELS_DIR/stop.default.frq
  Override STOPMODEL $MODELS_DIR/stop.$ProtName.frq
endif

if ($?SPLICE3MODEL == 0) then
  set SPLICE3MODEL = $MODELS_DIR/splice3.default.frq
  Override SPLICE3MODEL $MODELS_DIR/splice3.$ProtName.frq
endif

if ($?SPLICE5MODEL == 0) then
  set SPLICE5MODEL = $MODELS_DIR/splice5.default.frq
  Override SPLICE5MODEL $MODELS_DIR/splice5.$ProtName.frq
endif

#
# out files prefix
#

set base = $OutDir/$GenoName.$ProtName

#
# skip exonerate calculations if already done
#

if (-e $base.exo.raw) then
  Notify "  file $base.exo.raw found <exonerate skipped>"
  goto parse
endif

#
# speedup exonerate 
#

if ($PASS1_SPEEDUP != 0) then

  $PROG_DIR/do_filterbx.sh $GenoFile $ProtFile  \
            $PASS1_BLASTX_FILTER_IDMIN          \
            $PASS1_BLASTX_FILTER_NBMIN          \
            $PASS1_BLASTX_FILTER_NBMAX > D_$$

  set n = `egrep "^>" D_$$ | wc -l`
  if ($n > 0) then
    Notify "  $n sequences kept"
    set DbFile = D_$$
  else
    Notify "  no sequence match"
    if ($PASS1_SLOWDOWN != 0) then
      Notify "  reverting to original $ProtName"
      set DbFile = $ProtFile
    else
      echo "" > $base.exo.raw
      goto parse
    endif
  endif
else
  set DbFile = $ProtFile
endif

#
# run exonerate
#

Notify "  running exonerate of $GenoName on $ProtName"
exonerate --model protein2genome            \
          --percent $PASS1_PERCENT          \
          --showalignment TRUE              \
          --showvulgar TRUE                 \
          --showtargetgff TRUE              \
          --geneticcode $PASS1_GENETIC_CODE \
          --minintron $PASS1_MIN_INTRON     \
          --maxintron $PASS1_MAX_INTRON     \
          --bestn $PASS1_BESTN              \
          --frameshift $PASS1_FRAMESHIFT    \
          --proteinsubmat $PASS1_SUBMAT     \
          --splice3 $SPLICE3MODEL           \
          --splice5 $SPLICE5MODEL           \
          $DbFile $GenoFile > $base.exo.raw
CheckAbort 20 "exonerate failure"

#
# extract best clusters
#
parse:

$AwkCmd -v MAX_SPAN=$PASS1_MAX_SPAN         \
        -v ALLOW_STOP=$PASS1_ALLOW_STOP     \
        -v EXCLUDE=$GenoName                \
        -f $LIB_DIR/bestclust.awk $base.exo.raw > $base.exo.best

#
# get annotations
#

egrep "^$ProtName " $ProtDir/Annot.lst | awk '{print "c annot", $0}' > T_$$

#
# extend start/stop
#

$AwkCmd -f $LIB_DIR/libutil.awk -f $LIB_DIR/extend.awk             \
                                -v FASTA=$GenoFile                 \
                                -v START_MODEL=$STARTMODEL         \
                                -v STOP_MODEL=$STOPMODEL           \
                                -v START_WALK=$PASS1_START_WALK    \
                                -v STOP_WALK=$PASS1_STOP_WALK      \
                                    $base.exo.best >> T_$$
#
# translate
#

$AwkCmd -v FASTA=$GenoFile -f $LIB_DIR/libutil.awk \
        -f $LIB_DIR/translate.awk T_$$ > $base.iff

#
# convert to embl
#

$AwkCmd -f $LIB_DIR/toEmbl.awk $base.iff |\
  $AwkCmd -f $LIB_DIR/cutline.awk > $base.res

#
# end
#

Notify "  output file: $base.res"

(\rm -f ?_$$) >> /dev/null

Exit 0
