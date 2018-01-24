#!/usr/bin/env tcsh -f
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
# usage: do_exonerate.sh dna.fasta prot.fasta [model_dir [out_dir]]
#
unsetenv ORG_SOURCED

setenv ORG_HOME `dirname $0`/../../..
source $ORG_HOME/scripts/csh_init.sh

alias Override 'if (-e \!:2) set \!:1 = \!:2'

NeedArg 2

set GenoFile = $Argv[1]; Shift
set GenoName = `basename $GenoFile:r`

set ProtFile = $Argv[1]; Shift
set ProtDir  = `dirname $ProtFile`
set ProtName = `basename $ProtFile | $AwkCmd -F'.' '{print $1}'`
set ProtType = `basename $ProtDir`

NeedFile $GenoFile
NeedFile $ProtFile
NeedFile $ProtDir/Annot.lst

set ModelsDir = $PROG_DIR/../models
if ($#Argv > 0) then
  set ModelsDir = $Argv[1]; Shift
  Notify "  exonerate models : $ModelsDir"
else
  Warning "  using default exonerate models : $ModelsDir"
endif

NeedFile $ModelsDir/start.default.frq
NeedFile $ModelsDir/stop.default.frq
NeedFile $ModelsDir/splice3.default.frq
NeedFile $ModelsDir/splice5.default.frq

set OutDir = .
if ($#Argv > 0) then 
  set OutDir = $Argv[1]; Shift
endif

if (! -d $OutDir) mkdir $OutDir

set ParamsDir = $PROG_DIR/../params

NeedFile $ParamsDir/default

#
# general parameters
#

source $ParamsDir/default

#
# family specific parameters
#

if (-e $ParamsDir/$ProtName) then
  Notify "  override parameters with $ParamsDir/$ProtName"
  source $ParamsDir/$ProtName
endif

#
# start/stop/splice models
#

if ($?STARTMODEL == 0) then
  set STARTMODEL = $ModelsDir/start.default.frq
  Override STARTMODEL $ModelsDir/start.$ProtName.frq
endif

if ($?STOPMODEL == 0) then
  set STOPMODEL = $ModelsDir/stop.default.frq
  Override STOPMODEL $ModelsDir/stop.$ProtName.frq
endif

if ($?SPLICE3MODEL == 0) then
  set SPLICE3MODEL = $ModelsDir/splice3.default.frq
  Override SPLICE3MODEL $ModelsDir/splice3.$ProtName.frq
endif

if ($?SPLICE5MODEL == 0) then
  set SPLICE5MODEL = $ModelsDir/splice5.default.frq
  Override SPLICE5MODEL $ModelsDir/splice5.$ProtName.frq
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

  tcsh -f $PROG_DIR/do_filterbx.csh $GenoFile $ProtFile  \
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
exonerate \
	--model protein2genome            \
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
    --refine region                   \
    --refineboundary 5000             \
    --singlepass FALSE                \
    --dpmemory 1024                   \
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

egrep "^$ProtName " $ProtDir/Annot.lst | $AwkCmd '{print "c annot", $0}' > T_$$

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

echo "c pass pass1 $ProtType" > $base.iff

$AwkCmd -v FASTA=$GenoFile -f $LIB_DIR/libutil.awk \
        -f $LIB_DIR/translate.awk T_$$ >> $base.iff
        
#
# extract CDS
#

$AwkCmd -v FASTA=$GenoFile -f $LIB_DIR/libutil.awk \
        -f $LIB_DIR/cds.awk T_$$ >> $OutDir/$GenoName.cds.fasta


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
