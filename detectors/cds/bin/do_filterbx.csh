#!/usr/bin/env tcsh -f
#
# filter a DB thru BlastX
# usually to speedup further DB search 
#
# output on stdout
#
# usage: do_filterbx.sh dna.fasta prot.fasta [idmin nbmin nbmax]
#
unsetenv ORG_SOURCED

setenv ORG_HOME `dirname $0`/../../..
source $ORG_HOME/scripts/csh_init.sh

set ParamsDir = $PROG_DIR/../params
set ModelsDir = $PROG_DIR/../models

NeedArg 2

set GenoFile = $Argv[1]; Shift
set ProtFile = $Argv[1]; Shift

NeedFile $GenoFile
NeedFile $ProtFile
NeedFile $ParamsDir/default

#
# general parameters
#

source $ParamsDir/default


set IDMIN = 70
set NBMIN = 50
set NBMAX = 200

if ($#Argv >= 1) set IDMIN = $Argv[1]
if ($#Argv >= 2) set NBMIN = $Argv[2]
if ($#Argv >= 3) set NBMAX = $Argv[3]

#
# format ProtFile
#

if (! -e $ProtFile.phr) then
  Notify "  formatting $ProtFile"
  (makeblastdb -dbtype prot -in $ProtFile) >& /dev/null
  CheckAbort 10 "makeblastdb failure"
endif

#
# blastx
#

Notify "  blasting $ProtFile"
blastx 	-query $GenoFile -db $ProtFile -outfmt 7   \
		-query_gencode $PASS1_BLASTX_FILTER_GCODE  \
		-max_intron_length $PASS1_MAX_INTRON      |\
  $AwkCmd -v IDMIN=$IDMIN                          \
          -v NBMIN=$NBMIN                          \
          -v NBMAX=$NBMAX                          \
          -f $LIB_DIR/filterbx.awk > T_$$
#
# extract subdb
#

$AwkCmd -v FILE=T_$$ -f $LIB_DIR/subdb.awk $ProtFile

#
# end
#

(\rm -f ?_$$) >> /dev/null

Exit 0
