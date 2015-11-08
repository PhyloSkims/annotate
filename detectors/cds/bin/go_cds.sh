#!/bin/csh -f
#
#                           Annotate CDS 
#
#========================================================================================
#
#  Annotate CDS
#
#  go_cds.sh <FASTAFILE>
#
#		- <FASTAFILE> : The fasta file containing the genome to annotate
#
#  Results are printed to the standard output
#
#========================================================================================
# usage: go_cds.sh fasta
#
unsetenv ORG_SOURCED

setenv ORG_HOME `dirname $0`/../../..
source $ORG_HOME/scripts/csh_init.sh

NeedArg 1

set Fasta = $Argv[1]

NeedFile $Fasta

set Genome = `basename $Fasta:r`

NeedFile $CDS_DATA_DIR/chlorodb/core

#
# run everything into temporary place
#

set temp = $Genome.tmp
if (! -d $temp) then
  Notify "making directory $temp"
  mkdir $temp
endif

#
# pass1: run exonerate
#

set fams = `ls $CDS_DATA_DIR/chlorodb/core/*.fst`

Notify "running pass1: exonerate of $Genome"

foreach f ($fams)
  set prot = `basename $f:r`
  $PROG_DIR/go_pass1.sh $Fasta $prot $temp
end

#
# pass2: transsplicing
#

#
# pass3: prokov
#

#
# end : output on stdout
#

cat $temp/*.res

# cleanup everything

AssignUndef TMP_CLEANUP 1

if ($TMP_CLEANUP != 0) then
  Notify "  cleanup $temp"
  (\rm -r $temp) >& /dev/null
endif

Exit 0
