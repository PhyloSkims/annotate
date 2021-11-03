#!/usr/bin/env tcsh -f
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
# usage: go_cds.sh fasta [db_root]
#
unsetenv ORG_SOURCED

setenv ORG_HOME `dirname $0`/../../..
source $ORG_HOME/scripts/csh_init.sh

NeedArg 1

set Fasta = $Argv[1]; Shift

NeedFile $Fasta

set Genome = `basename $Fasta:r`

set DbRoot = $CDS_DATA_DIR/chlorodb

if ($#Argv > 0) then
  set DbRoot = $Argv[1]; Shift
endif

NeedDir $DbRoot/core
NeedFile $DbRoot/core/Annot.lst

NeedDir $DbRoot/models

#
# run everything into temporary place
#

set temp = `hostname`.$$.Genome.tmp
if (! -d $temp) then
  Notify "making directory $temp"
  mkdir $temp
endif

#
# find the absolute path of the fasta genome file
#
echo $Fasta | grep '^/' > /dev/null
if ( $status == 1 ) then
	set AbsGenoFile = `pwd`/$Fasta
	set DirGenoFile = `dirname $AbsGenoFile`
	set DirGenoFile = `(cd $DirGenoFile;pwd)`
	set AbsGenoFile = $DirGenoFile/`basename $AbsGenoFile`
else
	set AbsGenoFile = $Fasta
endif

pushd $temp >& /dev/null
ln -s $AbsGenoFile genome.fasta
popd >& /dev/null

set Fasta = $temp/genome.fasta


#
# pass1: run exonerate
#

foreach dir ("core" "shell" "dust")
  if (-d $DbRoot/$dir) then
    set fams = `ls $DbRoot/$dir/*.clean.fst`
    Notify "running pass1:$dir exonerate of $Genome on $DbRoot"
    foreach f ($fams)
      tcsh -f $PROG_DIR/do_exonerate.csh $Fasta $f $DbRoot/models $temp
    end
  endif
end

cp $temp/genome.cds.fasta $Genome.cds.fasta 

#
# pass2: transsplicing
#

$PROG_DIR/do_rps12.sh $Fasta > $temp/$Genome.rps12.res

#
# pass3: prokov
#

# $PROG_DIR/do_prokov.sh $Fasta $Genome.cds.fasta $temp

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
