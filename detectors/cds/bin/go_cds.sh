#!/bin/bash
#
#                           Annotate CDS 
#
#========================================================================================
#
#  Annotate CDS
#
#  go_cds.sh <FASTAFILE> [DBROOT]
#
#		- <FASTAFILE> : The fasta file containing the genome to annotate
#       - [DBROOT]    : optionnal argument allowing to specify database directory
#
#  Results are printed to the standard output
#
#========================================================================================
# usage: go_cds.sh fasta [db_root]
#

# -- CAUTION -- Works as long than the script 
#               is not called through a symlink
THIS_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${THIS_DIR}/../../../scripts/bash_init.sh"

needarg 1

Fasta=$1; shift

needfile $Fasta

# Genome names is set from the base
# name of the genome file without its extension
Genome=$(basename ${Fasta%.*})

# DbRoot is set to its default values except 
# if the second argument precise another DbRoot

DbRoot="$CDS_DATA_DIR/chlorodb"

if (( $# > 0)) ; then
  DbRoot="$1"; Shift
fi

needdir $DbRoot
needdir $DbRoot/core
needfile $DbRoot/core/Annot.lst
needdir $DbRoot/models

assignundef cdsdetection_pass1 yes
assignundef cdsdetection_pass2 yes

temp=$(mktempdir $(hostname))

AbsGenoFile=$(getAbsolutePath $Fasta)
pushd $temp >& /dev/null
ln -s $AbsGenoFile genome.fasta
popd >& /dev/null

Fasta="$temp/genome.fasta"

#
# pass1: run exonerate
#

if [[ "$cdsdetection_pass1" == "yes" ]] ; then
    for dir in "core" "shell" "dust" ; do
    if [[ -d $DbRoot/$dir ]] ; then
        fams=$(ls $DbRoot/$dir/*.clean.fst)
        loginfo "running pass1:$dir exonerate of $Genome on $DbRoot"
        for f in $fams ; do
        tcsh -f $PROG_DIR/do_exonerate.csh $Fasta $f $DbRoot/models $temp
        done
    fi
    done

    cp $temp/genome.cds.fasta $Genome.cds.fasta 

fi


#
# pass2: transsplicing
#

if [[ "$cdsdetection_pass2" == "yes" ]] ; then
    $PROG_DIR/do_rps12.sh $Fasta $temp
fi

#
# pass3: prokov
#

# $PROG_DIR/do_prokov.sh $Fasta $Genome.cds.fasta $temp

#
# end : output on stdout
#

cat $temp/*.res

# cleanup everything

assignundef TMP_CLEANUP 1

if (( $TMP_CLEANUP != 0 )) ; then
  loginfo "  cleanup $temp"
  (\rm -r $temp) >& /dev/null
fi

exit 0
