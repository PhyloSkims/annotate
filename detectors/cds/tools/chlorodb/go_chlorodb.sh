#!/bin/csh -f
#
# make ChloroDB's
#
# usage: copy genbank/embl files into 'DB_DIR/download'
# usage: [create a paramter.sh file in 'DB_DIR']
# usage: go_chlorodb [DB_DIR]
#
unsetenv ORG_SOURCED

setenv ORG_HOME `dirname $0`/../../../..
source $ORG_HOME/scripts/csh_init.sh

#
# which DB to process
#

set DB_BASE = $DATA_DIR/cds/chlorodb  # default location

if ($#Argv > 0) then
  set DB_BASE = $Argv[1]; Shift
endif

set DB_BASE = `cd $DB_BASE && pwd -P`

NeedDir $DB_BASE/download

if (! -d $DB_BASE/info)  mkdir $DB_BASE/info
if (! -d $DB_BASE/fasta) mkdir $DB_BASE/fasta

cd $DB_BASE/info

#
# params
#

if (! -e $DB_BASE/parameters.sh) then
  @ n = `find $DB_BASE/download -depth 1 -type f -print | wc -l`
  @ cor_cutoff = $n / 2
  @ atg_cutoff = $n / 10
  @ dbs_cutoff = $n / 4
  if ($cor_cutoff == 0) @ cor_cutoff = 1
  if ($atg_cutoff == 0) @ atg_cutoff = 1
  if ($dbs_cutoff == 0) @ dbs_cutoff = 1
  echo "# sourced file"                          >  $DB_BASE/parameters.sh
  echo ""                                        >> $DB_BASE/parameters.sh
  echo "set CORE_NCDS_CUTOFF      = $cor_cutoff" >> $DB_BASE/parameters.sh
  echo "set CORE_START_ATG_CUTOFF = $atg_cutoff" >> $DB_BASE/parameters.sh
  echo "set CORE_START_DFT_CUTOFF = $atg_cutoff" >> $DB_BASE/parameters.sh
  echo "set CORE_START_OTH_CUTOFF = 10"          >> $DB_BASE/parameters.sh
  echo "set CORE_STOP_CUTOFF      = $cor_cutoff" >> $DB_BASE/parameters.sh
  echo "set CORE_SPLICE_CUTOFF    = $atg_cutoff" >> $DB_BASE/parameters.sh
  echo ""                                        >> $DB_BASE/parameters.sh
  echo "set SHEL_NCDS_CUTOFF      = 10"          >> $DB_BASE/parameters.sh
  echo ""                                        >> $DB_BASE/parameters.sh
  echo "set CORE_DELTA            = Inf"         >> $DB_BASE/parameters.sh
  echo "set CORE_COVMIN           = 30"          >> $DB_BASE/parameters.sh
  echo "set CORE_PMAX             = 1e-6"        >> $DB_BASE/parameters.sh
  echo "set CORE_IDMIN            = 30"          >> $DB_BASE/parameters.sh
  echo "set CORE_SIZMIN           = $cor_cutoff" >> $DB_BASE/parameters.sh
  echo ""                                        >> $DB_BASE/parameters.sh
  echo "set SHEL_DELTA            = 0.5"         >> $DB_BASE/parameters.sh
  echo "set SHEL_COVMIN           = 30"          >> $DB_BASE/parameters.sh
  echo "set SHEL_PMAX             = 1e-6"        >> $DB_BASE/parameters.sh
  echo "set SHEL_IDMIN            = 30"          >> $DB_BASE/parameters.sh
  echo "set SHEL_SIZMIN           = $dbs_cutoff" >> $DB_BASE/parameters.sh
  echo ""                                        >> $DB_BASE/parameters.sh
  echo "set DUST_DELTA            = 0.5"         >> $DB_BASE/parameters.sh
  echo "set DUST_COVMIN           = 30"          >> $DB_BASE/parameters.sh
  echo "set DUST_PMAX             = 1e-6"        >> $DB_BASE/parameters.sh
  echo "set DUST_IDMIN            = 30"          >> $DB_BASE/parameters.sh
  echo "set DUST_SIZMIN           = 10"          >> $DB_BASE/parameters.sh
  
endif

source $DB_BASE/parameters.sh

##set CMIN_COD = 0
##set FMIN_COD = 0.01

#
# temporarily uncompress
#

set ff = `find $DB_BASE/download -depth 1 -name \*.gz -print`

if ($#ff != 0) then
  Notify "uncompressing $#ff entries"
  foreach f ($ff)
    gunzip -f $f
  end
endif

#
# convert gbk/embl to fasta
#

set ff = `find $DB_BASE/download -depth 1 \( -name \*.gbk -or -name \*.embl \) -print`

Notify "convert $#ff gbk/embl entries to fasta"

foreach f ($ff)
  set nom = `basename $f:r`
  set typ = $f:e
  $AwkCmd -f $LIB_DIR/$typ.tofasta.awk $f > $DB_BASE/fasta/$nom.fst
end

#
# get gbk/embl info
#

Notify "get gbk/embl info for $#ff entries"

echo "" | awk -v HEADONLY=1 -f $LIB_DIR/gbk.info.awk > db.info.txt  # just get header

foreach f ($ff)
  set nom = `basename $f:r`
  set typ = $f:e
  $AwkCmd -f $LIB_DIR/$typ.oneliner.awk  $f |\
  $AwkCmd -f $LIB_DIR/libutil.awk -f $LIB_DIR/$typ.info.awk |\
  egrep -v '^#' >> db.info.txt
end

#
# get cds info
#

Notify "get gbk/embl cds for $#ff entries"

echo "" | awk -v HEADONLY=1 -f $LIB_DIR/gbk.cds_long.awk > db.cds.txt  # just get header

foreach f ($ff)
  set nom = `basename $f:r`
  set typ = $f:e
  $AwkCmd -f $LIB_DIR/$typ.oneliner.awk  $f |\
  $AwkCmd -v FASTA=$DB_BASE/fasta/$nom.fst -f $LIB_DIR/libutil.awk \
          -f $LIB_DIR/$typ.cds_long.awk |\
  egrep -v '^#' >> db.cds.txt
end

#
# get fasta for  prots
#

Notify "get prots"
$AwkCmd -f $LIB_DIR/libutil.awk -f $LIB_DIR/cds2fasta.awk db.cds.txt > db.prot.fst

#
# get introns
#

Notify "get gbk/embl introns for $#ff entries"

echo "" | awk -v HEADONLY=1 -f $LIB_DIR/gbk.intron.awk > db.intron.txt  # just get header

foreach f ($ff)
  set nom = `basename $f:r`
  set typ = $f:e
  $AwkCmd -f $LIB_DIR/$typ.oneliner.awk  $f |\
  $AwkCmd -v FASTA=$DB_BASE/fasta/$nom.fst -f $LIB_DIR/libutil.awk \
          -f $LIB_DIR/$typ.intron.awk |\
  egrep -v '^#' >> db.intron.txt
end

#
# make models
#

Notify "Making models"

echo -n ""                                              >  db.models.params.txt
echo "CORE_NCDS_CUTOFF <- $CORE_NCDS_CUTOFF"            >> db.models.params.txt
echo "CORE_START_ATG_CUTOFF <- $CORE_START_ATG_CUTOFF"  >> db.models.params.txt
echo "CORE_START_DFT_CUTOFF <- $CORE_START_DFT_CUTOFF"  >> db.models.params.txt
echo "CORE_START_OTH_CUTOFF <- $CORE_START_OTH_CUTOFF"  >> db.models.params.txt
echo "CORE_STOP_CUTOFF <- $CORE_STOP_CUTOFF"            >> db.models.params.txt
echo "CORE_SPLICE_CUTOFF <- $CORE_SPLICE_CUTOFF"        >> db.models.params.txt
echo "SHEL_NCDS_CUTOFF <- $SHEL_NCDS_CUTOFF"            >> db.models.params.txt

$LIB_DIR/make.models.r |& Cat

GetStatus
OnError then 
  Error 2 "model parameter too stringent"
endif

#
# add matrices
#

cp -f $PROG_DIR/matrices/* models

#
# make subDBs
#

if (-e db.core.pat.txt) then
  Notify "Making core DB (take some time... please wait)"
  $PROG_DIR/subdb/go_subdb.sh db.prot.fst db.core.pat.txt \
    $CORE_DELTA $CORE_COVMIN $CORE_PMAX $CORE_IDMIN $CORE_SIZMIN
endif

if (-e db.shell.pat.txt) then
  Notify "Making shell DB (take some time... please wait)"
  $PROG_DIR/subdb/go_subdb.sh db.prot.fst db.shell.pat.txt \
    $SHEL_DELTA $SHEL_COVMIN $SHEL_PMAX $SHEL_IDMIN $SHEL_SIZMIN
endif

if (-e db.dust.pat.txt) then
  Notify "Making dust DB (take some time... please wait)"
  $PROG_DIR/subdb/go_subdb.sh db.prot.fst db.dust.pat.txt \
    $DUST_DELTA $DUST_COVMIN $DUST_PMAX $DUST_IDMIN $DUST_SIZMIN
endif

#
# recompress entries
#

set ff = `find $DB_BASE/download -depth 1 -type f -print`

if ($#ff != 0) then
  Notify "recompressing $#ff entries"
  foreach f ($ff)
    gzip -f $f
  end
endif

# compress fasta

set ff = `find $DB_BASE/fasta -depth 1 -name \*.fst -print`

if ($#ff != 0) then
  Notify "compressing $#ff fasta entries"
  foreach f ($ff)
    gzip -f $f
  end
endif

# install everything in proper directory

foreach dir ("core" "shell" "dust")
  if (-e $DB_BASE/$dir) \rm -r $DB_BASE/$dir
  if ((-d db.$dir.pat.db) && (-e db.$dir.pat.db/Annot.lst)) then
    Notify "installing $DB_BASE/$dir"
    \mv -f db.$dir.pat.db $DB_BASE/$dir
  endif
end

if (-e $DB_BASE/models) \rm -r $DB_BASE/models
if (-d models) \mv -f models $DB_BASE

Notify "Done"
exit 0

