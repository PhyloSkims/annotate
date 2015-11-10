#!/bin/csh -f
#
# compare CDS annotation in reference file to predicted file
# annotation file are in Genbank/Embl format
#
# usage: go_compare reference predicted
#
# output on stdout
#
unsetenv ORG_SOURCED

setenv ORG_HOME `dirname $0`/../../../..
source $ORG_HOME/scripts/csh_init.sh

NeedArg 2

set RefFile = $Argv[1]
set PrdFile = $Argv[2]

NeedFile $RefFile
NeedFile $PrdFile

set RefType = $RefFile:e
set PrdType = $PrdFile:e

if ((! -e $LIB_DIR/$RefType.oneliner.awk) || (! -e $LIB_DIR/$PrdType.oneliner.awk)) then
  Error 1 "file extension should be 'gbk' or 'embl'"
endif

#
# parse ref and prediction
#

Notify "get genome info from $RefFile"

$AwkCmd -f $LIB_DIR/$RefType.oneliner.awk $RefFile |\
$AwkCmd -f $LIB_DIR/libutil.awk -f $LIB_DIR/$RefType.cds.awk > R_$$

Notify "get prediction info from $PrdFile"

$AwkCmd -f $LIB_DIR/$PrdType.oneliner.awk $PrdFile |\
$AwkCmd -f $LIB_DIR/libutil.awk -f $LIB_DIR/$PrdType.cds.awk > P_$$

#
# compare
#

Notify "compare bank to predictions"

$AwkCmd -f $LIB_DIR/libnws.awk      \
        -f $LIB_DIR/compareCds.awk  \
        R_$$ P_$$ > S_$$

# base statistics

egrep "^MATCH" S_$$ | tr '.' ' ' | awk '{print $5}' |\
sort | uniq -c | sort -nr | awk '{print "#",$0}' > U_$$

# add chlorodb/core statistics

if (-d $DATA_DIR/cds/chlorodb/core) then

  ls $DATA_DIR/cds/chlorodb/core/*.fst |\
  sed -e 's@^.*core/@@1' | sed -e 's/.fst$//g' |\
  sort > C_$$

  egrep "^MATCH" S_$$ | grep "MISSED" | awk '{print $2}' | sort | uniq > D_$$

  join D_$$ C_$$ > E_$$
  @ nc = `cat C_$$ | wc -l`
  @ mt = `cat D_$$ | wc -l`
  @ mc = `cat E_$$ | wc -l`
  @ mn = $mt - $mc
  set LC = `cat E_$$`

  echo "#"                                   >> U_$$
  echo "# $mc MISSED in ChloroDB-Core ($LC)" >> U_$$
  echo "# $mn MISSED not in ChloroDB-Core"   >> U_$$
  echo "#"                                   >> U_$$
  echo ""                                    >> U_$$
endif

cat S_$$ >> U_$$

cat U_$$

#
# end
#

(\rm -f ?_$$) >> /dev/null

Exit 0
