#!/bin/csh -f

setenv ORG_HOME `dirname $0`/../../..
source $ORG_HOME/scripts/csh_init.sh

echo "+ testing CDS"

setenv TMP_CLEANUP 0
setenv PASS1_SPEEDUP  1
setenv PASS1_SLOWDOWN 0
setenv PASS1_BLASTX_FILTER_NBMAX 10

`dirname $0`/../bin/go_cds.sh test.fst > test.bak

diff -q test.bak test.ref >& /dev/null

set stat = $status

if ($stat == 0) then
  echo "+ $VTC[3]CDS test Ok$VTC[1]"
  \rm -r test.bak test.tmp
else
  echo "* $VTC[2]CDS test Failure$VTC[1]"
endif

exit $stat
