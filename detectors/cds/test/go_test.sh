#!/bin/csh -f

echo "+ [testing CDS]"

setenv PASS1_SPEEDUP  1
setenv PASS1_SLOWDOWN 0
setenv PASS1_BLASTX_FILTER_NBMAX 50

../bin/go_cds.sh test.fst > test.bak

diff -q test.bak test.ref >& /dev/null

set stat = $status

if ($stat == 0) then
  echo '+[0;32m CDS test Ok[m'
  \rm -r test.bak
else
  echo '*[0;32m CDS test Failure[m'
endif

exit $stat
