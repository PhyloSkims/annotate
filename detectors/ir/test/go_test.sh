#!/bin/csh -f

echo "+ [testing IR]"

../bin/go_ir.sh test.fst > test.bak

diff -q test.bak test.ref >& /dev/null

set stat = $status

if ($stat == 0) then
  echo '+[0;32m IR test Ok[m'
  \rm -r test.bak
else
  echo '+[0;31m IR test Failure[m'
endif

exit $stat
