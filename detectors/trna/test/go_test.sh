#!/bin/csh -f

echo "+ [testing tRNA]"

../bin/go_trna.sh test.fst > test.bak

diff -q test.bak test.ref >& /dev/null

set stat = $status

if ($stat == 0) then
  echo '+[0;32m tRNA test Ok[m'
  \rm -r test.bak
else
  echo '+[0;31m tRNA test Failure[m'
endif

exit $stat
