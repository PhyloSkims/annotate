#!/bin/csh -f

echo "+ [testing Normalize]"

../bin/go_normalize.sh test.fst > test.bak

diff -q test.bak test.ref >& /dev/null

set stat = $status

if ($stat == 0) then
  echo '+[0;32m Normalize test Ok[m'
  \rm -r test.bak
else
  echo '+[0;31m Normalize test Failure[m'
endif

exit $stat
