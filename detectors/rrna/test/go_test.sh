#!/bin/csh -f

echo "+ [testing rRNA]"

../bin/go_rrna.sh test.fst > test.bak

diff -q test.bak test.ref >& /dev/null

set stat = $status

if ($stat == 0) then
  echo '+[0;32m rRNA test Ok[m'
  \rm -r test.bak
else
  echo '+[0;31m rRNA test Failure[m'
endif

exit $stat
