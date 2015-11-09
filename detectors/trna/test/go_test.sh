#!/bin/csh -f

setenv ORG_HOME `dirname $0`/../../..
source $ORG_HOME/scripts/csh_init.sh

echo "+ testing tRNA"

`dirname $0`/../bin/go_trna.sh test.fst > test.bak

diff -q test.bak test.ref >& /dev/null

set stat = $status

if ($stat == 0) then
  echo "+ $VTC[3]tRNA test Ok$VTC[1]"
  \rm -r test.bak
else
  echo "* $VTC[2]tRNA test Failure$VTC[1]"
endif

exit $stat
