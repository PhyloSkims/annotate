#!/bin/csh -f

setenv ORG_HOME `dirname $0`/../../..
source $ORG_HOME/scripts/csh_init.sh

echo "+ testing IR"

`dirname $0`/../bin/go_ir.sh test.fst > test.bak

diff -q test.bak test.ref >& /dev/null

set stat = $status

if ($stat == 0) then
  echo "+ $VTC[3]IR test Ok$VTC[1]"
  \rm -r test.bak
else
  echo "* $VTC[2]IR test Failure$VTC[1]"
endif

exit $stat
