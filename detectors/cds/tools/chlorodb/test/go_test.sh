#!/bin/csh -f

setenv ORG_HOME `dirname $0`/../../../../..
source $ORG_HOME/scripts/csh_init.sh

echo "+ testing go_chlorodb.sh"

\cp -r test.db TMP
`dirname $0`/../go_chlorodb.sh TMP


find TMP -name \*.fst -print | sort > test.bak

diff -q test.bak test.ref >& /dev/null

set stat = $status

if ($stat == 0) then
  echo "+ $VTC[3]test Ok$VTC[1]"
  \rm -r TMP test.bak
else
  echo "* $VTC[2]test Failure$VTC[1]"
endif

exit $stat
