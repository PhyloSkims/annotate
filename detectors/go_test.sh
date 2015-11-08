#!/bin/csh -f
#

set dirs = ("normalize" "ir" "rrna" "trna" "cds")

echo -n "" > test.log

@ nerr = 0

foreach d ($dirs)
  (cd $d/test && go_test.sh) |& tee -a test.log |& egrep '^\+|\*'
  @ nerr += $status
end

exit $nerr


