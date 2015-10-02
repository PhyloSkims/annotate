#!/bin/csh -f
#
# prokov tests
#

#
# learn
#

set test = "learning - 1"
echo "test: $test"

../src/prokov_learn -a -k 4 -o test.ascii.mat.bak bacil.learn

diff test.ascii.mat.bak test.ascii.mat.ref
if ($status != 0) goto failure

set test = "learning - 2"
echo "test: $test"

../src/prokov_learn -k 4 -o test.bin.mat.bak bacil.learn
diff test.bin.mat.bak test.bin.mat.ref
if ($status != 0) goto failure

set test = "learning - 3"
echo "test: $test"

../src/prokov_learn -k 4 -K 2 -o test.bin.hol.mat.bak bacil.learn bacil.holes
diff test.bin.hol.mat.bak test.bin.hol.mat.ref
if ($status != 0) goto failure


#
# score
#

set test = "scoring - 1"
echo "test: $test"

../src/prokov_score -m test.bin.mat.ref bacil.check > ! test.score.bak
./diff.csh test.score.bak test.score.ref
if ($status != 0) goto failure


#
# curve
#

set test = "curve - 1"
echo "test: $test"

../src/prokov_curve -m test.bin.mat.ref bacil.test > ! test.curve.bak
./diff.csh test.curve.bak test.curve.ref
if ($status != 0) goto failure

#
# curve frame tests
#

set test = "curve - 2"
echo -n "test: $test "

foreach f (*S*.fasta)
   echo -n "$f:r "
   ../src/prokov_curve -m test.bin.mat.ref $f > ! $f:r.curve.bak
  ./diff.csh $f:r.curve.bak $f:r.curve.ref
  if ($status != 0) goto failure
end

echo ""

#
# cds
#

set test = "cds - 1"
echo "test: $test"

../src/prokov_cds -m test.bin.mat.ref bacil.test > ! test.cds.bak
./diff.csh test.cds.bak test.cds.ref
if ($status != 0) goto failure

#
# test environnement
#

set test = "envir - 1"
echo "test: $test"

mkdir TEST_MAT
cp test.bin.mat.ref TEST_MAT/toto
setenv PROKOV_MATDIR ./TEST_MAT
../src/prokov_cds -m toto bacil.test > ! test.cds.bak
if ($status != 0) goto failure

#
# Success
#

echo " "
echo " All tests were succesfull "
echo " "

set stat = 0

goto fin

#
# Failure
#

failure:

echo "--------------------------------------------------"
echo " ! FAILURE at Test $test"
echo "--------------------------------------------------"

set stat = 1

goto fin

#
# fin
#

fin:

if ($stat == 0) \rm -f *.bak
\rm -r TEST_MAT

exit $stat
