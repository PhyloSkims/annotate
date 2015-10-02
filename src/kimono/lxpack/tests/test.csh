#!/bin/csh -f

set path = (. ../src $path); rehash

set out = "tst"
#set out = "ref"  # for recording

echo borbu1.0.$out
kimono borbu1.fst | awk '{print $1, $2, $3, $4, $5}' > borbu1.0.$out
diff borbu1.0.$out ref/borbu1.0.ref
if (($out != "ref") && ($status != 0)) exit $status

echo borbu1.1.$out
kimono -i borbu1.lst -p 1 borbu1.fst  | awk '{print $1, $2, $3, $4, $5}' > borbu1.1.$out
diff borbu1.1.$out ref/borbu1.1.ref
if (($out != "ref") && ($status != 0)) exit $status

echo borbu1.2.$out
kimono -i borbu1.lst -p 2 borbu1.fst  | awk '{print $1, $2, $3, $4, $5}' > borbu1.2.$out
diff borbu1.2.$out ref/borbu1.2.ref
if (($out != "ref") && ($status != 0)) exit $status

echo borbu1.3.$out
kimono -i borbu1.lst -p 123 borbu1.fst  | awk '{print $1, $2, $3, $4, $5}' > borbu1.3.$out
diff borbu1.3.$out ref/borbu1.3.ref
if (($out != "ref") && ($status != 0)) exit $status

echo borbu1.4.$out
kimono -x borbu1.lst borbu1.fst  | awk '{print $1, $2, $3, $4, $5}' > borbu1.4.$out
diff borbu1.4.$out ref/borbu1.4.ref
if (($out != "ref") && ($status != 0)) exit $status

echo borbu1.5.$out
kimono -a gc -c cumul borbu1.fst  | awk '{print $1, $2, $3, $4, $5}' > borbu1.5.$out
diff borbu1.5.$out ref/borbu1.5.ref
if (($out != "ref") && ($status != 0)) exit $status

echo borbu1.6.$out
kimono -c codon -e 0 -i borbu1.lst borbu1.fst  | awk '{print $1, $2, $3, $4, $5}' > borbu1.6.$out
diff borbu1.6.$out ref/borbu1.6.ref
if (($out != "ref") && ($status != 0)) exit $status

echo borbu1.7.$out
kimono -c genes -i borbu1.lst borbu1.fst  | awk '{print $1, $2, $3, $4, $5}' > borbu1.7.$out
diff borbu1.7.$out ref/borbu1.7.ref
if (($out != "ref") && ($status != 0)) exit $status

echo "+ all test ok"

if ($out != "ref") \rm -f *.$out

exit 0
