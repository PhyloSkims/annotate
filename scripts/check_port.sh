#!/bin/csh -f
#
# check unix tools version and port installation
#

set HERE = `dirname $0`
set ORG_HOME = $HERE/..
set LIB = $HERE/lib
set ORG_PORTNAME = `$ORG_HOME/config/guess_port`

set path = ($path $ORG_HOME/ports/$ORG_PORTNAME/bin)

@ nerr = 0

# gnu gawk version >= 4.0.1

set ref = 4.0.1
$LIB/check_version.sh 'gawk --version | head -1 | awk '"'"'{print $3}'"'" $ref
@ nerr += $status

# gnu make version >= 3.81

set ref = 3.81
$LIB/check_version.sh 'make --version | head -1 | awk '"'"'{print $3}'"'" $ref
@ nerr += $status

# (gnu) tar version >= 1.15

set ref = 1.15
$LIB/check_version.sh 'tar --version | head -1 | awk '"'"'{print $NF}'"'" $ref
@ nerr += $status

# gcc version >= 3.4.6

set ref = 3.4.6
$LIB/check_version.sh 'gcc --version | head -1 | awk '"'"'{print $NF}'"'" $ref
@ nerr += $status

# check errors

if ($nerr > 0) then
  echo "[0;31m* invalid unix tools version(s)[m"
  exit 2
endif

# check binaries architecture

set prog = $ORG_HOME/ports/$ORG_PORTNAME/bin/blastp

if (-e $prog) then 
  echo -n "+ checking port compilation"
  ($prog -h) >& /dev/null
  set Stat = $status
  if ($Stat != 0) then
    echo "[0;31m Wrong architecture[m"
    exit 2
  else
    echo "[0;32m OK[m"
  endif
else
  echo "+ port not yet compiled"
endif

exit 0

