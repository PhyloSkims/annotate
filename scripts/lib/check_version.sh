#!/bin/csh -f
#
# check (gnu) software version
#
# usage: check_version.sh cmd version
#
#
set HERE = `dirname $0`
set ORG_HOME = $HERE/../..
set ORG_PORTNAME = `$ORG_HOME/config/guess_port`

set path = ($path $ORG_HOME/ports/$ORG_PORTNAME/bin)

if ($#argv != 2) then
  egrep "^# *usage:" $0 | sed -e 's/^# *//1'
  exit 1
endif

set name = `echo "$1" | awk '{print $1}'`

echo -n "+ checking $name"

set ok = `which $name`

if ($status != 0) then
  echo ""
  echo "* $name : command not found"
  echo "* please consider installing $name from src/_unix_tools_"
  exit 1
endif

set num = `eval "$1" | tr -d -C '[0-9].'`

echo -n " version $num"

awk -v NUM=$num -v REF=$2 -f $HERE/version.awk

if ($status != 0) then
  echo ""
  echo "[0;31m* $name version $num (should be >= $2)[m"
  echo "[0;31m* please consider installing $name from src/_unix_tools_[m"
  exit 2
endif

echo " (>= $2) [0;32mOK[m"

exit 0
              



