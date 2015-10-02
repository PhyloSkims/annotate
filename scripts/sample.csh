#!/bin/csh
#
# a sample script to show how to retrieve proper path for port
#

set ORG_HOME = `dirname $0`/..

set PORTNAME = `$ORG_HOME/config/guess_port`

set path = ($path $ORG_HOME/ports/$PORTNAME/bin)

#

which aragorn

exit 0



