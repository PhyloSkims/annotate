
$Id: README.txt 1825 2013-02-26 09:39:47Z viari $

This directory contains Makefile machine specific configuration files
(and default targets to help you writing Makefile's)

These headers should be used with GNU make or compatible

#
# portname
#

To check your port, issue :

  ./guess_port
  
  if output is 'unknown <mach>:<sys>:<rel>' then you should :
    - add a port entry in guess_port for <mach>:<sys>:<rel>
    - create a ports/<port>.conf configuration file
      (the best is to start from another port file,
       choose whatever looks closest)

#
# configuration flags
#

auto.conf : the main configuration file :
             - determine the machine port thru 'guess_port' shell
             - include 'default.conf' file
             - include the machine specific 'ports/<port>.conf' file

default.conf : default configuration (included by 'auto.conf')

ports/<port>.conf  : machine specific configuration (included by 'auto.conf')

#
# utility targets
#

targets/help.targ : target for standard help

targets/propagate.targ : target for propagating targets to subdirectories

targets/package.targ : default targets for standard package with 'configure'

targets/empty.targ : default empty targets (defined as double colon rules)

targets/lxbin.targ : default make targets for standard lx binary (without libraries)

targets/debug.targ : target to print debug information (for dev.)

