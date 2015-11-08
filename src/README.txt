#
# ------------------------------------------
# ORG.Annot - Organelle Annotator
# version 1.0.0 - Sept 2015
# ------------------------------------------
#
# this directory contains sources for tools used by ORG.Annot
# this is only useful if you need to recompile binaries
#
# -------------------------
# -0- Requirements
# -------------------------
#
# gnu tar     :  version >= 1.15
# gnu make    :  version >= 3.81
# gcc and g++ :  version >= 3.4.6
# gawk        :  version >= 4.0.1
#
# note: valid versions of gnu tar and make are in _unix_tools_
#
#
# -------------------------
# -1- Checking configuration files
# -------------------------
#
# issue: 
# $ config/guess_port
#
# if output is 'unknown <mach>:<sys>:<rel>' 
# then you should :
#   - add a port entry in guess_port for <mach>:<sys>:<rel>
#   - create a ports/<port>.conf configuration file
#     (the best is to start from another port file,
#      choose whatever looks closest)
#
# if output is something else (e.g. i386-darwin)
# then the configuration files are already defined for your port
#
# -------------------------
# -2- Compiling
# -------------------------
#
# just run make from 'src':
#
# $ cd src
# $ make
#
#
# if everything runs fine, you may then issue:
#
# $ make clean
#
# to cleanup unecessary files
#
# at this point all binaries should be located in ports/<PORTNAME>/bin
#
# -------------------------
# -3- Recompiling a port
# -------------------------
#
# if you need to fully recompile a port
#
# $ cd src
# $ make portclean
# $ make
# $ make clean
#




