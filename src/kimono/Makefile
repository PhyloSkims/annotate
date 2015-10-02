# ---------------------------------------------------------------
# $Id: $
# ---------------------------------------------------------------
# @file: Makefile
# @desc: makefile for lxpack
#
# @history:
# @history:
# @+	    <Gloup> : Apr 97 : Created
# @+	    <Gloup> : Mar 02 : Updated for LXxware
#
# @note: should be processed with gnu compatible make
# @note: helixware_compatible
#
# @end:
# ---------------------------------------------------------------
#
include ../../config/auto.conf

DIRS = lxpack

include ../../config/targets/propagate.targ

include ../../config/targets/help.targ

all::
	$(MAKE) ACTION=$@ _action
	test -d $(PRTDIR) || mkdir $(PRTDIR)
	test -d $(BINDIR) || mkdir $(BINDIR)
	\cp -f lxpack/ports/$(PORTNAME)/bin/* $(BINDIR)

clean::
	$(MAKE) -C lxpack portclean
