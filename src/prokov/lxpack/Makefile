# ---------------------------------------------------------------
# $Id: $
# ---------------------------------------------------------------
# @file: Makefile
# @desc: makefile for prokov
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
include ../../../config/auto.conf

DIRS = src \
       tests

include ../../../config/targets/propagate.targ

include ../../../config/targets/help.targ

portclean::
	$(MAKE) ACTION=$@ _action
	(! test -d $(PRTDIR)) || \rm -r $(PRTDIR)
