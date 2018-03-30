# ---------------------------------------------------------------
# $Id: $
# ---------------------------------------------------------------
# @file: Makefile
# @desc: makefile for aragorn
#
# @history:
# @history:
# @+	    <Gloup> : Sep 12: Created
#
# @note: should be processed with gnu compatible make
# @note: helixware_compatible
#
# @end:
# ---------------------------------------------------------------
#
include ../../config/auto.conf

OSRC = aragorn1.2.38.c

PROGS = aragorn

include ../../config/targets/lxbin.targ
include ../../config/targets/help.targ

INCDIR = .

aragorn: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
