# ---------------------------------------------------------------
# $Id: Makefile 1825 2013-02-26 09:39:47Z viari $
# ---------------------------------------------------------------
# @file: Makefile
# @desc: makefile for src
#
# @history:
# @+	    <Gloup> : Mar 92 : first draft
# @+	    <Gloup> : Jul 93 : revised version
# @+	    <Gloup> : Feb 01 : Adapted to LXware
#
# @note: should be processed with gnu compatible make
# @note: helixware_compatible
#
# @end:
# ---------------------------------------------------------------
#
include ../config/auto.conf

DIRS = aragorn   \
       hmmer3    \
       exonerate \
       kimono    \
       ncbiblast \
       prokov    \
       repseek   \
       sequtils  \
       sumaclust \
       sumatra   \
       getopt    \
       ITSx
#       muscle    \
#       clustalo  \
#       cd-hit	 \

include ../config/targets/propagate.targ

include ../config/targets/help.targ

portclean::
	$(MAKE) ACTION=$@ _action
	(! test -d $(PRTDIR)) || \rm -r $(PRTDIR)
