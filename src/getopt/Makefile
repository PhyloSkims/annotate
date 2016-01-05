# ---------------------------------------------------------------
# $Id: $
# ---------------------------------------------------------------
# @file: Makefile
# @desc: makefile for package exonerate
#
# @history:
# @+	    <Gloup> : Sept 15 : Adapted to ORG.Annot
#
# @note: should be processed with gnu compatible make
# @note: helixware_compatible
#
# @end:
# ---------------------------------------------------------------
#

include ../../config/auto.conf

PKG = getopt-1.1.6
MAKEOPTIONS=prefix=$(PRTPATH) 

ifeq ($(CFGPRT),i386-darwin)
	LDFLAGS:=$(LDFLAGS) -lintl
endif

include $(CFGDIR)targets/package.targ

include $(CFGDIR)targets/help.targ

