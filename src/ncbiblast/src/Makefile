# ---------------------------------------------------------------
# $Id: $
# ---------------------------------------------------------------
# @file: Makefile
# @desc: makefile for package ncbi blast
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

include ../../../config/auto.conf

PKG = ncbi-blast-2.2.31

PKGPATH := build.$(PORTNAME)

# catchup : test and clean rules
# ncbi blast don't know how to do this...

test:: PKGDIR = _sink_
clean:: PKGDIR = _sink_
portclean:: PKGDIR = _sink_

include $(CFGDIR)targets/package.targ

include $(CFGDIR)targets/help.targ

portclean::
	(! test -d $(PKGPATH)) || \rm -r $(PKGPATH)
