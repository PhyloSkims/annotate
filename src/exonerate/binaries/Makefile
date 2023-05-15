# ---------------------------------------------------------------
# $Id: $
# ---------------------------------------------------------------
# @file: Makefile
# @desc: makefile for package ncbi blast / binaries
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

PKG = exonerate-2.2.0

PKGTAR = $(PKG)-$(PORTNAME).tar.gz

PKGDIR = bin-$(PORTNAME)

PROGS = exonerate

BINPATH = $(abspath $(BINDIR))

include $(CFGDIR)targets/empty.targ
include $(CFGDIR)targets/help.targ

.PHONY: extract install

all:: install

extract:
	test -d $(PKGDIR) || mkdir $(PKGDIR)
	test -d $(PKGDIR)/bin || $(TAR) zxf $(PKGTAR) -C $(PKGDIR) --strip-components 1
	
install: extract
	test -d $(BINDIR) || mkdir $(BINDIR)
	cd $(PKGDIR)/bin && \cp -f $(PROGS) $(BINPATH)

clean::
	(! test -d $(PKGDIR)) || \rm -r $(PKGDIR)

portclean::
	cd $(BINDIR) && \rm -f $(PROGS)

