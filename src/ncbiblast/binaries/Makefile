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

PKG = ncbi-blast-2.2.31

PKGTAR = $(PKG).$(PORTNAME).tgz

PKGDIR = ${PKG}.$(PORTNAME)

PROGS = blast_formatter makeblastdb \
        blastn blastp blastx

BINPATH = $(abspath $(BINDIR))

include $(CFGDIR)targets/empty.targ
include $(CFGDIR)targets/help.targ

.PHONY: extract install

all:: install
	
install: extract
	cd $(PKGDIR)/bin && \cp -f $(PROGS) $(BINPATH)

clean::
	(! test -d $(PKGDIR)) || \rm -r $(PKGDIR)

portclean::
	cd $(BINDIR) && \rm -f $(PROGS)

