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

PKG = ITSx_1.0.11
PKGTAR ?= $(PKG).tgz

PKGDIR ?= build.$(PORTNAME)

PRTPATH = $(abspath $(PRTDIR))
DATADIR = $(CFGDIR)../data
DATAITS = $(DATADIR)/its

#
# Rules
#

.PHONY: all clean test portclean pkg pkg.expand pkg.install

all:: pkg

pkg.expand::
	test -d $(PKGDIR) || mkdir $(PKGDIR)
	$(TAR) zxf $(PKGTAR) -C $(PKGDIR) --strip-components 1

pkg.install:: pkg.expand
	@mkdir -p $(BINDIR)
	@cp $(PKGDIR)/ITSx $(BINDIR)
	@mkdir -p $(DATAITS)
	@cp -r $(PKGDIR)/ITSx_db $(DATAITS)
	@echo "+++++++++++ package $(PKG) done"

pkg:: pkg.install

test::
	echo No test available
	
clean::
	echo Done

portclean::
	(! test -d $(PKGDIR)) || \rm -r $(PKGDIR)


include $(CFGDIR)targets/help.targ

