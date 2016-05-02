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

HMMPRESS= $(BINDIR)/hmmpress

HMMDIR  = $(PKGDIR)/ITSx_db/HMMs
HMMS    = $(wildcard $(HMMDIR)/*.hmm)
HMMP    = $(patsubst %.hmm,%.hmm.h3p,$(HMMS))
HMMI    = $(patsubst %.hmm,%.hmm.h3i,$(HMMS))
HMMM    = $(patsubst %.hmm,%.hmm.h3m,$(HMMS))

#
# Rules
#

.PHONY: all clean test portclean pkg pkg.expand pkg.install

%.hmm.h3i: %.hmm
	echo $(HMMPRESS)
	(! test -s $< ) || $(HMMPRESS) -f $<

all:: pkg.install

pkg.expand::
	test -d $(PKGDIR) || mkdir $(PKGDIR)
	$(TAR) zxf $(PKGTAR) -C $(PKGDIR) --strip-components 1

pkg.install:: pkg
	@mkdir -p $(BINDIR)
	@cp $(PKGDIR)/ITSx $(BINDIR)
	@mkdir -p $(DATAITS)
	@cp -r $(PKGDIR)/ITSx_db $(DATAITS)
	@echo "+++++++++++ package $(PKG) done"

pkg:: pkg.expand clean $(HMMI) 

test::
	echo No test available
	
clean::
	\rm -f $(HMMP) $(HMMI) $(HMMM)
	echo Done

portclean::
	(! test -d $(PKGDIR)) || \rm -r $(PKGDIR)


include $(CFGDIR)targets/help.targ

