# ---------------------------------------------------------------
# @pckg: SeqUtils / V1.0 / Gloup Jan 92
#
# @file: Makefile
# @desc: makefile for Utilsrc /  Biological utilities
#
# @history:
# @+	    <Gloup> : Jan 96 : Created
# @+	    <Gloup> : Jul 93 : revised version
# @+	    <Gloup> : Feb 01 : Adapted to helixware
#
# @note: should be processed with gnu compatible make
# @note: helixware_compatible
#
# @end:
# ---------------------------------------------------------------
#

include ../../config/auto.conf

#
# Machine independant flags
#

PROGS = util_complinv	\
	util_cut	        \
	util_translate	    \
	util_tab_aa	        \
	util_tab_codon	    \
	util_tab_dinuc	    \
	util_tab_nuc	    \
	util_skew	        \
	util_skew_teta	    \
	util_skew_chi_mono	\
	util_skew_mono		\
	util_codon_skew		\
	util_codon_skew_corr

OSRC     = $(PROGS:=.c)

include ../../config/targets/lxbin.targ
include ../../config/targets/help.targ

INCDIR	  =	../../include

LIBS	  =	-lfasta -lbio -laabi

