# ---------------------------------------------------------------
# $Id: $
# ---------------------------------------------------------------
# @file: Makefile
# @desc: makefile for package ebi exonerate
#
# @history:
# @+	    <Gloup> : Apr 97 : Created
# @+	    <Gloup> : Mar 02 : Updated for LXxware
# @+	    <EC> : May 2023 : Adapted to ORG.Annot
#
# @note: should be processed with gnu compatible make
# @note: helixware_compatible
#
# @end:
# ---------------------------------------------------------------
#
#
# Due to the long compilation time of ncbi blast sources (~1h),
# the default build actually consists in extracting pre-compiled
# binaries.
#
# To revert to full compilation from sources, then change
# the following line

DIRS = binaries

# to
# 
# DIRS = src

include ../../config/targets/propagate.targ
include ../../config/targets/help.targ

all::
	@echo "+++++++++++ ebi exonerate package ($(DIRS)) done"
