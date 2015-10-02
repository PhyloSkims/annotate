# ---------------------------------------------------------------
# @file: Makefile
# @desc: makefile prokov
#
# @history:
# @+	    <Gloup> : Sun Oct 31 1999 : Created
# @+   	    <Gloup> : Dec 00 : HelixWare port
# @+	    <Gloup> : Mar 01 : Corrected <if> bug
#
# @note: should be processed with gnu compatible make
# @note: helixware_compatible
#
# @end:
# ---------------------------------------------------------------
#

include ../../../../config/auto.conf

PROGS = prokov_learn prokov_curve \
	    prokov_score prokov_cds   \
	    prokov_orf

OSRC = fasta_io.c	    \
       libkov_util.c	\
	   libkov_io.c	    \
	   libkov_proba.c

include $(CFGDIR)targets/lxbin.targ
include $(CFGDIR)targets/help.targ

#
# file dependencies
#

prokov_learn : $(OBJ) prokov_learn.c
	$(CC) $(CFLAGS) -o $@ $@.c $(OBJ) $(LDFLAGS) $(LDLIBS)

prokov_curve : $(OBJ) prokov_curve.c
	$(CC) $(CFLAGS) -o $@ $@.c $(OBJ) $(LDFLAGS) $(LDLIBS)

prokov_score : $(OBJ) prokov_score.c
	$(CC) $(CFLAGS) -o $@ $@.c $(OBJ) $(LDFLAGS) $(LDLIBS)

prokov_cds : $(OBJ) prokov_cds.c
	$(CC) $(CFLAGS) -o $@ $@.c $(OBJ) $(LDFLAGS) $(LDLIBS)

prokov_orf : $(OBJ) prokov_orf.c
	$(CC) $(CFLAGS) -o $@ $@.c $(OBJ) $(LDFLAGS) $(LDLIBS)

