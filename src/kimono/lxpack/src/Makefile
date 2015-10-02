# ---------------------------------------------------------------
# $Id: $
#
# @file: Makefile
# @desc: makefile for kimono/src
#
# @history:
#
# @note: should be processed with gnu compatible make
# @note: helixware_compatible
#
# @end:
# ---------------------------------------------------------------
#

include ../../../../config/auto.conf

PROGS  = kimono kimfit

OSRC   = fasta_io.c	 \
	 kim_util.c	 \
	 kim_alphabet.c	 \
	 kim_sequence.c	 \
	 kim_zone.c	 \
	 kim_indexer.c	 \
	 kim_counter.c	 \
	 kim_mixer.c	 \
	 kim_distance.c	 \
	 kim_genetic.c	 \
	 kim_codonskew.c

include ../../../../config/targets/lxbin.targ
include ../../../../config/targets/help.targ

#
# file dependencies
#

kimono:	$(OBJ) kimono_main.c kimono_help.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

kimfit:	$(OBJ) kimfit_main.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)
