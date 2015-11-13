#---------------------------------------------------------------
# $Id: $
# ---------------------------------------------------------------
# @file: Makefile
# @desc: makefile for muscle
# ---------------------------------------------------------------
#
include ../../../config/auto.conf

all: muscle install 

muscle:
	chmod +x ./mk
	(export CXX=$(CXX) && ./mk)

test:

install:
	cp muscle $(BINDIR)
	
clean:
	-\rm -f *.o muscle

portclean: clean
	-\rm -f $(BINDIR)/muscle
