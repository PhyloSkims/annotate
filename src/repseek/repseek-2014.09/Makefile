#
# Design for Repseek
# compile all src and output one binary
# <amikezor>
# April 2004
# 

# MACHINE = SUN
# MACHINE = SGI
# MACHINE = LINUX
# MACHINE = OSF

include ../../../config/auto.conf



# SGI: PROTO= PROTO=1
# others: PROTO= PROTO=0
PROTO= PROTO=0

##### defined

SHELL = bash

SRC  = sort.c\
       help.c\
       output.c\
       filter2seq.c \
       filter.c \
       families.c \
       families_2seqs.c\
       memory.c\
       memory_align.c\
       sequence.c\
       readfst.c\
       lmin.c \
       smin.c\
       KMRK.c \
       KMRK_Seeds.c \
       KMRK_merge_seeds.c\
       KMRK_mask.c\
       align_matrix.c\
       align_seeds.c\
       align_di.c\
       align_blast2like.c\
       read_seeds.c


OBJ  = $(SRC:.c=.o)


## Rules

all: repseek install

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<;


repseek:  $(OBJ) main_repseek.c
	$(CC) $(CFLAGS) -o $@ $(OBJ) main_repseek.c $(LDFLAGS);


install: repseek
	cp  repseek $(BINDIR)


clean:
	+rm -f *.o repseek

archive:
	\rm -rf RepSeek; mkdir RepSeek
	cp -r $(SRC) *.h Makefile main_repseek.c *.txt RepSeek;
	tmp=`date | awk '{print $$3$$2$$NF}'`;tar czvf RepSeek.$$tmp.tgz  RepSeek;
	\rm -rf RepSeek

