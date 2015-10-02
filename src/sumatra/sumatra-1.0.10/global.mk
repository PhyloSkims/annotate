include ../../../config/auto.conf

LIBFASTAPATH = -L./sumalibs/libfasta
LIBLCSPATH   = -L./sumalibs/liblcs
LIBFILEPATH  = -L./sumalibs/libfile
LIBUTILSPATH = -L./sumalibs/libutils

LIBFASTA = ./sumalibs/libfasta/libfasta.a
LIBLCS   = ./sumalibs/liblcs/liblcs.a
LIBFILE  = ./sumalibs/libfile/libfile.a
LIBUTILS = ./sumalibs/libutils/libutils.a



#ifeq ($(CC),gcc)
#        CFLAGS = -O3 -s -DOMP_SUPPORT -fopenmp -w
#else
#        CFLAGS = -O3 -w
#endif


default: all

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $< 


########
#
# libraries compilation
#
########

./sumalibs/libfasta/libfasta.a:
	$(MAKE) -C ./sumalibs/libfasta

./sumalibs/liblcs/liblcs.a:
	$(MAKE) -C ./sumalibs/liblcs

./sumalibs/libfile/libfile.a:
	$(MAKE) -C ./sumalibs/libfile

./sumalibs/libutils/libutils.a:
	$(MAKE) -C ./sumalibs/libutils
	
