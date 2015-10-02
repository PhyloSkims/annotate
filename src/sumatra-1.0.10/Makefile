EXEC=sumatra

SUMATRA_SRC= sumatra.c \
			 mtcompare_sumatra.c
			              
SUMATRA_OBJ= $(patsubst %.c,%.o,$(SUMATRA_SRC))


SRCS= $(SUMATRA_SRC)

LIB= -lfasta -llcs -lfile -lutils -lz -lm

include ./global.mk

all: $(EXEC)


########
#
# sumatra compilation
#
########

# executable compilation and link

sumatra: $(SUMATRA_OBJ) $(LIBFASTA) $(LIBLCS) $(LIBFILE) $(LIBUTILS)
	$(CC) $(LDFLAGS) -o $@ -pthread $(SUMATRA_OBJ) $(LIBFASTAPATH) $(LIBLCSPATH) $(LIBFILEPATH) $(LIBUTILSPATH) $(LIB)
	
########
#
# project management
#
########

clean:
	rm -f *.o	
	rm -f *.P
	rm -f $(EXEC)
	$(MAKE) -C ./sumalibs/libfasta clean
	$(MAKE) -C ./sumalibs/liblcs clean
	$(MAKE) -C ./sumalibs/libfile clean
	$(MAKE) -C ./sumalibs/libutils clean

		
