EXEC=sumaclust

SUMACLUST_SRC= sumaclust.c \
               mtcompare_sumaclust.c

             
SUMACLUST_OBJ= $(patsubst %.c,%.o,$(SUMACLUST_SRC))


SRCS= $(SUMACLUST_SRC)

LIB= -lfasta -llcs -lfile -lutils -lm #-ll


include ./global.mk

all: $(EXEC) install


########
#
# sumaclust compilation
#
########

# executable compilation and link

#ifeq ($(CC),gcc)
#        LFLAGS = -fopenmp 
#else
#        LFLAGS = 
#endif

sumaclust: $(SUMACLUST_OBJ) $(LIBFASTA) $(LIBLCS) $(LIBFILE) $(LIBUTILS)
	$(CC) $(LIBFASTAPATH) $(LIBLCSPATH) $(LIBFILEPATH) $(LIBUTILSPATH) $(LDFLAGS) -o $@ $(LFLAGS) $(SUMACLUST_OBJ) $(LIB) 
	
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

portclean: clean
	-rm -f $(BINDIR)/$(EXEC)

install:
	cp $(EXEC) $(BINDIR)

		
