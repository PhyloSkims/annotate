include ../../../../config/auto.conf

all: muscle install 

muscle:
	chmod +x ./mk
	(export CXX=$(CXX) && ./mk)

install:
	cp muscle $(BINDIR)
	
