
CC=gcc
LDFLAGS=

CFLAGS = -O3 -w

default: all

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $< $(LIB)
