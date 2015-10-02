include ../../../../../config/auto.conf


default: all

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $< 
