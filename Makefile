CC = gcc
CFLAGS = -Wall -O2 -lm
DEPS = routines.h system.h
OBJ = main.o routines.o system.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

solver: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f *.o *.csv solver