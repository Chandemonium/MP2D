CC=g++
CFLAGS=-g
CFLAGS=-O3
DEPS = MP2D.h
OBJ = Main.o MP2D.o


%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

MP2D: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) 

.PHONY: clean

clean: 
	rm -f *.o MP2D
