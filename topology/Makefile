CC=gcc
CFLAGS=-lm -I.
OBJ = topology_algorithm.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

topology_parbrain: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)


