# Default values
MPICC ?= mpicc

CFLAGS ?= -O2
CS = -lcxsparse
MYLDFLAGS = $(LDFLAGS) $(LIB) $(CS) -lm
OBJ = matops.o brain.o nvu.o adjacency.o solver.o diffusion.o
EXE = parBrainSim
SRC = src

all: tags $(OBJ) $(EXE)

tags: $(SRC)/*.h $(SRC)/*.c
	ctags $(SRC)/*.h $(SRC)/*.c

brain.o: $(SRC)/brain.c $(SRC)/brain.h $(SRC)/diffusion.h Makefile
	$(MPICC) $(CFLAGS) $(CPPFLAGS) -c $(SRC)/brain.c 

adjacency.o: $(SRC)/adjacency.c $(SRC)/brain.h Makefile
	$(MPICC) $(CFLAGS) $(CPPFLAGS) -c $(SRC)/adjacency.c 

nvu.o: $(SRC)/nvu.h $(SRC)/nvu.c Makefile
	$(MPICC) $(CFLAGS) $(CPPFLAGS) -c $(SRC)/nvu.c

matops.o: $(SRC)/matops.c $(SRC)/matops.h Makefile
	$(MPICC) $(CFLAGS) $(CPPFLAGS) -c $(SRC)/matops.c
	
solver.o: $(SRC)/solver.c $(SRC)/solver.h Makefile
	$(MPICC) $(CFLAGS) $(CPPFLAGS) -c $(SRC)/solver.c	
	
diffusion.o: $(SRC)/diffusion.c $(SRC)/diffusion.h Makefile
	$(MPICC) $(CFLAGS) $(CPPFLAGS) -c $(SRC)/diffusion.c		

parBrainSim: $(SRC)/run_parbrain.c matops.o brain.o nvu.o adjacency.o solver.o diffusion.o Makefile
	$(MPICC) $(CFLAGS) $(CPPFLAGS) -o parBrainSim $(SRC)/run_parbrain.c $(OBJ) $(MYLDFLAGS)

clean:
	rm -rf tags $(OBJ) $(EXE)
