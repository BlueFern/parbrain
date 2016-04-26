UNAME_M = $(shell uname -m)
UNAME_S = $(shell uname -s)

CFALL = $(CFLAGS) $(TARGET_ARCH)  
# by default leave those empty. We don't need an entry on power7
# on mac /usr/local is in the standard path by default
INC=
LIB=
RPATH=
ifeq ($(UNAME_M), ppc64)
	CFARCHDEP = -m64 -mtune=power7 -mcpu=power7 -pthread -std=c99 -O2
	# using poe on ppc64
	MPCC = mpcc -compiler gcc 
        RPATH = -Wl,-rpath=$(LD_RUN_PATH)
else
	CFARCHDEP = -Wall -std=c99 -g  
	# regular wrapper everywhere else
	MPCC = mpicc
	# Need to adjust INC on ubuntu
	ifeq ($(UNAME_S), Linux)
		INC= -I/usr/include/suitesparse
	endif
endif

CF = $(CFALL) $(CFARCHDEP)
CS = -lcxsparse
OBJ = matops.o brain.o nvu.o adjacency.o solver.o diffusion.o
EXE = parBrainSim
SRC = src

all: tags $(OBJ) $(EXE)

tags: $(SRC)/*.h $(SRC)/*.c
	ctags $(SRC)/*.h $(SRC)/*.c

brain.o: $(SRC)/brain.c $(SRC)/brain.h $(SRC)/diffusion.h Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -c $(SRC)/brain.c 

adjacency.o: $(SRC)/adjacency.c $(SRC)/brain.h Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -c $(SRC)/adjacency.c 

nvu.o: $(SRC)/nvu.h $(SRC)/nvu.c Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -c $(SRC)/nvu.c

matops.o: $(SRC)/matops.c $(SRC)/matops.h Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -c $(SRC)/matops.c
	
solver.o: $(SRC)/solver.c $(SRC)/solver.h Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -c $(SRC)/solver.c	
	
diffusion.o: $(SRC)/diffusion.c $(SRC)/diffusion.h Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -c $(SRC)/diffusion.c		

#testmatops: testmatops.c matops.o Makefile
#	$(MPCC) $(CF) $(INC) $(LIB) -o testmatops testmatops.c matops.o $(CS) -lm 

#testmat2: testmat2.c matops.o Makefile
#	$(MPCC) $(CF) $(INC) $(LIB) -o testmat2 testmat2.c matops.o $(CS) -lm

#testbrain: testbrain.c matops.o brain.o adjacency.o nvu.o Makefile
#	$(MPCC) $(CF) $(INC) $(LIB) -o testbrain testbrain.c $(OBJ) $(CS) -lm

parBrainSim: $(SRC)/run_parbrain.c matops.o brain.o nvu.o adjacency.o solver.o diffusion.o Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -o parBrainSim $(SRC)/run_parbrain.c $(OBJ) $(RPATH) $(CS) -lm

clean:
	rm -r tags $(OBJ) $(EXE)
