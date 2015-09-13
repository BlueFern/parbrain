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
OBJ = matops.o brain.o nvu.o adjacency.o 
#EXE = testmatops testmat2 simulate testbrain
EXE = simulate

all: tags $(OBJ) $(EXE)

tags: *.h *.c
	ctags *.h *.c

brain.o: brain.c brain.h Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -c brain.c 

adjacency.o: adjacency.c brain.h Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -c adjacency.c 

nvu.o: nvu.h nvu.c Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -c nvu.c

matops.o: matops.h matops.c Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -c matops.c

testmatops: testmatops.c matops.o Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -o testmatops testmatops.c matops.o $(CS) -lm 

testmat2: testmat2.c matops.o Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -o testmat2 testmat2.c matops.o $(CS) -lm

testbrain: testbrain.c matops.o brain.o adjacency.o nvu.o Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -o testbrain testbrain.c $(OBJ) $(CS) -lm

simulate: simulate.c matops.o brain.o nvu.o adjacency.o Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -o simulate simulate.c $(OBJ) $(RPATH) $(CS) -lm

clean:
	rm -r tags $(OBJ) $(EXE) *.dSYM
