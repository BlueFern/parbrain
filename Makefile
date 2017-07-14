# Optimised for use on brats01 machine

## Do this in terminal if it can't find libcxsparse
## LD_LIBRARY_PATH=/opt/pgi/linux86-64/2017/SuiteSparse/4.5.5
## export LD_LIBRARY_PATH

## Do this every log in
# module load pgi-openmpi and suitesparse

# Default values if they are not set, optimised for use on brats01 (change for local computer)
SUITESPARSE ?= /opt/pgi/linux86-64/2017/SuiteSparse/4.5.5
MPICC ?= mpicc
CFLAGS ?= -O2
RPATH ?= -Wl,-rpath=$(SUITESPARSE)/lib

LIB = -lcxsparse -lm
MYCPPFLAGS = -I$(SUITESPARSE)/include
MYLDFLAGS = -L$(SUITESPARSE)/lib

# If LIBRARY_PATH is set, use it to overwrite RPATH and ignore MYCPPFLAGS and MYLDFLAGS
ifneq ($(LIBRARY_PATH),)
	RPATH = -Wl,-rpath=$(LIBRARY_PATH)
	MYCPPFLAGS =
	MYLDFLAGS =
endif

CPPFLAGS += $(MYCPPFLAGS)
LDFLAGS += $(MYLDFLAGS) $(LIB)

OBJ = run_parbrain.o matops.o brain.o nvu.o adjacency.o solver.o diffusion.o
EXE = parBrainSim
SRC = src

all: tags $(OBJ) $(EXE)

tags: $(SRC)/*.h $(SRC)/*.c
	ctags $(SRC)/*.h $(SRC)/*.c

brain.o: $(SRC)/brain.c $(SRC)/brain.h $(SRC)/diffusion.h
	$(MPICC) $(CFLAGS) $(CPPFLAGS) -c $(SRC)/brain.c 

adjacency.o: $(SRC)/adjacency.c $(SRC)/brain.h
	$(MPICC) $(CFLAGS) $(CPPFLAGS) -c $(SRC)/adjacency.c 

nvu.o: $(SRC)/nvu.h $(SRC)/nvu.c
	$(MPICC) $(CFLAGS) $(CPPFLAGS) -c $(SRC)/nvu.c

matops.o: $(SRC)/matops.c $(SRC)/matops.h
	$(MPICC) $(CFLAGS) $(CPPFLAGS) -c $(SRC)/matops.c
	
solver.o: $(SRC)/solver.c $(SRC)/solver.h
	$(MPICC) $(CFLAGS) $(CPPFLAGS) -c $(SRC)/solver.c	
	
diffusion.o: $(SRC)/diffusion.c $(SRC)/diffusion.h
	$(MPICC) $(CFLAGS) $(CPPFLAGS) -c $(SRC)/diffusion.c

run_parbrain.o: $(SRC)/run_parbrain.c $(SRC)/solver.h
	$(MPICC) $(CFLAGS) $(CPPFLAGS) -c $(SRC)/run_parbrain.c	

parBrainSim: run_parbrain.o matops.o brain.o nvu.o adjacency.o solver.o diffusion.o
	$(MPICC) $(CFLAGS) $(CPPFLAGS) $(RPATH) -o parBrainSim $(OBJ) $(LDFLAGS) 

clean:
	rm -rf tags $(OBJ) $(EXE)
