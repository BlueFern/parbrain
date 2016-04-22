#ifndef DIFFUSION_H_
#define DIFFUSION_H_

#include "nvu.h"


// Number of variables stored in diffusion structs.
static const int NUM_DIFF_VARS = 1;

// Enumerator to keep track of the diffusion variables positions.
enum diff_idx
{
	DIFF_K
};

// Ghost block to store diffusion variables. Ghost blocks are placed around
// the 'perimeter' of the tissue blocks allocated to an MPI process.
typedef struct ghost_block
{
	double *vars;
} ghost_block;


// Get the indices for the four immediate neighbours of a given block.
void set_neighbours(int idx, int m, int n, int *neighbours);

// Same as set_neighbours?
void set_domain_neighbours(int idx, int m, int n, int *neighbours);

// Get the indices for all neighbours for all tissue blocks in the given MPI domain.
void set_block_neighbours(int nlocal, int mlocal, nvu_workspace* w);

// Set the indices of the edges
void set_edge_indices(int nlocal, int mlocal, nvu_workspace *w);

// Allocate the space for ghost blocks.
void init_ghost_blocks(int nlocal, int mlocal, nvu_workspace *w);

// Calculate diffusion for tissue blocks within given MPI domain.
void diffusion(int block_number, double t, double *u, double *du, nvu_workspace *w);



#endif
