#ifndef DIFFUSION_H_
#define DIFFUSION_H_

#include "nvu.h"
#include "brain.h"

// Forward declaration to avoid compilation errors.
typedef struct workspace workspace;

// Number of neighbours for each ghost block.
static const int NUM_NEIGHBOURS = 4;

// Enumerator to keep track of the diffusion variables positions.
enum diff_idx
{
	DIFF_K , DIFF_NA, DIFF_Kk, DIFF_Vk, DIFF_NAk, DIFF_HCO3k, DIFF_CAk, DIFF_CLk
};

// Ghost block to store diffusion variables. Ghost blocks are placed around
// the perimeter of the tissue blocks allocated to an MPI process.
typedef struct ghost_block
{
	double *vars;
} ghost_block;


// Get the indices for the neighbours of a block in a grid of blocks.
void set_neighbours(int idx, int m, int n, int *neighbours);

// Get the indices for the four neighbouring MPI domains of a given domain.
void set_domain_neighbours(int idx, int m, int n, int *neighbours);

// Get the indices for the four immediate neighbours of a given block within an MPI domain.
void set_block_neighbours(int nlocal, int mlocal, nvu_workspace* w);

// Set the indices of the edges
void set_edge_indices(int nlocal, int mlocal, nvu_workspace *w);

// Allocate the space for ghost blocks.
void init_ghost_blocks(int nlocal, int mlocal, nvu_workspace *w);

// Calculate diffusion for tissue blocks within given MPI domain.
void diffusion(int block_number, double t, double *u, double *du, nvu_workspace *w);

void update_ghost_blocks(workspace *W, double *y);

#endif
