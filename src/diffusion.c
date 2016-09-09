#include "diffusion.h"

void diffusion(int block_number, double t, double *u, double *du, nvu_workspace *w)
{
    const double tau = 4.3; // (sec) The diffusion rate - characteristic time scale for ion to travel one cell length.

    // Offset into the neighbours array to get the indices
    // of the neighbours for the current tissue block.
    int neigh_offset = block_number * NUM_NEIGHBOURS;

    double state_K_e = u[K_e];	// Extra-cellular potassium in our block.

    // Iterate over all neighbours.
    for(int neigh_idx = 0; neigh_idx < NUM_NEIGHBOURS; neigh_idx++)
    {
    	// Get the index of the neighbour.
    	int neigh_idx_val = w->neighbours[neigh_offset + neigh_idx];

    	// Neighbour's potassium value to be retrieved.
    	double neighbour_K_e = 0;

    	if(neigh_idx_val < 0)
    	{
    		// If the neighbour is a ghost block, get the value from it.
    		neighbour_K_e = w->ghost_blocks[-neigh_idx_val - 1].vars[DIFF_K];
    	}
    	else
    	{
    		// Otherwise it is a real tissue block.
    		neighbour_K_e = u[(neigh_idx_val - block_number) * w->neq + K_e];
    	}

    	// Calculate the flux.
    	double flu_diff_K = (neighbour_K_e - state_K_e) / tau;

    	// Update the derivative.
    	du[K_e] += flu_diff_K;
    }
}

/*
 * Populate the neighbours array with the indices of neighbours for the
 * give block number.
 *
 * Indexing in a grid of blocks starts at the lower left corner.
 * Block indices inside a single domain are organised in column-major order.
 *
 * Ghost blocks around the perimeter are numbered with a 1-base negative
 * index in the following order:
 *
 * West side, bottom-up.
 * North side, left-to-right.
 * East side, bottom-up.
 * South side, left-to-right.
 *
 * Yes, this is a bit wacky.
 *
 *                boing         boing         boing
 *  e-e           . - .         . - .         . - .
 * (\_/)\       '       `.   ,'       `.   ,'       .
 *  `-'\ `--.___,         . .           . .          .
 *    '\( ,_.-'
 *        \\               "             "            a:f
 *        ^'
 */

void set_neighbours(int idx, int m, int n, int *neighbours)
{
	// Calculate row and column for the given block index.
	int i = idx % m;
	int j = idx / m;

	// West neighbour.
	if(j == 0)
	{
		// No actual West neighbour if we are in the first column.
		// Ghost block number is the row number, converted to 1-based negative index.
		neighbours[0] = -(i + 1);
	}
	else
	{
		// The neighbour is the length of the column behind us.
		neighbours[0] = idx - m;
	}

	// North neighbour.
	if(i == (m - 1))
	{
		// No actual North neighbour if we are in the last row.
		// Ghost block number is West side length (m) plus column number, converted to 1-based negative index.
		neighbours[1] = -(m + j + 1);
	}
	else
	{
		// The neighbour is the next block.
		neighbours[1] = idx + 1;
	}

	// East neighbour.
	if(j == (n - 1))
	{
		// No actual East neighbour if we are in the last column.
		// Ghost block number is West and North sides lengths (m + n) plus row number, converted to 1-based negative index.
		neighbours[2] = -(m + n + i + 1);
	}
	else
	{
		// The neighbour is the length of the column ahead of us.
		neighbours[2] = idx + m;
	}

	// South neighbour.
	if(i == 0)
	{
		// No actual South neighbour if we are in the first row.
		// Ghost block number is West, North and East sides (m + n + m) plus column number, converted to 1-based negative index.
		neighbours[3] = -(m + n + m + j + 1);
	}
	else
	{
		// The neighbour is the previous block.
		neighbours[3] = idx -1;
	}

#if 0
	int u = 0;
	printf("index %d: ", idx);
	for (u=0; u < 4; u++)
	{
		printf("%d ", neighbours[u]);
	}
	printf("\n");
#endif

}

// Get the indices for all neighbours for all tissue blocks in the given MPI domain.
void set_block_neighbours(int nlocal, int mlocal, nvu_workspace* w)
{
	// Each block has four neighbours.
	// TODO: This also needs to be freed somewhere along with ghost blocks.

	w->neighbours = malloc(nlocal * mlocal * 4 * sizeof(int));

	int block_offset = 4; // neighbours per block

	// neighbours array gets filled in the following way:
	// W0, N0, E0, S0, W1, N1, ... etc.

	for (int idx = 0; idx < nlocal * mlocal; idx++)
	{
		set_neighbours(idx, mlocal, nlocal, w->neighbours + block_offset * idx);
	}
}

void set_domain_neighbours(int idx, int m, int n, int *neighbours)
{
	set_neighbours(idx, m, n, neighbours);
}

void set_edge_indices(int nlocal, int mlocal, nvu_workspace *w)
{
	// TODO: Free this memory.
	w->edge_indices = malloc(nlocal * mlocal * sizeof(int));

	for(int i = 0; i < mlocal; i++)
	{
		// West side edge.
		w->edge_indices[i] = i;
		// East side edge.
		w->edge_indices[mlocal + nlocal + i] = mlocal * (nlocal - 1) + i;
	}

	for(int j = 0; j < nlocal; j++)
	{
		// North side edge.
		w->edge_indices[mlocal + j] = j * mlocal + (mlocal - 1);
		// South side edge.
		w->edge_indices[2 * mlocal + nlocal + j] = j * mlocal;
	}
}

void init_ghost_blocks(int nlocal, int mlocal, nvu_workspace *w)
{
	// Ghost block surround the perimeter of the MPI domain.
	w->num_ghost_blocks = 2 * (nlocal + mlocal);

	// TODO: Free the space in a function deallocating the space for nvu_workspace, which is to be written and called at the right place.
    w->ghost_blocks = malloc(w->num_ghost_blocks * sizeof(ghost_block));

    for(int i = 0; i < w->num_ghost_blocks; i++)
    {
    	w->ghost_blocks[i].vars = malloc(NUM_DIFF_VARS * sizeof(double));

    	w->ghost_blocks[i].vars[DIFF_K] = 3e3;
    }
}

void update_ghost_blocks(workspace *W, double *y)
{

	int odd_even = 1;
	int idx_offset = 0;

	// Exchange ghost block values separately for each side.
	for(int side_idx = 0; side_idx < 4; side_idx++)
	{
		// Choose length of the side depending on the odd or even value.
		int side_length = odd_even ? W->mlocal : W->nlocal;

		// If neighbour exists.
		if(W->domain_neighbours[side_idx] >= 0)
		{
			// TODO: This likely will have to be a derived type when we have more than one diffusion variable to post.
			double *send_array = malloc(sizeof(double) * side_length);
			double *receive_array = malloc(sizeof(double) * side_length);

			// Populate send array with values from the state variables array.
			for(int i = 0; i < side_length; i++)
			{
				// Get the state variable based on the indices of the edge block for this side (and all other sides).
				send_array[i] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + K_e];
			}

			// Simultaneously send and receive the values for the ghost blocks.
			MPI_Sendrecv(
					send_array, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					receive_array, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			// Copy received data into the ghost blocks.
			for(int i = 0; i < side_length; i++)
			{
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_K] = receive_array[i];
			}

			free(send_array);
			free(receive_array);
		}
		// If there is no neighbour, copy the current block's value into the ghost block.
		else
		{
			for(int i = 0; i < side_length; i++)
			{
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_K] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + K_e];
			}
		}
		idx_offset += side_length;
		odd_even = !odd_even;
	}
}
