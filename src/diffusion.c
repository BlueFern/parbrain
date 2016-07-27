#include "diffusion.h"


void diffusion(int block_number, double t, double *u, double *du, nvu_workspace *w)
{
	// TODO: Declare fluxes as a vector to calculate them in a loop.
	double flu_diff_K_0, flu_diff_K_1, flu_diff_K_2, flu_diff_K_3;

    // Offset into the neighbours array to get the indices
    // of the neighbours for the current tissue block.
    int neigh_offset = block_number * 4; // TODO: Declare the constant where appropriate.
    int W_neighbour = w->neighbours[neigh_offset + 0];
    int W_neighbour_offset = 0;
    int W_ghost_block_idx = 0;

    const double tau = 0.7; 	// (sec) the diffusion rate - characteristic time scale for ion to travel one cell length
    double state_K_e = u[K_e];	// variable to diffuse

    double W_ghost_block_val = 0;
    double W_neighbour_val = 0;

    // If the West neighbour is a ghost block.
    if(W_neighbour < 0)
    {
    	W_ghost_block_idx = - W_neighbour - 1;
    	W_ghost_block_val = w->ghost_blocks[W_ghost_block_idx].vars[DIFF_K]; // value set in init_ghost_blocks
    	flu_diff_K_0 = (W_ghost_block_val - state_K_e) / tau;
    }
    // Otherwise it is a proper tissue block.
    else
    {
    	W_neighbour_offset =  W_neighbour - block_number;
    	W_neighbour_val = u[W_neighbour_offset * w->neq + K_e];
    	flu_diff_K_0 = (W_neighbour_val - state_K_e) / tau;
    }

    int N_neighbour = w->neighbours[neigh_offset + 1];
    int N_neighbour_offset = 0;
    int N_ghost_block_idx = 0;
    double N_ghost_block_val = 0;
    double N_neighbour_val = 0;

    if(N_neighbour < 0)
    {
    	N_ghost_block_idx = - N_neighbour - 1;
    	N_ghost_block_val = w->ghost_blocks[N_ghost_block_idx].vars[DIFF_K];
    	flu_diff_K_1 = (N_ghost_block_val - state_K_e) / tau;
    }
    else
    {
    	N_neighbour_offset = N_neighbour - block_number;
    	N_neighbour_val = u[N_neighbour_offset * w->neq + K_e];
    	flu_diff_K_1 = (N_neighbour_val - state_K_e) / tau;
    }

    int E_neighbour = w->neighbours[neigh_offset + 2];
    int E_neighbour_offset = 0;
    int E_ghost_block_idx = 0;
    double E_ghost_block_val = 0;
    double E_neighbour_val = 0;

    if(E_neighbour < 0)
	{
		E_ghost_block_idx = - E_neighbour - 1;
		E_ghost_block_val = w->ghost_blocks[E_ghost_block_idx].vars[DIFF_K];
		flu_diff_K_2 = (E_ghost_block_val - state_K_e) / tau;
	}
	else
	{
		E_neighbour_offset = E_neighbour - block_number;
		E_neighbour_val = u[E_neighbour_offset * w->neq + K_e];
		flu_diff_K_2 = (E_neighbour_val - state_K_e) / tau;
	}

    int S_neighbour = w->neighbours[neigh_offset + 3];
    int S_neighbour_offset = 0;
    int S_ghost_block_idx = 0;
    double S_ghost_block_val = 0;
    double S_neighbour_val = 0;

	if(S_neighbour < 0)
	{
		S_ghost_block_idx = - S_neighbour - 1;
		S_ghost_block_val = w->ghost_blocks[S_ghost_block_idx].vars[DIFF_K];
		flu_diff_K_3 = (S_ghost_block_val - state_K_e) / tau;
	}
	else
	{
		S_neighbour_offset = S_neighbour - block_number;
		S_neighbour_val = u[S_neighbour_offset * w->neq + K_e];
		flu_diff_K_3 = (S_neighbour_val - state_K_e) / tau;
	}

/*
	printf("block number: %d, neigh_offset: %d, state_K_e: %f\n", block_number, neigh_offset, state_K_e);
	printf("W_neighbour: %d, W_ghost_block_idx: %d, W_ghost_block_val: %f, W_neighbour_offset: %d, W_neighbour_val: %f, flu_diff_K_0: %f\n", W_neighbour, W_ghost_block_idx, W_ghost_block_val, W_neighbour_offset, W_neighbour_val, flu_diff_K_0);
	printf("N_neighbour: %d, N_ghost_block_idx: %d, N_ghost_block_val: %f, N_neighbour_offset: %d, N_neighbour_val: %f, flu_diff_K_1: %f\n", N_neighbour, N_ghost_block_idx, N_ghost_block_val, N_neighbour_offset, N_neighbour_val, flu_diff_K_1);
	printf("E_neighbour: %d, E_ghost_block_idx: %d, E_ghost_block_val: %f, E_neighbour_offset: %d, E_neighbour_val: %f, flu_diff_K_2: %f\n", E_neighbour, E_ghost_block_idx, E_ghost_block_val, E_neighbour_offset, E_neighbour_val, flu_diff_K_2);
	printf("S_neighbour: %d, S_ghost_block_idx: %d, S_ghost_block_val: %f, S_neighbour_offset: %d, S_neighbour_val: %f, flu_diff_K_3: %f\n", S_neighbour, S_ghost_block_idx, S_ghost_block_val, S_neighbour_offset, S_neighbour_val, flu_diff_K_3);
*/
	//printf("Block %d | Fluxes W: %f, N: %f, E: %f, S: %f\n", block_number, flu_diff_K_0, flu_diff_K_1, flu_diff_K_2, flu_diff_K_3);

	// Add diffusion terms to K equation
	du[K_e] += flu_diff_K_0 + flu_diff_K_1 + flu_diff_K_2 + flu_diff_K_3 ;
}

void set_neighbours(int idx, int m, int n, int *neighbours)
{
	// Calculate row and column for the given block index.
	int i = idx % m;
	int j = idx / m;

	// West neighbour.
	if(j == 0)
	{
		// No West neighbour if we are in the first column.
		neighbours[0] = -1;
	}
	else
	{
		// The neighbour is the length of the column behind us.
		neighbours[0] = idx - m;
	}

	// North neighbour.
	if(i == (m - 1))
	{
		// No North neighbour if we are in the last row.
		neighbours[1] = -1;
	}
	else
	{
		// The neighbour is the next block.
		neighbours[1] = idx + 1;
	}

	// East neighbour.
	if(j == (n - 1))
	{
		// No East neighbour if we are in the last column.
		neighbours[2] = -1;
	}
	else
	{
		// The neighbour is the length of the column ahead of us.
		neighbours[2] = idx + m;
	}

	// South neighbour.
	if(i == 0)
	{
		// No South neighbour if we are in the first row.
		neighbours[3] = -1;
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
	//printf("Entering %s...\n", __FUNCTION__);

	// TODO: Free this memory.
	w->edge_indices = malloc(nlocal * mlocal * sizeof(int));

	for(int i = 0; i < mlocal; i++)
	{
		// West side edge.
		w->edge_indices[i] = i;
		// East side edge.
		w->edge_indices[mlocal + nlocal + i] = mlocal * (nlocal - 1) + i;

		// W_gb[i] = i;
		// E_gb[i] = mlocal * nlocal - mlocal + i;
	}

	for(int j = 0; j < nlocal; j++)
	{
		// North side edge.
		w->edge_indices[mlocal + j] = j * mlocal + (mlocal - 1);
		// South side edge.
		w->edge_indices[2 * mlocal + nlocal + j] = j * mlocal;

		// N_gb[j] = mlocal + j;
		// S_gb[j] = 2 * mlocal + nlocal + j;
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

	// All-side loop.
	int odd_even = 1;
	int idx_offset = 0;

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

			// TODO: Catch situations when there are no neighbours on the North side.
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



