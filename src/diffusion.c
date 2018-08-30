#include "diffusion.h"

void diffusion(int block_number, double t, double *u, double *du, nvu_workspace *w)
{

    // Offset into the neighbours array to get the indices
    // of the neighbours for the current tissue block.
    int neigh_offset = block_number * NUM_NEIGHBOURS;

    double state_K_e = u[i_K_e];	// Extracellular potassium in our block.
    double state_Na_e = u[i_Na_e];	// Extracellular sodium in our block.
    double state_K_k = u[i_K_k];	// Astrocytic potassium in our block.
    double state_v_k = u[i_v_k];	// Astrocytic membrane potential in our block.
    double state_Na_k = u[i_Na_k];	// Astrocytic sodium in our block.
    double state_HCO3_k = u[i_HCO3_k];	// Astrocytic HCO3 in our block.
    double state_Cl_k = u[i_Cl_k];	// Astrocytic chlorine in our block.
    double state_Ca_k = u[i_ca_k];	// Astrocytic calcium in our block.
    
    // Curvature ***********
    double diffusion_scaling;
    if (t> theta_t_on && THETA_SWITCH)  // To make the whole area the same theta value after theta_t_on (for wave segment or constant curvature simulations)
    {
		double theta = theta_all_space;
		diffusion_scaling = (M_PI/a_th) * ( cosh(eta_th) - cos(theta) );
//		printf("C = %1.2f\n", diffusion_scaling);
	}
	else
	{
	    // *** REMOVE!
	    if (u[i_coup] < 0.9 && TEMP_SWITCH)
        {
            diffusion_scaling = 0.1;
//            diffusion_scaling = u[i_coup] - 0.3; // if C < 0.8?
        }
        else
        {
            diffusion_scaling = u[i_coup];
        }
	}

    // Iterate over all neighbours.
    for(int neigh_idx = 0; neigh_idx < NUM_NEIGHBOURS; neigh_idx++)
    {
    	// Get the index of the neighbour.
    	int neigh_idx_val = w->neighbours[neigh_offset + neigh_idx];

    	// Neighbour's values to be retrieved.
    	double neighbour_K_e = 0;
    	double neighbour_Na_e = 0;
    	double neighbour_K_k = 0;
    	double neighbour_v_k = 0;
    	double neighbour_Na_k = 0;
    	double neighbour_HCO3_k = 0;
    	double neighbour_Cl_k = 0;
    	double neighbour_Ca_k = 0;

    	if(neigh_idx_val < 0)
    	{
    		// If the neighbour is a ghost block, get the value from it.
    		neighbour_K_e = w->ghost_blocks[-neigh_idx_val - 1].vars[DIFF_K];
    		neighbour_Na_e = w->ghost_blocks[-neigh_idx_val - 1].vars[DIFF_NA];
    		neighbour_K_k = w->ghost_blocks[-neigh_idx_val - 1].vars[DIFF_Kk];
    		neighbour_v_k = w->ghost_blocks[-neigh_idx_val - 1].vars[DIFF_Vk];
    		neighbour_Na_k = w->ghost_blocks[-neigh_idx_val - 1].vars[DIFF_NAk];
    		neighbour_HCO3_k = w->ghost_blocks[-neigh_idx_val - 1].vars[DIFF_HCO3k];
    		neighbour_Cl_k = w->ghost_blocks[-neigh_idx_val - 1].vars[DIFF_CLk];
    		neighbour_Ca_k = w->ghost_blocks[-neigh_idx_val - 1].vars[DIFF_CAk];
    	}
    	else
    	{
    		// Otherwise it is a real tissue block.
    		neighbour_K_e = u[(neigh_idx_val - block_number) * w->neq + i_K_e];
    		neighbour_Na_e = u[(neigh_idx_val - block_number) * w->neq + i_Na_e];
    		neighbour_K_k = u[(neigh_idx_val - block_number) * w->neq + i_K_k];
    		neighbour_v_k = u[(neigh_idx_val - block_number) * w->neq + i_v_k];
    		neighbour_Na_k = u[(neigh_idx_val - block_number) * w->neq + i_Na_k];
    		neighbour_HCO3_k = u[(neigh_idx_val - block_number) * w->neq + i_HCO3_k];
    		neighbour_Cl_k = u[(neigh_idx_val - block_number) * w->neq + i_Cl_k];
    		neighbour_Ca_k = u[(neigh_idx_val - block_number) * w->neq + i_ca_k];
    	}

    	// Calculate the flux.
    	double flu_diff_K 	= (D_Ke / pow(delta_x,2)) * (neighbour_K_e - state_K_e);
    	double flu_diff_Na 	= (D_Nae / pow(delta_x,2)) * (neighbour_Na_e - state_Na_e);
    	
//    	double flu_GJ_K 	= (D_gap / pow(delta_x,2)) * (neighbour_K_k - state_K_k);
//    	double flu_GJ_K 	= (D_gap / pow(delta_x,2)) * ( ((z_K * Farad)/(R_gas * Temp)) 	* ( (neighbour_K_k + state_K_k)/2 ) 		* (neighbour_v_k - state_v_k) );
    	    	
    	double flu_GJ_K 	= (D_gap / pow(delta_x,2)) * (  (neighbour_K_k - state_K_k) 		+ ((z_K * Farad)/(R_gas * Temp)) 	* ( (neighbour_K_k + state_K_k)/2 ) 		* (neighbour_v_k - state_v_k) );
    	double flu_GJ_Na 	= (D_gap / pow(delta_x,2)) * (  (neighbour_Na_k - state_Na_k) 		+ ((z_Na * Farad)/(R_gas * Temp)) 	* ( (neighbour_Na_k + state_Na_k)/2 ) 		* (neighbour_v_k - state_v_k) );
    	double flu_GJ_HCO3 	= (D_gap / pow(delta_x,2)) * (  (neighbour_HCO3_k - state_HCO3_k) 	+ ((z_HCO3 * Farad)/(R_gas * Temp)) * ( (neighbour_HCO3_k + state_HCO3_k)/2 ) 	* (neighbour_v_k - state_v_k) );
    	double flu_GJ_Cl	= (D_gap / pow(delta_x,2)) * (  (neighbour_Cl_k - state_Cl_k) 		+ ((z_Cl * Farad)/(R_gas * Temp)) 	* ( (neighbour_Cl_k + state_Cl_k)/2 ) 		* (neighbour_v_k - state_v_k) );
    	double flu_GJ_Ca 	= (D_gap / pow(delta_x,2)) * (  (neighbour_Ca_k - state_Ca_k) 		+ ((z_Ca * Farad)/(R_gas * Temp)) 	* ( (neighbour_Ca_k + state_Ca_k)/2 ) 		* (neighbour_v_k - state_v_k) );
    	    	    	    			
    	// Update the derivatives.
    	if (DIFFUSION_SWITCH == 1)
		{
        	// Calculate the flux.
        	double flu_diff_K 	= (D_Ke / pow(delta_x,2)) * (neighbour_K_e - state_K_e);
        	double flu_diff_Na 	= (D_Nae / pow(delta_x,2)) * (neighbour_Na_e - state_Na_e);
    		
			du[i_K_e] += diffusion_scaling * flu_diff_K;
			du[i_Na_e] += diffusion_scaling * flu_diff_Na;
	//    	// As these two variables contain the term "SC_coup * du[K_e] * 1e3", they must also be updated!
			du[i_K_s] += SC_coup * diffusion_scaling * flu_diff_K * 1e3;
			du[i_Na_s] += -SC_coup * diffusion_scaling * flu_diff_K * 1e3;
		}
    	else if (DIFFUSION_SWITCH == 2)
    	{
    		double Delta_K_e = neighbour_K_e - state_K_e;
    		double Delta_Na_e = neighbour_Na_e - state_Na_e;
    		double average_K_e = (neighbour_K_e + state_K_e)/2;
    		double average_Na_e = (neighbour_Na_e + state_Na_e)/2;
    		double Delta_v_e =  -(R_gas * Temp) / Farad * ( z_K * D_Ke * Delta_K_e + z_Na * D_Nae * Delta_Na_e ) / ( pow(z_K,2) * D_Ke * average_K_e + pow(z_Na,2) * D_Nae * average_Na_e );
    		  	
//    		if (fmod(t,1)<1e-5)
//    		{
//    		printf("t %f, Dv %f, v_k %f\n", t, Delta_v_e, state_v_k);
//    		}
    		
        	// Calculate the flux.
        	double flu_diff_K 	= (D_Ke / pow(delta_x,2)) * (  Delta_K_e + ((z_K * Farad)/(R_gas * Temp)) * ( average_K_e * Delta_v_e ) );
        	double flu_diff_Na 	= (D_Nae / pow(delta_x,2)) * (  Delta_Na_e + ((z_Na * Farad)/(R_gas * Temp)) * ( average_Na_e * Delta_v_e ) );
    		
			du[i_K_e] += diffusion_scaling * flu_diff_K;
			du[i_Na_e] += diffusion_scaling * flu_diff_Na;
	//    	// As these two variables contain the term "SC_coup * du[K_e] * 1e3", they must also be updated!
			du[i_K_s] += SC_coup * diffusion_scaling * flu_diff_K * 1e3;
			du[i_Na_s] += -SC_coup * diffusion_scaling * flu_diff_K * 1e3;
    	}
    	if (GJ_SWITCH == 1) // Just K gap junctions
    	{
    		du[i_K_k] 		+= flu_GJ_K;
    		du[i_Cl_k] 		+= z_K * flu_GJ_K;
//    		printf("t %f, before dv %f\n", t, du[i_v_k]);
    		du[i_v_k] 		+= gam * (z_K * flu_GJ_K);
    	}
    	else if (GJ_SWITCH == 2) // All ion gap junctions
    	{	
    		du[i_K_k] 		+= flu_GJ_K;
//    		du[i_Na_k] 		+= flu_GJ_Na;
//    		du[i_HCO3_k]	+= flu_GJ_HCO3;
    		du[i_ca_k]		+= flu_GJ_Ca;
    		du[i_Cl_k] 		+= z_K * flu_GJ_K + z_Ca * flu_GJ_Ca;//
    		du[i_v_k] 		+= gam * ( z_K * flu_GJ_K + z_Ca * flu_GJ_Ca);
//    		du[i_Cl_k] 		+= z_Na * flu_GJ_Na + z_HCO3 * flu_GJ_HCO3 + z_Ca * flu_GJ_Ca;//
//    		du[i_v_k] 		+= gam * ( z_Na * flu_GJ_Na + z_HCO3 * flu_GJ_HCO3 + z_Cl * flu_GJ_Cl + z_Ca * flu_GJ_Ca);
    		
    	}
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
 *        \\               "             "
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
	w->edge_indices = malloc((nlocal + mlocal) * 2 * sizeof(int));

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

    	w->ghost_blocks[i].vars[DIFF_K] = 3.493;
    	w->ghost_blocks[i].vars[DIFF_NA] = 150;
    	w->ghost_blocks[i].vars[DIFF_Kk] = 92708;
    	w->ghost_blocks[i].vars[DIFF_Vk] = -88.9;
    	w->ghost_blocks[i].vars[DIFF_NAk] = 18268;
    	w->ghost_blocks[i].vars[DIFF_HCO3k] = 9131;
    	w->ghost_blocks[i].vars[DIFF_CLk] = 7733;
    	w->ghost_blocks[i].vars[DIFF_CAk] = 0.1612;
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
			// Allocate arrays for sending and receiving different variables (currently 2, change manually)
			double *send_array_Ke = malloc(sizeof(double) * side_length);
			double *send_array_Nae = malloc(sizeof(double) * side_length);
			double *send_array_Kk = malloc(sizeof(double) * side_length);
			double *send_array_vk = malloc(sizeof(double) * side_length);
			double *send_array_Nak = malloc(sizeof(double) * side_length);
			double *send_array_HCO3k = malloc(sizeof(double) * side_length);
			double *send_array_Clk = malloc(sizeof(double) * side_length);
			double *send_array_Cak = malloc(sizeof(double) * side_length);
			
			double *receive_array_Ke = malloc(sizeof(double) * side_length);
			double *receive_array_Nae = malloc(sizeof(double) * side_length);
			double *receive_array_Kk = malloc(sizeof(double) * side_length);
			double *receive_array_vk = malloc(sizeof(double) * side_length);
			double *receive_array_Nak = malloc(sizeof(double) * side_length);
			double *receive_array_HCO3k = malloc(sizeof(double) * side_length);
			double *receive_array_Clk = malloc(sizeof(double) * side_length);
			double *receive_array_Cak = malloc(sizeof(double) * side_length);

			// Populate send arrays with values from the state variables array.
			for(int i = 0; i < side_length; i++)
			{
				// Get the state variables based on the indices of the edge block for this side (and all other sides).
				send_array_Ke[i] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + i_K_e];
				send_array_Nae[i] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + i_Na_e];
				send_array_Kk[i] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + i_K_k];
				send_array_vk[i] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + i_v_k];
				send_array_Nak[i] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + i_Na_k];
				send_array_HCO3k[i] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + i_HCO3_k];
				send_array_Clk[i] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + i_Cl_k];
				send_array_Cak[i] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + i_ca_k];
			}

			// Simultaneously send and receive the values for the ghost blocks.
			MPI_Sendrecv(
					send_array_Ke, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					receive_array_Ke, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Sendrecv(
					send_array_Nae, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					receive_array_Nae, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Sendrecv(
					send_array_Kk, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					receive_array_Kk, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Sendrecv(
					send_array_vk, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					receive_array_vk, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Sendrecv(
					send_array_Nak, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					receive_array_Nak, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Sendrecv(
					send_array_HCO3k, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					receive_array_HCO3k, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Sendrecv(
					send_array_Clk, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					receive_array_Clk, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Sendrecv(
					send_array_Cak, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					receive_array_Cak, side_length, MPI_DOUBLE, W->domain_neighbours[side_idx], 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			// Copy received data into the ghost blocks.
			for(int i = 0; i < side_length; i++)
			{
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_K] = receive_array_Ke[i];
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_NA] = receive_array_Nae[i];
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_Kk] = receive_array_Kk[i];
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_Vk] = receive_array_vk[i];
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_NAk] = receive_array_Nak[i];
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_HCO3k] = receive_array_HCO3k[i];
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_CLk] = receive_array_Clk[i];
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_CAk] = receive_array_Cak[i];
			}

			free(send_array_Ke);
			free(receive_array_Ke);
			free(send_array_Nae);
			free(receive_array_Nae);
			free(send_array_Kk);
			free(receive_array_Kk);
			free(send_array_vk);
			free(receive_array_vk);
			free(send_array_Nak);
			free(receive_array_Nak);
			free(send_array_HCO3k);
			free(receive_array_HCO3k);
			free(send_array_Clk);
			free(receive_array_Clk);
			free(send_array_Cak);
			free(receive_array_Cak);
		}
		// If there is no neighbour, copy the current block's value into the ghost block.
		else
		{
			for(int i = 0; i < side_length; i++)
			{
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_K] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + i_K_e];
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_NA] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + i_Na_e];
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_Kk] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + i_K_k];
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_Vk] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + i_v_k];
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_NAk] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + i_Na_k];
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_HCO3k] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + i_HCO3_k];
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_CLk] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + i_Cl_k];
				W->nvu_w->ghost_blocks[idx_offset + i].vars[DIFF_CAk] = y[W->neq * W->nvu_w->edge_indices[idx_offset + i] + i_ca_k];
			}
		}
		idx_offset += side_length;
		odd_even = !odd_even;
	}
}
