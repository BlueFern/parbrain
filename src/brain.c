#include "brain.h"
#include "diffusion.h"

const int NDEFAULT      = 9;
const int NSUBDEFAULT   = 3;
const int NSYMBOLS      = 4;

const double RMIN  = 20e-6 ;  							// m, radius of smallest vessel
const double BIFURCATION_SCALE = 1.4142135623730951;	// = sqrt(2), amount the radius decreases by when going down a level
const double L0    = 200e-6;  							// m (for nondimensionalising), length characteristic value
const double LRR   = 10    ;  							// Nondimensional, length to radius ratio
const double MU    = 3.5e-3;  							// Pa s, blood viscosity

// External variables (initialised in nvu.h, also used in nvu.c)
const double R0    = 20e-6 ;  // m (for nondimensionalising)
const double P0    = 8000  ;  // Pa (scaling factor for nondim)
const double PCAP  = 4000  ;  // Pa (capillary bed pressure)


/* Methods for external use: essentially function evaluation and Jacobian
 * updates. evaluate, solve and jacupdate both have the side effect of updating
 * conductance, pressure, and flow values in the workspace data structure,
 * so these values should not be relied upon to remain static
 * */

// Initialises the workspace
workspace * workspace_init(int argc, char **argv)
{
    workspace *W;
    W = malloc(sizeof *W);

    W->jacupdates = 0;
    W->fevals     = 0;
    W->QglobalPos = 0;
    W->PglobalPos = 0;

    init_parallel(W, argc, argv);   // Initialise splitting into subtrees and MPI stuff
    init_subtree(W);                // Init adjacency matrix and workspace for subtree
    init_roottree(W);               // Same, but for root tree

    set_spatial_coordinates(W);

    // Set the indices of the MPI domain neighbours.
    set_domain_neighbours(W->rank, W->mglobal, W->nglobal, W->domain_neighbours);

	compute_symbol_cholesky(W);            		// Precompute symbolic factorisations
    W->nvu_w = nvu_init();           				// Initialise ODE parameter workspace

    set_block_neighbours(W->nlocal, W->mlocal, W->nvu_w); // Calculate neighbours indices for each block.
    set_edge_indices(W->nlocal, W->mlocal, W->nvu_w); 	// Calculate the indices of all four edges.
    init_ghost_blocks(W->nlocal, W->mlocal, W->nvu_w); 	// Initialise the ghost block structs.

    W->neq = W->nvu_w->neq;
    W->nu  = W->neq * W->nblocks;   // no of state variables per rank
    set_conductance(W, 0, 1);       // set scaled conductances
    set_length(W);                  // Initialise the vessel lengths
    init_jacobians(W);              // Initialise Jacobian data structures
    init_io(W);                     // Initialise output files
    write_info(W);                  // Write summary information to disk

    return W;
}

void evaluate(workspace *W, double t, double *y, double *dy)
{
    W->fevals++;
    double r, l;

    // printf("time: %f", t);

    // Update conductances of autoregulating vessels
    for (int i = 0; i < W->nblocks; i++)
    {
        r = y[W->neq*i]; // Radius is always the first variable
        l = W->l[i];
        W->g[i] = pow(r, 4) / l;
    }
    // Solve for pressure and flow: takes p0 and pcap (boundary conditions) as input
    solve(W, nvu_p0(t), W->nvu_w->pcap);
    // Evaluate the right hand side equations
    rhs(W, t, y, W->p, dy);

    update_ghost_blocks(W, y);

    // Calculate diffusion for every block.
    int istart;
    for (int i = 0; i < W->nblocks; i++)
    {
    	istart = W->neq * i;
    	diffusion(i, t, y + istart, dy + istart, W->nvu_w);
    }

    // printf(", again, time: %f\n", t);
}

void solve(workspace *W, double p0, double pcap)
{
    compute_uv(W, pcap); // compute the intermediate variables u and v
    if (W->n_procs > 1)
    {
        communicate(W);
        compute_root(W, p0);
    }
    compute_sub(W, p0, pcap); // compute p, w and q
}

void jacupdate(workspace *W, double t, double *u)
{
    W->jacupdates++;
    double *f;
    double eps = 1e-6;

    // Evaluate the right hand side
    f = malloc(W->nu * sizeof (*f ));
    evaluate(W, t, u, f);
    //rhs(W, t, y, W->p, f);

    // This order is important. dpdg depends on dgdx, and on the correct
    // value of W->w from the call to evaluate. dfdx and dfdp will modify
    // this value
    eval_dgdx(W, t, u);
    eval_dpdg(W, t, u);
    eval_dfdx(W, t, u, f, eps);
    eval_dfdp(W, t, u, f, eps); 

    if (W->isjac) cs_spfree(W->J);
    cs *Y;
    // Compute the product dfdx + dfdp dpdg dgdx

    Y = cs_multiply(W->dfdp->A, W->dpdgneg);
    W->J = cs_add(W->dfdx->A, Y, 1, -1);
    W->isjac = 1;
    cs_spfree(Y);

    free(f);
}

/* Private functions */
void init_parallel(workspace *W, int argc, char **argv)
{
    // Sort out MPI configuration
    MPI_Comm_size(MPI_COMM_WORLD, &W->n_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &W->rank);
    assert(is_power_of_two(W->n_procs));

    /*
    Parse input parameters and set tree sizes. There are four N values
    that matter: 
        N is the number of levels in the tree (total)
        N0: number of levels in the root subtree
        Np: number of levels in the subtrees corresponding to each
        core (N0 + Np = N)
        Nsub: number of levels in the small scale subtrees for Jacobian
        computation
    */
    W->N    = NDEFAULT; 
    W->Nsub = NSUBDEFAULT;

    if (argc > 1)
    {
        W->N = atoi(argv[1]); // N has been specified at command line
    }
    if (argc > 2)
    {
    	W->Nsub = atoi(argv[2]); // Nsub has been specified at command line
    }

    W->N0 = (int) round(log2((double) W->n_procs)); 
    W->Np = W->N - W->N0; 		// Dependent on number of levels and number of cores running


    /* Check that the user input etc. makes sense. This catches the two bad
     * cases of too many workers, or Nsub being too large */
    assert(W->Np > W->Nsub);


    // Configure buffer and parameters for MPI_Allgather
    W->buf  = malloc(W->n_procs * NSYMBOLS * sizeof(*W->buf));
    W->flag = malloc(W->n_procs * sizeof (*W->flag));
    W->n_writes = 0;
}

void init_io(workspace *W)
{
    int sizes[2];
    int subsizes[2];
    int starts[2];
    int res=0;
    int i=0;

    // Initialise files for MPI I/O. Requires init_parallel to have been
    // called first

    char tSuffix[] = "/tissueBlocks.dat";
    char qSuffix[] = "/flow.dat";
    char pSuffix[] = "/pressure.dat";

    W->dirName = malloc(FILENAMESIZE/2 * sizeof(*W->dirName));

    if (W->rank == 0)
    {
        sprintf(W->dirName, "np%02d_nlev%02d_sbtr%02d", W->n_procs, W->N, W->Nsub);

        res = mkdir(W->dirName, S_IRWXU | S_IRWXG | S_IRWXO); //res == 0 then the directory doesn't exist

        //if dirName already exists, add a suffix (eg. _1, _2 etc)
        while (res == -1)
		{
			i++;
			sprintf(W->dirName, "np%02d_nlev%02d_sbtr%02d_%d", W->n_procs, W->N, W->Nsub, i);
			res = mkdir(W->dirName, S_IRWXU | S_IRWXG | S_IRWXO);
		}
    }

    MPI_Bcast(W->dirName, FILENAMESIZE/2, MPI_CHAR, 0, MPI_COMM_WORLD);

    W->Toutfilename = malloc(FILENAMESIZE * sizeof(*W->Toutfilename));
    sprintf(W->Toutfilename, "%s%s", W->dirName, tSuffix);
    MPI_File_open(MPI_COMM_WORLD, W->Toutfilename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &W->Toutfile);

    W->Qoutfilename = malloc(FILENAMESIZE*sizeof(*W->Qoutfilename));
    sprintf(W->Qoutfilename, "%s%s", W->dirName, qSuffix);
    MPI_File_open(MPI_COMM_WORLD, W->Qoutfilename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &W->Qoutfile);

    W->Poutfilename = malloc(FILENAMESIZE*sizeof(*W->Poutfilename));
    sprintf(W->Poutfilename, "%s%s",W->dirName,pSuffix);
    MPI_File_open(MPI_COMM_WORLD, W->Poutfilename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &W->Poutfile);

    // Create subarray data type. This assumes column major ordering
    // (MPI_ORDER_FORTRAN). The factor of W->neq should be moved to
    // subsizes[1] if row major (MPI_ORDER_C) is desired
    subsizes[0] = W->mlocal * W->neq;
    subsizes[1] = W->nlocal;

    sizes[0] = subsizes[0] * W->mglobal;
    sizes[1] = subsizes[1] * W->nglobal;

    starts[0] = subsizes[0] * (W->rank % W->mglobal);
    starts[1] = subsizes[1] * (W->rank / W->mglobal);

    // Create subarray for writing state variables
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE, &W->subarray);
    MPI_Type_commit(&W->subarray);

    // Create subarray for writing single double per block
    subsizes[0] = W->mlocal;
    sizes[0] = subsizes[0] * W->mglobal;
    starts[0] = subsizes[0] * (W->rank % W->mglobal);
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE, &W->subarray_single);
    MPI_Type_commit(&W->subarray_single);

    W->n_writes = 0;
    W->displacement = 0; // bytes
    W->displacement_per_write = (sizeof(double)) * (1 + W->nu * W->n_procs); // bytes 

}

void close_io(workspace *W)
{
    // Close data files
    MPI_Type_free(&W->subarray);
    MPI_File_close(&W->Toutfile);
    MPI_File_close(&W->Qoutfile);
    MPI_File_close(&W->Poutfile);
    free(W->Toutfilename);
    free(W->Qoutfilename);
    free(W->Poutfilename);
}

void write_data(workspace *W, double t, double *y)
{
    // Set view for writing of timestamps
    MPI_File_set_view(W->Toutfile, W->displacement, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

    if (W->rank == 0)
    {
        MPI_File_write(W->Toutfile, &t, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_set_view(W->Toutfile, W->displacement + sizeof(t), MPI_DOUBLE, W->subarray, "native", MPI_INFO_NULL);
    MPI_File_write_all(W->Toutfile, y, W->nu, MPI_DOUBLE, MPI_STATUS_IGNORE);

    W->displacement += W->displacement_per_write;
    W->n_writes++;
}

void write_flow(workspace *W, double t, double *q, double *q0)
{
    int displ0, displ1;
   
    int xbranch = 1;  
    int pos = W->QglobalPos;  	// rank-specific position in Qtot
    int pos_q = 0;     			// position in q vector

    int nlocal = W->nlocal; 		// because the values will get updated here
    int mlocal = W->mlocal;
    int nglobal = W->nglobal;
    int mglobal = W->mglobal;

    MPI_File_set_view(W->Qoutfile, pos*sizeof(double), MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

    if (W->rank == 0)
    {
        MPI_File_write(W->Qoutfile, &t, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }

    pos += 1; //increase due to timestamp written above
    MPI_Barrier(MPI_COMM_WORLD); //Just in case, leave the barrier here. May not be necessary.

    // Subtrees
    for (int level = 0; level < W->Np; level++)
    {
        //skip all elements until we reach the portion of data we are interested in
        displ0 = (W->rank/mglobal) * mglobal * mlocal * nlocal + (W->rank % mglobal) * mlocal;
        displ1 = 0;

        for (int j = 0; j < nlocal; j++)
        {
            MPI_File_set_view(W->Qoutfile, (pos + displ0 + displ1)*sizeof(double), MPI_DOUBLE, MPI_DOUBLE,"native", MPI_INFO_NULL);
            MPI_File_write_all(W->Qoutfile, &q[pos_q], mlocal, MPI_DOUBLE, MPI_STATUS_IGNORE);
            pos_q = pos_q + mlocal;
            displ1 += mglobal * mlocal; //to jump to the next chunk of data in this level
        }
        pos += mglobal * nglobal * mlocal * nlocal;

        if (xbranch)
            mlocal = mlocal / 2;
        else
            nlocal = nlocal / 2;

        xbranch = !xbranch;
    }

    // Write data from root of tree

//    MPI_Barrier(MPI_COMM_WORLD); //Just in case, leave the barrier here. May not be necessary.
    MPI_File_set_view(W->Qoutfile, pos*sizeof(double), MPI_DOUBLE, MPI_DOUBLE,"native", MPI_INFO_NULL);
    int roottreeSize = POW_OF_2(W->N0) - 1;
    if (W->rank == 0)
    {
        MPI_File_write(W->Qoutfile, &q0[0], roottreeSize, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }
    pos += roottreeSize;
    W->QglobalPos = pos;
}

void write_pressure(workspace *W, double t, double *p, double *p0)
{
    int displ0, displ1;
   
    int xbranch = 0;  // different from q, because we don't have leaf level
    int pos = W->PglobalPos;   // rank-specific position in Ptot
    int pos_p = 0;     // position in p vector 

    int nlocal = W->nlocal;
    int mlocal_half = W->mlocal/2;  // half as many rows, because p only holds internal node pressure values
    int nglobal = W->nglobal;
    int mglobal = W->mglobal;

    MPI_File_set_view(W->Poutfile, pos*sizeof(double), MPI_DOUBLE, MPI_DOUBLE,"native", MPI_INFO_NULL);
    if (W->rank == 0)
    {
        MPI_File_write(W->Poutfile, &t, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }

    pos += 1; //increase due to timestamp written above
    MPI_Barrier(MPI_COMM_WORLD); //Just in case, leave the barrier here. May not be necessary.

    // Subtrees
    for (int level = 0; level < W->Np; level++)
    {
        displ0 = (W->rank/mglobal) * mglobal * nlocal * mlocal_half + (W->rank % mglobal) * mlocal_half; //skip all elements until we reach the portion of data we are interested in
        displ1 = 0;

        for (int j = 0; j < nlocal; j++)
        {
            MPI_File_set_view(W->Poutfile, (pos + displ0 + displ1)*sizeof(double), MPI_DOUBLE, MPI_DOUBLE,"native", MPI_INFO_NULL);
            MPI_File_write_all(W->Poutfile, &p[pos_p], mlocal_half, MPI_DOUBLE, MPI_STATUS_IGNORE);
            pos_p = pos_p + mlocal_half;
            displ1 += mglobal * mlocal_half; //to jump to the next chunk of data in this level
        }
    pos += mglobal * nglobal * mlocal_half * nlocal;

		if (xbranch)
			mlocal_half = mlocal_half / 2;
		else
			nlocal = nlocal / 2;

	xbranch = !xbranch;
    }
    
    // write data from roottree
    MPI_File_set_view(W->Poutfile, pos*sizeof(double), MPI_DOUBLE, MPI_DOUBLE,"native", MPI_INFO_NULL);
    int roottreeSize = POW_OF_2(W->N0) - 1;
    if (W->rank == 0)
    {
        MPI_File_write(W->Poutfile, &p0[0], roottreeSize, MPI_DOUBLE, MPI_STATUS_IGNORE);
        
    }
    pos = pos+ roottreeSize;
    
    W->PglobalPos = pos;
}

void write_info(workspace *W)
{
    // Write the summary info to disk
    char * infofilename;
    if (W->rank == 0)
    {
        FILE *fp;
        // Write the data file
		char iSuffix[] = "/info.dat";
		infofilename = malloc(FILENAMESIZE * sizeof(*infofilename));
        sprintf(infofilename, "%s%s", W->dirName, iSuffix);

        fp = fopen(infofilename, "w");
        fprintf(fp, "n_processors    n_blocks_per_rank        n_state_vars   m_local         n_local         m_global        n_global\n");
        fprintf(fp, "%-16d", W->n_procs);
        fprintf(fp, "%-16d", W->nblocks);
        fprintf(fp, "%-16d", W->neq);
        fprintf(fp, "%-16d", W->mlocal);
        fprintf(fp, "%-16d", W->nlocal);
        fprintf(fp, "%-16d", W->mglobal);
        fprintf(fp, "%-16d", W->nglobal);
        fprintf(fp, "\n");
        fclose(fp);

        free(infofilename);
    }

    // Write the x and y coordinates to the file
    MPI_File_set_view(W->Toutfile, W->displacement, MPI_DOUBLE, W->subarray_single, "native", MPI_INFO_NULL);
    MPI_File_write_all(W->Toutfile, W->x, W->nblocks, MPI_DOUBLE, MPI_STATUS_IGNORE);

    W->displacement += sizeof(*W->x) * W->nblocks * W->n_procs;

    MPI_File_set_view(W->Toutfile, W->displacement, MPI_DOUBLE, W->subarray_single, "native", MPI_INFO_NULL);
    MPI_File_write_all(W->Toutfile, W->y, W->nblocks, MPI_DOUBLE, MPI_STATUS_IGNORE);

    W->displacement += sizeof(*W->y) * W->nblocks * W->n_procs;

    //MPI_File_write(W->Toutfile, W->x, W->nblocks, MPI_DOUBLE, MPI_STATUS_IGNORE);
    //MPI_File_write(W->Toutfile, W->y, W->nblocks, MPI_DOUBLE, MPI_STATUS_IGNORE);
}

void set_spatial_coordinates(workspace *W)
{
    // work out some stuff
    double procs = (double) W->n_procs;
    int log2P = (int) log2(procs);
    int mlocal, nlocal, mglobal, nglobal;
    int iglobal, jglobal;

    int m, n; // number of rows / cols of blocks globally
    // if rectangular, make it so there are more rows than columns (i.e. if N is even, then N-1 % 2 = 1 and m is bigger)
    m = POW_OF_2((W->N - 1)/2 + (W->N - 1) %2);
    n = POW_OF_2((W->N - 1)/2);

    // Work out arrangement of workers, again, if rectangular, set more rows than columns
    if (W->N % 2 == 1) // odd number of levels - square
    {
        mglobal = POW_OF_2(log2P / 2);
        nglobal = POW_OF_2( (log2P / 2) + log2P % 2);
    }
    else // even number of levels - rectangular
    {
        mglobal = POW_OF_2( (log2P / 2) + log2P % 2);
        nglobal = POW_OF_2(log2P / 2);
    }
    
    // Work out how many rows / columns of blocks we have for each core
    mlocal = m / mglobal;
    nlocal = n / nglobal;

    iglobal = W->rank % mglobal;
    jglobal = W->rank / mglobal;


    double xoffset, yoffset;
    xoffset = (double) ((2*jglobal - (nglobal-1)) * nlocal) * L0;
    yoffset = (double) ((2*iglobal - (mglobal-1)) * mlocal) * L0;

    for (int j = 0; j < nlocal; j++)
    {
        for (int i = 0; i < mlocal; i++)
        {
            W->x[i + mlocal * j] = xoffset + L0 * (double) (2*j - (nlocal - 1));
            W->y[i + mlocal * j] = yoffset + L0 * (double) (2*i - (mlocal - 1));
        }
    }

    W->mlocal = mlocal;
    W->nlocal = nlocal;
    W->mglobal = mglobal;
    W->nglobal = nglobal;
}
 
int is_power_of_two (unsigned int x)
{
      return ((x != 0) && !(x & (x - 1)));	// &: bitwise AND, numbers that are powers of two are of binary form 10...0
}


// TODO: Is there function to free the subtree somewhere?

void init_subtree(workspace *W)
{
    // First construct the local subtree 
    W->A    = adjacency(W->Np); // where Np = N - N0, dependent on the number of levels (N) and cores (2^N0)
    W->At   = cs_transpose(W->A, 1);
    W->G    = speye(W->A->n); 	// creates identity matrix I_n
    W->g    = W->G->x;			// numerical values of G
    W->level = malloc(W->A->n * sizeof (*W->level));

    for (int i = 0; i < W->A->n; i++)
    {
        W->level[i] = (int) floor(log2(W->A->n - i)) + W->N0;
    }

    W->nblocks = W->A->m + 1;

    // Initialise workspace variables for solving
    W->l = malloc (W->nblocks * (sizeof *W->l));
    W->x = malloc (W->nblocks * (sizeof *W->x));
    W->y = malloc (W->nblocks * (sizeof *W->y));
    W->b = malloc (W->A->n * (sizeof *W->b));
    W->u = malloc (W->A->m * (sizeof *W->u));
    W->v = malloc (W->A->m * (sizeof *W->v));
    W->p = malloc (W->A->m * (sizeof *W->p));
    W->q = malloc (W->A->n * (sizeof *W->q));
    W->w = malloc (W->A->n * (sizeof *W->w));
    W->ucomm = malloc (W->n_procs * (sizeof *W->ucomm));
    W->vcomm = malloc (W->n_procs * (sizeof *W->vcomm));
    W->gcomm = malloc (W->n_procs * (sizeof *W->gcomm));


    // Initialise general purpose vectors for workspaces.
    // These vectors seem to be sitting there without use.
    W->xn = malloc(W->A->n * sizeof(*W->xn));
    W->xm = malloc(W->A->m * sizeof(*W->xm));
}

void compute_symbol_cholesky(workspace *W)
{
    cs *X, *Y;
    X = cs_multiply(W->A, W->G);
    Y = cs_multiply(X, W->At);
    cs_spfree(X);

    W->symbchol = cs_schol(0, Y);  // symbchol = AGAt
    cs_spfree(Y);

    X = cs_multiply(W->A0, W->G0);
    Y = cs_multiply(X, W->A0t);
    cs_spfree(X);

    W->symbchol0 = cs_schol(0, Y);	// symbchol0 = AGAt for the root subtree
    cs_spfree(Y);
}

void init_roottree(workspace *W)
{
    cs *T;
    int *ikeep;

    // Construct the full adjacency matrix and then remove the first m
    // columns
    T = adjacency(W->N0 + 1);
    ikeep = malloc(T->m * (sizeof *ikeep));

    for (int i = 0; i < T->m; i++)
    {
        ikeep[i] = T->n - T->m + i;
    }

    W->A0 = subsref(T, NULL, ikeep, -1, T->m);
    cs_spfree(T);

    W->A0t    = cs_transpose(W->A0, 1);
    W->G0     = speye(W->A0->n);
    W->g0     = W->G0->x;
    W->level0 = malloc(W->A->n * sizeof (*W->level0));

    for (int i = 0; i < W->A0->n; i++)
    {
        W->level0[i] = (int) floor(log2(W->A0->n - i));
    }

    // Initialise workspace variables for solving 
    W->b0  = malloc (W->A0->n * sizeof (*W->b0));
    W->p0  = malloc (W->A0->m * sizeof (*W->b0));
    W->q0  = malloc (W->A0->n * sizeof (*W->b0));
    W->xn0 = malloc (W->A0->n * sizeof (*W->xn0));

    free(ikeep);
}


void set_conductance(workspace *W, int unscaled, int computeroot)
{
    // if unscaled is true, we can compute the conductances for an unscaled
    // version of the problem.
    double r, l;

    if (unscaled)
    {
        for (int i = 0; i < W->A->n; i++)
        {
            r = compute_radius(W->level[i], W->N);
            l = compute_length(W->level[i], W->N);
            W->g[i] = M_PI * pow(r, 4) / (8.0 * MU * l);
        }
    }
    else
    {
        for (int i = 0; i < W->A->n; i++)
        {
            r = compute_radius(W->level[i], W->N) / R0;
            l = compute_length(W->level[i], W->N) / L0;
            W->g[i] = pow(r, 4) / l;
        }
    }

    if (computeroot)
    {
        if (unscaled)
        {
            for (int i = 0; i < W->A0->n; i++)
            {
                r = compute_radius(W->level0[i], W->N);
                l = compute_length(W->level0[i], W->N);
                W->g0[i] = M_PI * pow(r, 4) / (8.0 * MU * l);
            }
        }
        else
        {
            for (int i = 0; i < W->A0->n; i++)
            {
                r = compute_radius(W->level0[i], W->N) / R0;
                l = compute_length(W->level0[i], W->N) / L0;
                W->g0[i] = pow(r, 4) / l;
            }
        }
    }
}

void set_length(workspace *W)
{
    // Set lengths of autoregulating vessels
    for (int i = 0; i < W->nblocks; i++)
    {
        W->l[i] = compute_length(W->level[i], W->N) / L0;
    }
}

double compute_length(int level, int n_levels)
{
    double length;
    length = LRR * RMIN * (double) POW_OF_2((n_levels - level - 1)/ 2);
    return length;
}

double compute_radius(int level, int n_levels)
{
    double radius;
    //r = RMIN * pow(2., ((double) (n_levels - level - 1)) / 2.);
    radius = RMIN * (double) POW_OF_2((n_levels - level - 1)/ 2);
    return radius;
}

void init_jacobians(workspace *W)
{
    // Create the data structures for each of the Jacobians
    W->isjac = 0;
    init_dgdx(W);
    init_dpdg(W);
    init_dfdx(W);
    init_dfdp(W);
}

void init_dgdx(workspace *W)
{
    //conductance depends on the vessel scaled length and radius only. This
    //function simply sets up the n * (nblocks * neq)

	// standard triplet matrix requirements
    cs *T;
    int *Ti, *Tj;
    double *Tx;

    int neq = W->nvu_w->neq; // number of equations per block

    // matrix is 
    T = cs_spalloc(W->A->n, neq * W->nblocks, W->nblocks, 1, 1);
    Ti = T->i;
    Tj = T->p;
    Tx = T->x;

    for (int i = 0; i < W->nblocks; i++)
    {
        Ti[i] = i;
        Tj[i] = neq*i; // radius term is always the first in the block
        Tx[i] = 1.0;
    }

    T->nz = W->nblocks;
    W->dgdx = cs_compress(T);
    cs_spfree(T);
}

void init_dpdg(workspace *W)
{
    // In this function we'll figure out A_A, etc, and
    // compute the symbolic factorisation of A_A G A_A^T

    // Compute the indices B of the nodes to remove from A
    int Nr = W->Np - W->Nsub;
    int m = W->A->m;
    int B0 = POW_OF_2(W->Np - 1) - POW_OF_2(Nr);
    int B1 = POW_OF_2(W->Np - 1) - POW_OF_2(Nr - 1);

    // Construct the temporary matrix to store the projector in
    cs *T;
    int *Ti, *Tj;
    double *Tx;

    T = cs_spalloc(m - (B1 - B0), m, m - (B1 - B0), 1, 1);
    Ti  = T->i;
    Tj = T->p;
    Tx = T->x;

    int k = 0;

    for (int i = 0; i < B0; i++)
    {
        Ti[k] = k;
        Tj[k] = i;
        Tx[k++] = 1.0;
    }

    for (int i = B1; i < m; i++)
    {
        Ti[k] = k;
        Tj[k] = i;
        Tx[k++] = 1.0;
    }

    T->nz = k;

    cs *Pt;
    Pt = cs_compress(T);
    cs_spfree(T);
    
    // Compute A_A, A_A.', and the matrix Proj used to construct dp/dg
    W->A_A  = cs_multiply(Pt, W->A);
    W->A_At = cs_transpose(W->A_A, 1);
    W->Proj = cs_transpose(Pt, 1);
    cs_spfree(Pt);

    // Compute symbolic Cholesky factorisation of A_A G A_A^T
    T = cs_multiply(W->A_A, W->A_At);
    W->symbchol_reduced = cs_schol(0, T);
    cs_spfree(T);
}

void init_dfdx(workspace *W)
{
    // Load sparsity pattern for one block (dfdx_pattern) and add to make Jacobian for all blocks (dfdx)
    int nblocks = W->nblocks;
    cs *J;

    J = blkdiag(W->nvu_w->dfdx_pattern, nblocks, nblocks);
    W->dfdx = numjacinit(J);

    cs_spfree(J);
}

void init_dfdp(workspace *W)
{
	// Creates large Jacobian for all NVU blocks - essentially the dfdp_pattern added together heaps
    // The ordering of the blocks is such that the first two see the first
    // pressure, the second two see the second, etc.
    cs *B, *J;
    B = vertcat(W->nvu_w->dfdp_pattern, W->nvu_w->dfdp_pattern);
    J = blkdiag(B, W->nblocks / 2, W->A->m);

    W->dfdp = numjacinit(J);
    
    cs_spfree(J);
    cs_spfree(B);
}

void compute_uv(workspace *W, double pcap)
{
	// u = -AGb, v = [g,0]
    cs *AG, *B;
    csn *Nu;

    // Set up boundary conditions b - put them in q for now
    for (int i = 0; i < W->A->m + 1; i++)
    {
        W->q[i] = -pcap;
    }

    for (int i = W->A->m + 1; i < W->A->n; i++)
    {
    	W->q[i] = 0;
    }

    // Solving A*G*At p = -A*G b (b is in q for now)
    // Define the matrices AG = A*G, and B = A*G*At
    AG = cs_multiply(W->A, W->G);
    B  = cs_multiply(AG, W->At);

    // Numerical Cholesky factorisation using precomputed symbolic, symbchol = AGAt. Maybe replaces symbchol with B?
    Nu = cs_chol(B, W->symbchol);

    if (!Nu) printf("Numerical cholesky decomposition failed in compute_uv");

    cs_spfree(B); // B is (potentially) big, so free it

    // Define the RHS of the u equations
    for (int i = 0; i < W->A->m; i++)
    {
        W->u[i] = 0.0;                       /* u = 0 */
    }

    cs_gaxpy(AG, W->q, W->u);               /* u = AGb + u = AGb */

    for (int i = 0; i < W->A->m; i++)
    {
        W->u[i] = -W->u[i];                 /* u = -u = -AGb */
    }

    cs_spfree(AG); 

    // And solve, using our computed factorisations - what is this outputting..?
    cholsoln(Nu, W->symbchol, W->A->m, W->u, W->xm);

    // Define RHS of v equations
    W->v[W->A->m-1] = W->g[W->A->n-1];

    for (int i = 0; i < W->A->m - 1; i++)
    {
        W->v[i] = 0.0;
    }

    cholsoln(Nu, W->symbchol, W->A->m, W->v, W->xm);
    cs_nfree(Nu); 

    // Put the end entries into communication buffers
    W->ucomm[W->rank] = W->u[W->A->m-1];
    W->vcomm[W->rank] = W->v[W->A->m-1];
    W->gcomm[W->rank] = W->g[W->A->n-1];
}

void communicate(workspace *W)
{
    double send_buf[NSYMBOLS], *recv_buf;

    recv_buf = W->buf;

    // Put local variables into the communication buffer
    send_buf[0] = W->ucomm[W->rank];
    send_buf[1] = W->vcomm[W->rank];
    send_buf[2] = W->gcomm[W->rank];
    send_buf[3] = (double) W->flag[W->rank];

    // Fill up the buffer with MPI all-to-all communication
    MPI_Allgather(send_buf, NSYMBOLS, MPI_DOUBLE, recv_buf, NSYMBOLS, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // Populate u1, v1, g1
    for (int i = 0; i < W->n_procs; i++)
    {
        if (i != W->rank)
        {
            W->ucomm[i] = recv_buf[NSYMBOLS * i    ];
            W->vcomm[i] = recv_buf[NSYMBOLS * i + 1];
            W->gcomm[i] = recv_buf[NSYMBOLS * i + 2];
            W->flag[i]  = (int) recv_buf[NSYMBOLS * i + 3]; 
        }
    }
}




void compute_root(workspace *W, double p0)
{
    cs *AG, *B, *X, *D;
    csn *Nu;
    double *d;

    d = W->xn0;

    for (int i = 0; i < W->A0->n; i++)
    {
        d[i] = 0.0;
    }

    // Define the diagonal modification D to be added to AGA'
    for (int i = 0; i < W->n_procs; i++)
    {
        d[i/2] += W->gcomm[i] * (1 - W->vcomm[i]); // Parent is simply i/2 
    }

    D = spdiags(d, W->A0->n);

    // Define the matrices AG = A_0 G_0 and B = A_0 G_0 A_0^T + D 
    AG = cs_multiply(W->A0, W->G0);
    X = cs_multiply(AG, W->A0t);
    B = cs_add(X, D, 1.0, 1.0);

    cs_spfree(D); 
    cs_spfree(X);

    // Perform numerical Cholesky factorisation 
    Nu = cs_chol(B, W->symbchol0);
    if (!Nu) printf("Numerical Cholesky failed in compute_root\n");
    cs_spfree(B);

    // Define the boundary condition b (stick it in q0 for now)
    W->q0[W->A0->n-1] = p0;
    for (int i = 0; i < (W->A0->n - 1); i++)
    {
        W->q0[i] = 0.0;
    }

    // Compute the right hand side of the p equation
    for (int i = 0; i < W->A0->m; i++)
    {
        W->p0[i] = 0.0;                              // p0 = 0
    }
    cs_gaxpy(AG, W->q0, W->p0);                     // p0 = AGb + p0 = AGb 

    for (int i = 0; i < W->A0->m; i++)
    {
        W->p0[i] = -W->p0[i];
    }
    // p = -p = -AGb */
    // p = -AGb + sum g_i u_i e_ki 

    for (int i = 0; i < W->n_procs; i++)
    {
    	W->p0[i/2] += W->gcomm[i] * W->ucomm[i];
    }

    cs_spfree(AG);

    // And solve, using our numerical Cholesky
    cholsoln(Nu, W->symbchol0, W->A0->m, W->p0, W->xn0);          
    // W->p0 is overwritten with p 

    // Now compute q0: W->q0 is currently p_0 e_n
    cs_gaxpy(W->A0t, W->p0, W->q0);       // q = A'p + q = A'p + p_0 en

    for (int i = 0; i < W->A0->m; i++)
    {
        W->q0[i] *= W->g0[i];
    }

    cs_nfree(Nu);
}

// The function that actually solves for p, w and q!
void compute_sub(workspace *W, double p0, double pcap)
{
    double pk;

    if (W->n_procs == 1) // the whole tree
    {
        pk = p0;
    }
    else // a sub tree
    {
        pk = W->p0[W->rank / 2];
    }

    for (int i = 0; i < W->A->m; i++)
    {
        W->p[i] = W->u[i] + pk * W->v[i];
    }

    // figure out w
    W->w[W->A->n - 1] = pk;

    // Set w to boundary conditions b
    for (int i = 0; i < W->A->m + 1 ; i++)
    {
        W->w[i] = -pcap;
    }
    for (int i = W->A->m + 1; i < (W->A->n - 1); i++)
    {
        W->w[i] = 0;
    }

    cs_gaxpy(W->At, W->p, W->w); // w = At*p + w (w set to b)

    for (int i = 0; i < W->A->n; i++)
    {
        W->q[i] = W->w[i] * W->g[i]; // q = w g
    }
}

void eval_dgdx(workspace *W, double t, double *y)
{
    double r, l;

    for (int i = 0; i < W->nblocks; i++)
    {
        r = y[i * W->neq];
        l = W->l[i];
        W->dgdx->x[i] = 4.0*pow(r, 3) / l;		// g = r^4 / l so dg/dx = 4 r^3 / l
    }
}

void eval_dpdg(workspace *W, double t, double *y)
{
    // Free the old one
    if (W->isjac) cs_spfree(W->dpdgneg);

    // Update the conductance matrix
    double r;

    for (int i = 0; i < W->nblocks; i++)
    {
        r = y[W->neq*i];
        W->g[i] = pow(r, 4) / W->l[i];			// g = r^4 / l
    }

    // Form the matrix A_A G A_A.'
    cs *X, *Y;
    Y = cs_multiply(W->A_A, W->G);
    X = cs_multiply(Y, W->A_At);

    cs_spfree(Y);

    // Form the right hand side matrix
    cs *D, *B, *C;

    D = spdiags(W->w, W->A->n); // w is set by the last evaluate
    C = cs_multiply(D, W->dgdx);
    cs_spfree(D);

    B = cs_multiply(W->A_A, C);
    cs_spfree(C);

    // Solve to form dpdg_A
    cs *PA;

    PA = mldivide_chol(X, W->symbchol_reduced, B);
    W->dpdgneg = cs_multiply(W->Proj, PA);

    cs_spfree(X);
    cs_spfree(PA);
    cs_spfree(B);
}

void eval_dfdx(workspace *W, double t, double *y, double *f, double eps)
{
    int i, j;
    double *y1, *h, *f1;

    y1 = malloc(W->nu * sizeof (*y1));
    h  = malloc(W->nu * sizeof (*h));
    f1 = malloc(W->nu * sizeof (*f1));
    
    for (int igrp = 0; igrp < W->dfdx->ng; igrp++)
    {
        for (int k = 0; k < W->nu; k++)
        {
            y1[k] = y[k];
        }

        for (int k = W->dfdx->r[igrp]; k < W->dfdx->r[igrp + 1]; k++)
        {
            j = W->dfdx->g[k];
            h[j] = eps; // * fabs(y[j]); 
            y1[j] += h[j];
        }

        rhs(W, t, y1, W->p, f1);  // f(t, y + dy, p)

        for (int k = W->dfdx->r[igrp]; k < W->dfdx->r[igrp+1]; k++)
        {
            j = W->dfdx->g[k];

            for (int ip = W->dfdx->A->p[j]; ip < W->dfdx->A->p[j+1]; ip++)
            {
                i = W->dfdx->A->i[ip];
                W->dfdx->A->x[ip] = (f1[i] - f[i]) / h[j];
            }
        }
    }

    free(y1); free(h); free(f1);
}

void eval_dfdp(workspace *W, double t, double *y, double *f, double eps)
{
    int i, j;
    double *h, *p1, *f1;

    h =  malloc(W->A->m  * sizeof (*h));
    p1 = malloc(W->A->m  * sizeof (*p1));
    f1 = malloc(W->nu * sizeof (*f1));

    for (int igrp = 0; igrp < W->dfdp->ng; igrp++)
    {
        // Set p1 back to p
        for (int k = 0; k < W->A->m; k++)
            p1[k] = W->p[k];

        // Increment entry of h for every entry column in the group 
        for (int k = W->dfdp->r[igrp]; k < W->dfdp->r[igrp+1]; k++)
        {
            j = W->dfdp->g[k];
            h[j] = eps; // * fabs(W->p[j]);
            p1[j] += h[j];
        }

        // Evaluate the right hand side
        rhs(W, t, y, p1, f1);

        // Iterate over the columns in the group
        for (int k = W->dfdp->r[igrp]; k < W->dfdp->r[igrp+1]; k++)
        {
            j = W->dfdp->g[k];

            // and then the rows in the column
            for (int ip = W->dfdp->A->p[j]; ip < W->dfdp->A->p[j+1]; ip++)
            {
                i = W->dfdp->A->i[ip];
                W->dfdp->A->x[ip] = (f1[i] - f[i]) / h[j];
            }
        }
    }

    free(h); free(f1); free(p1);
}

void rhs(workspace *W, double t, double *u, double *p, double *du)
{
    // Evaluate the right hand sides. Pressures etc have already been done
    int istart;
    for (int i = 0; i < W->nblocks; i++)
    {
        istart = W->neq * i;
        // Evaluate the individual right hand side
        nvu_rhs(t, W->x[i], W->y[i], p[i/2], u + istart, du + istart, W->nvu_w); // p[i/2]? p[0], p[1/2], p[1], ...? ints so round down - 0, 0, 1, 1 etc?
    }
}

void set_initial_conditions(workspace *W, double *u)
{
    int istart;
    for (int i = 0; i < W->nblocks; i++)
    {
        istart = W->neq * i;
        nvu_ics(u + istart, W->x[i], W->y[i], W->nvu_w);
    }
}

