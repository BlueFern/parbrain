#ifndef NVU_H
#define NVU_H
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#include <cs.h>

static const int i_radius  = 0; // radius has to be 0, this is assumed elsewhere
// AC
static const int R_k       = 1;
static const int N_Na_k    = 2;
static const int N_K_k     = 3;
static const int N_HCO3_k  = 4;
static const int N_Cl_k    = 5;
static const int N_Na_s    = 6;
static const int N_K_s     = 7;
static const int N_HCO3_s  = 8;
static const int K_p       = 9;
static const int w_k       = 10;
// SMC
static const int ca_i      = 11;
static const int ca_sr_i   = 12;
static const int v_i       = 13;
static const int w_i       = 14;
static const int ip3_i     = 15;
static const int K_i       = 16;
// EC
static const int ca_j      = 17;
static const int ca_er_j   = 18;
static const int v_j       = 19;
static const int ip3_j     = 20;
// Mech
static const int Mp        = 21;
static const int AMp       = 22;
static const int AM        = 23;
static const int PLC_i     = 24;
static const int K_df_i    = 25;
static const int K_flux_i  = 26;
// NO pathway
//static const int NOi       = 24;
//static const int NOj       = 25;
//static const int NOn       = 26;
//static const int cGMP      = 27;
//static const int eNOS      = 28;
//static const int nNOS      = 29;
//static const int ca_n      = 30;
//static const int E_b       = 31;
//static const int E_6c      = 32;
//static const int E_5c      = 33;

// Constants we may want to use that are defined in brain.c. 
extern const double RMIN;   // radius of smallest vessel
extern const double R0;     // radius scaling characteristic value
extern const double L0;     // length characteristic value
extern const double LRR;    // length to radius ratio
extern const double MU;     // blood viscosity (Pa s)
extern const double P0;     // pressure characteristic value (Pa)
extern const double PROOT;  // nominal tree root pressure (Pa)
extern const double PCAP;   // nominal capillary bed pressure (Pa)

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

// nvu_workspace gets created once per node (not once per block). So the
// parameters therein are generic for each nvu. Spatial inhomogeneity
// should be implemented by making the RHS explicitly dependent on the
// spatial coordinates of the block. 
//
// The struct needs to have the fields
//      neq: the number of differential equations per block
//      dfdp_pattern: a sparse neq x 1 matrix, where entries are 1 if the
//          equation depends on p (or q), but 0 otherwise
//      dfdx_pattern: sparse neq x neq matrix - the sparsity pattern of the
//          Jacobian of the block
typedef struct nvu_workspace {
    // Mandatory fields (accessed externally). Initialised in nvu_init
    int neq;
    cs *dfdp_pattern; // neq * 1 matrix indicating dependence on p
    cs *dfdx_pattern; // neq * neq matrix indicating Jacobian structure of nvu 

    // Other NVU parameters.
    double a1, a2, a3, a4, a5;
    double b1, d1, d2, g1, g2;
    double l;
    double gamma, cstar;
    double pcap;

    // Indices of neighbours for every tissue block. Ghost blocks are numbered
    // with negative indices in counter-clockwise direction.
    int *neighbours;

    // Indices of all edge tissue blocks for one MPI domain.
    int *edge_indices;

    // Ghost blocks for diffusion.
    // The number of ghost blocks is calculated on the basis of the number of blocks.
    int num_ghost_blocks;
    // Array of ghost blocks of size num_ghost_blocks.
    ghost_block *ghost_blocks;

} nvu_workspace;

// Initialisation routine. Gets called once before simulation
nvu_workspace* nvu_init(void); 

// Right hand side routine for one block
void   nvu_rhs(double t, double x, double y, double p, double *u, double *du, nvu_workspace *w);

// Tidy up routine. Free anything allocated in nvu_init here
void  *nvu_free(nvu_workspace *w); 

// Time-varying input pressure function
double nvu_p0(double t);

//time- and space-dependent Glu input
//double nvu_Glu(double t, double x, double y);

//time- and space-dependent K+ input
double K_input(double t, double x, double y);
	
//time- and space-dependent flux_ft input
double flux_ft(double t, double x, double y);
	
//time- and space-dependent PLC input
double PLC_input(double t, double x, double y);
	
//factorial
double factorial(int c);

// Initial conditions
void nvu_ics(double *u0, double x, double y, nvu_workspace *w);

// Get the indices for the four immediate neighbours of a given block.
void set_neighbours(int idx, int m, int n, int *neighbours);

// Get the indices for all neighbours for all tissue blocks in the given MPI domain.
void set_block_neighbours(int nlocal, int mlocal, nvu_workspace* w);

// Allocate the space for ghost blocks.
void init_ghost_blocks(int nlocal, int mlocal, nvu_workspace *w);

// Calculate diffusion for tissue blocks within given MPI domain.
void diffusion(int block_number, double t, double *u, double *du, nvu_workspace *w);


#endif
