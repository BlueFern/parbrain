#ifndef NVU_H
#define NVU_H

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#include <cs.h>

// Forward declaration to avoid errors
typedef struct ghost_block ghost_block;

// State variable indexing (TODO: shouldn't be initialised here but is needed in various files - make external or something?)
static const int i_radius  = 0; // radius has to be 0, this is assumed elsewhere

// AC
static const int R_k       = 1;
static const int N_Na_k    = 2;
static const int N_K_k     = 3;
static const int N_HCO3_k  = 4;
static const int N_Cl_k    = 5;
static const int w_k       = 10;

// SC
static const int N_Na_s    = 6;
static const int N_K_s     = 7;
static const int N_HCO3_s  = 8;

// PVS
static const int K_p       = 9;

// SMC
static const int ca_i      = 11; // or use enum
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

// Neuron - ions
static const int v_sa	   = 24;
static const int v_d	   = 25;
static const int K_sa	   = 26;
static const int Na_sa	   = 27;
static const int Cl_sa	   = 28;
static const int K_d	   = 29;
static const int Na_d	   = 30;
static const int Cl_d	   = 31;
static const int K_e	   = 32;
static const int Na_e	   = 33;
static const int Cl_e	   = 34;

// Neuron - other
static const int Buff_e	   = 35;
static const int Buff_s	   = 36;
static const int O2		   = 37;
static const int CBV	   = 38;
static const int K_buff	   = 39;

// Neuron Gating Variables
static const int m1	   	   = 40;
static const int m2	   	   = 41;
static const int m3	   	   = 42;
static const int m4	   	   = 43;
static const int m5	   	   = 44;
static const int m6	   	   = 45;
static const int m7	   	   = 46;
static const int m8	   	   = 47;
static const int h1	   	   = 48;
static const int h2	   	   = 49;
static const int h3	   	   = 50;
static const int h4	   	   = 51;
static const int h5	   	   = 52;
static const int h6	   	   = 53;

// Other
static const int PLC_i     = 54; //!
static const int current_i    = 55; //!

// Constants we may want to use that are defined in brain.c. 
extern const double RMIN;   // radius of smallest vessel
extern const double R0;     // radius scaling characteristic value
extern const double P0;     // pressure characteristic value (Pa)
extern const double PCAP;   // nominal capillary bed pressure (Pa)

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

    // Other NVU parameters for radius and pressure. TODO: rename
    double a1, a2, a3, a4, a5;
    double b1, d1, d2, g1, g2;
    double l;
    double pcap;		// pressure at capillaries (min)

    // Indices of neighbours for every tissue block. Ghost blocks are numbered
    // with negative indices in counter-clockwise flu_Kactivation_idirection.
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
void   nvu_rhs(double t, double x, double y, double p, double *u, double *du, nvu_workspace *nvu_w);

// Tidy up routine. Free anything allocated in nvu_init here
void  *nvu_free(nvu_workspace *nvu_w);

// Time-varying input pressure function
double nvu_p0(double t);

//time- and space-dependent Glu input
//double nvu_Glu(double t, double x, double y);

//time- and space-dependent current input
double current_input(double t, double x, double y);
	
//time- and space-dependent PLC input
double PLC_input(double t, double x, double y);

//time- and space-dependent ECS K+ input
double ECS_input(double t, double x, double y);
	
//factorial
double factorial(int c);

// Initial conditions
void nvu_ics(double *u0, double x, double y, nvu_workspace *nvu_w);

#endif
