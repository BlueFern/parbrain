#ifndef NVU_H
#define NVU_H

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#include <cs.h>

// TODO: make not external cos it's BAD!!! But check first

// External constant variables defined in brain.c.
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

    // Other NVU parameters for radius and pressure.
    double a1, a2, a3, a4, a5;
    double b1, d1, d2, g1, g2;
    double l;
//    double gamma, cstar;
    double pcap;

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

//time- and space-dependent K+ input
double K_input(double t, double x, double y);
	
//time- and space-dependent flux_ft input
double flux_ft(double t, double x, double y);
	
//time- and space-dependent PLC input
double PLC_input(double t, double x, double y);
	
//factorial
double factorial(int c);

// Initial conditions
void   nvu_ics(double *u0, double x, double y, nvu_workspace *nvu_w);


#endif
