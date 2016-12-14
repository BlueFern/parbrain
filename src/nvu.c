// possible mistakes in equations:
// / unitcon instead of *
// switch stretch channels on

#include "nvu.h"
#include <math.h>

// Pressure constants
static const double HRR        = 0.1   ;  // Nondimensional (thickness to radius ratio)
static const double RSCALE     = 0.6   ;  // Dimensionless
static const double E0         = 66e3  ;  // Pa
static const double EPASSIVE   = 66e3  ;  // Pa
static const double EACTIVE    = 233e3 ;  // Pa
static const double ETA        = 2.8e2 ;  // Pa s
static const double T0         = 1     ;  // s
static const double PA2MMHG    = 0.00750061683;

// nvu_init: this user-supplied function does any precomputation required
// for the model.
nvu_workspace *nvu_init(void)
{
    nvu_workspace *nvu_w;

    // Specify the sparsity patterns of the equations with respect to
    // pressure and state variables here. An equation is considered to be
    // dependent on pressure if it contains any pressure (transmural or
    // drop) or flow term
    int dfdp_pattern[] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 	// neq * 1 column vector

    // Column: variables the eq depends on
    // Row: variables that depend on that eq variable
    int dfdx_pattern[] = {1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0, 		//remember indexing starts at 0
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,0,0,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
						  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // column major, neq*neq

    // Initialise the workspace
    nvu_w = malloc(sizeof *nvu_w);
    nvu_w->neq = 56; 

    // Construct sparse matrices containing the sparsity patterns
    // TODO: modify dense2sparse so we can just use two calls to that,
    // rather than this messy code
    //
    // If you just define the integer arrays dfdp_pattern and dfdx_pattern
    // as above, you can leave the following two blocks as they are.

    // Takes the two patterns defined above and puts them into sparse matrices (cs*)
    cs *T;
    T = cs_spalloc(nvu_w->neq, 1, 1, 1, 1);

    for (int i = 0; i < nvu_w->neq; i++)
    {
        if (dfdp_pattern[i])
        {
            cs_entry(T, i, 0, 1.0);
        }
    }

    nvu_w->dfdp_pattern = cs_compress(T);
    cs_spfree(T);

    T = cs_spalloc(nvu_w->neq, nvu_w->neq, 1, 1, 1);

    for (int j = 0; j < nvu_w->neq; j++)
    {
        for (int i = 0; i < nvu_w->neq; i++)
        {
            if (dfdx_pattern[nvu_w->neq*j + i])
            {
                cs_entry(T, i, j, 1.0);
            }
        }
    }

    nvu_w->dfdx_pattern = cs_compress(T);
    cs_spfree(T);

    // Allocate other nvu workspace parameters
    double Rstar = R0;
    double hstar = HRR * Rstar;
    nvu_w->a1 = E0 * T0 * Rstar / (ETA * R0);
    nvu_w->a2 = P0 * Rstar * T0 / (ETA * hstar);
    nvu_w->a3 = Rstar / R0;
    nvu_w->a4 = 1 - RSCALE;
    nvu_w->a5 = EACTIVE / EPASSIVE - 1;
    nvu_w->pcap  = PCAP / P0;
    nvu_w->l  = 1; // normalised away

    return nvu_w;
}

// This frees the nvu_workspace structure. If you have allocated any memory
// within this struct, here is where you free it
void *nvu_free(nvu_workspace *nvu_w)
{
    cs_spfree(nvu_w->dfdp_pattern);
    cs_spfree(nvu_w->dfdx_pattern);
    free(nvu_w);
    return nvu_w;
}

// right hand side evaluation function. 
//      t       time, 
//      x,y     spatial coordinates, 
//      p       the pressure at the top of the vessel, 
//      u       state variables, the first of which is the vessel radius
//      du      output vector, in the same order (already allocated)
void nvu_rhs(double t, double x, double y, double p, double *u, double *du, nvu_workspace *nvu_w)
{

// general constants:
	const double F       	 = 96500         ;// [C mol-1] Faradays constant.
	const double R_gas       = 8.315         ;// [J mol-1K-1]
	const double Temp        = 300           ;// [K]
    const double unitcon     = 1e3           ;// [-] Factor to convert equations to another unit.

// NE & AC constants:
    const double L_p         = 2.1e-9        ;// [m uM-1s-1]
    const double R_tot       = 8.79e-8       ;// [m]   total volume surface area ratio AC+SC
    const double X_k         = 12.41e-3      ;// [uMm]
    const double z_Na        = 1             ;// [-]
    const double z_K         = 1             ;// [-]
    const double z_Cl        = -1            ;// [-]
    const double z_NBC       = -1            ;// [-]
    const double g_K_k       = 40            ;// [ohm-1m-2]
    const double g_KCC1_k    = 1e-2          ;// [ohm-1m-2]
    const double g_NBC_k     = 7.57e-1       ;// [ohm-1m-2]
    const double g_Cl_k      = 8.797e-1      ;// [ohm-1m-2]
    const double g_NKCC1_k   = 5.54e-2       ;// [ohm-1m-2]
    const double g_Na_k      = 1.314         ;// [ohm-1m-2]
    const double J_NaK_max   = 1.42e-3       ;// [uMm s-1]
    const double K_Na_k      = 10e3          ;// [uM]
    const double K_K_s       = 1.5e3         ;// [uM]
    const double k_C         = 7.35e-5       ;// [muM s-1]

// Perivascular Space constants:
	const double R_decay     = 0.05; 	// s^-1
	const double K_p_min 	 = 3e3; 	// uM

// BK channel constants:
    const double A_ef_k      = 3.7e-9        ; 				// m2       Area of an endfoot of an astrocyte, equal to Area astrocyte at synaptic cleft
    const double v_6         = 22e-3         ; 				// V        Constant
    const double v_4         = 14.5e-3       ;				// V        A measure of the spread of the distrubution
    const double psi_w       = 2.664         ; 				// s-1      A characteristic time
    const double G_BK_k      = 4.3e3         ; 				// pS      Constant estimation based on Ermentrout
    const double g_BK_k      = G_BK_k * 1e-12 / A_ef_k ;	// ohm-1m-2  Specific capacitance of the BK-Channel in units of Ostby
    const double VR_pa       = 0.001       	 ; 				// [-]       The estimated volume ratio of perivascular space to astrocyte: Model estimation
    const double VR_ps       = 0.001         ; 				// [-]       The estimated volume ratio of perivascular space to SMC: Model Estimation

// SMC constants:
    const double F_il 		= 7.5e2;		//[-] scaling factor to fit the experimental data of Filosa
    const double z_1 		=4.5;			//[-] parameter fitted on experimental data of Filosa
    const double z_2 		=-1.12e2;		//[-] parameter fitted on experimental data of Filosa
    const double z_3 		=4.2e-1;		//[-] parameter fitted on experimental data of Filosa
    const double z_4 		=-1.26e1;		//[-] parameter fitted on experimental data of Filosa
    const double z_5 		=-7.4e-2; 		//[-] parameter fitted on experimental data of Filosa
    const double Fmax_i		= 0.23;			// (microM/s)
    const double Kr_i 		= 1; 			// (microM) Half saturation constant for agonist-dependent Ca entry
    const double G_Ca		= 0.00129;		// (microM/mV/s)
    const double v_Ca1		= 100;			// (mV)
    const double v_Ca2		= -24;			// (mV)
    const double R_Ca		= 8.5;			// (mV)
    const double G_NaCa		= 0.00316;		// (microM/mV/s)
    const double c_NaCa		= 0.5;			// (microM)
    const double v_NaCa		= -30;
    const double B_i		= 2.025;
    const double cb_i		= 1;
    const double C_i		= 55;
    const double sc_i		= 2;
    const double cc_i		= 0.9;
    const double D_i		= 0.24;
    const double vd_i		= -100;
    const double Rd_i		= 250;
    const double L_i		= 0.025;
    const double gam		= 1970; 		// mVmicroM-1 The change in membrane potential by a scaling factor
    const double F_NaK		= 0.0432;
    const double G_Cl		= 0.00134;
    const double v_Cl		= -25;
    const double G_K		= 0.00446;
    const double vK_i		= -94;
    const double lam 		= 45;
    const double c_w		= 0;
    const double bet		= 0.13;
    const double v_Ca3		= -27; 			// correct
    const double R_K		= 12;
    const double k_i		= 0.1;

// Stretch-activated channels
    const double G_stretch   = 0.0061;       // uM mV-1 s-1   (stretch activated channels)
//    const double P_str       = 30;  		// Pressure will now come from vasculature
    const double Esac        = -18;          // mV
    const double alpha1      = 0.0074;
    const double sig0        = 500;

// EC constants:
    const double Fmax_j		= 0.23;			// [microM/s]
    const double Kr_j		= 1;
    const double B_j 		= 0.5;
    const double cb_j		= 1;
    const double C_j		= 5;
    const double sc_j		= 2;
    const double cc_j		= 0.9;
    const double D_j		= 0.24;
    const double L_j		= 0.025;
    const double G_cat 		= 0.66e-3; 		//!
    const double E_Ca		= 50;
    const double m3cat		= -0.18; 		//-6.18 changed value!
    const double m4cat 		= 0.37;
    const double JO_j 		= 0.029; 		//constant Ca influx (EC)
    const double C_m 		= 25.8;
    const double G_tot		= 6927;
    const double vK_j 		= -80;
    const double a1			= 53.3;
    const double a2			= 53.3;
    const double b			= -80.8;
    const double c 			= -0.4; 		//-6.4 changed value!
    const double m3b		= 1.32e-3;
    const double m4b		= 0.3;
    const double m3s		= -0.28;
    const double m4s		= 0.389;
    const double G_R		= 955;
    const double v_rest		= -31.1;
    const double k_j		= 0.1;

    const double J_PLC 		= 0.18;	//0.18 *****************************

    const double g_hat      = 0.5;
    const double p_hat      = 0.05;
    const double p_hatIP3   = 0.05;
    const double C_Hillmann = 1;
    const double K2_c        = 0.5 * C_Hillmann;
    const double K3_c        = 0.4 * C_Hillmann;
    const double K4_c        = 0.1 * C_Hillmann;
    const double K5_c        = 0.5 * C_Hillmann;
    const double K7_c        = 0.1 * C_Hillmann;
    const double gam_cross   = 17 * C_Hillmann;
	
// Neuron constants
	const double SC_coup	 = 9.5;
	const double tns 		 = 3;
	const double g 			 = 20;
	const double Farad 		 = 96.485;
	const double E_Cl_sa	 = -70;
	const double E_Cl_d		 = -70;
	const double Ra			 = 1.83e5;
	const double dhod 		 = 4.5e-2;
	const double As 		 = 1.586e-5;
	const double Ad			 = 2.6732e-4;
	const double Vs			 = 2.16e-9;
	const double Vd			 = 5.614e-9;
	const double fe			 = 0.15;
	const double Cm			 = 7.5e-7;
	const double ph			 = 26.6995;
	const double Mu			 = 8e-4;
	const double B0 		 = 500;
	const double gNaP_GHk	 = 2e-6;
	const double gKDR_GHk	 = 10e-5;
	const double gKA_GHk	 = 1e-5;
	const double gNMDA_GHk	 = 1e-5;
	const double gNaT_GHk    = 10e-5;
	const double gNaleak_sa	 = 4.1333e-5;
	const double gKleak_sa	 = 1.4623e-4;
	const double gClleak_sa	 = 41.333e-5;
	const double gNaleak_d	 = 4.1915e-5;
	const double gKleak_d	 = 1.4621e-4;
	const double gClleak_d	 = 41.915e-5;
	const double Imax		 = 0.052;
	const double O2_0		 = 2e-2;
	const double alph 		 = 0.05;
	const double D_Na 		 = 1.33e-5;
	const double D_K 		 = 1.96e-5;
	const double D_Cl 		 = 2.03e-5;
	const double K_init_e 	 = 2.9;
	const double Na_init_sa  = 10;
	const double Na_init_d 	 = 10;
	const double LU_R_init	 = 1.9341e-5;
	const double CBF_init 	 = 3.2e-2;
	const double O2_b 		 = 4e-2;
	const double gamm		 = 0.1;
	const double Mg			 = 1.2;
	

    double pt, e, r0, q, g, P_str; // pressure stuff

    // Initialise state variables
    double state_r, AMp_AM;
    double state_R_k,   state_N_Na_k, state_N_K_k, state_N_HCO3_k, state_N_Cl_k, state_N_Na_s, state_N_K_s, state_N_HCO3_s, state_K_p, state_w_k; // AC state
    double state_ca_i, state_ca_sr_i, state_v_i, state_w_i, state_ip3_i, state_K_i; // SMC state
    double state_ca_j, state_ca_er_j, state_v_j, state_ip3_j; // EC state
    double state_Mp, state_AM, state_AMp; // Mech state
    double state_v_sa, state_v_d, state_K_sa, state_Na_sa, state_Cl_sa, state_K_d, state_Na_d, state_Cl_d, state_K_e, state_Na_e, state_Cl_e, state_Buff_e, state_Buff_s, state_O2, state_CBV, state_K_buff, state_m1, state_m2, state_m3, state_m4, state_m5, state_m6, state_m7, state_m8, state_h1, state_h2, state_h3, state_h4,, state_h5, state_h6;

    // Fluxes
    double R_s, flu_N_Cl_s, flu_Na_k, flu_K_k, flu_HCO3_k, flu_Cl_k, flu_Na_s, flu_K_s, flu_HCO3_s, flu_Cl_s, flu_E_Na_k, flu_E_K_k, flu_E_Cl_k, flu_E_NBC_k, flu_E_BK_k, flu_J_NaK_k, flu_v_k, flu_J_KCC1_k, flu_J_NBC_k, flu_J_NKCC1_k, flu_J_Na_k, flu_J_K_k, flu_J_BK_k, flu_w_inf, flu_phi_w; // AC fluxes
    double flu_M, flu_E_K_i, flu_h_r, flu_v_cpl_i, flu_c_cpl_i, flu_I_cpl_i, flu_rho_i, flu_ip3_i, flu_SRuptake_i, flu_CICR_i, flu_extrusion_i, flu_leak_i, flu_VOCC_i, flu_NaCa_i, flu_NaK_i, flu_Cl_i, flu_K_i, flu_Kactivation_i, flu_degrad_i, flu_v_KIR_i, flu_G_KIR_i, flu_J_KIR_i, flu_J_stretch_i; // SMC fluxes
    double flu_v_cpl_j, flu_c_cpl_j, flu_I_cpl_j, flu_rho_j, flu_O_j, flu_ip3_j, flu_ERuptake_j, flu_CICR_j, flu_extrusion_j, flu_leak_j, flu_cation_j, flu_BKCa_j, flu_SKCa_j, flu_K_j, flu_R_j, flu_degrad_j, flu_J_stretch_j; // EC fluxes
    double flu_K1_c, flu_K6_c; // Mech fluxes   
	double J_KDR_sa, J_KA_sa, J_Kleak_sa, J_Kpump_sa, J_NaP_sa, J_Naleak_sa, J_Napump_sa, J_NaT_sa, J_KDR_d, J_KA_d, J_Kleak_d, J_Kpump_d, J_NMDA_K_d, J_NaP_d, J_Naleak_d, J_Napump_d, J_NMDA_Na_d, J_pump1_sa, J_pump1_d, J_pump1init_sa, J_pump1init_d;
	double 	flu_J_tot_sa, flu_J_tot_d, flu_J_K_tot_sa, flu_J_Na_tot_sa, flu_J_Cl_tot_sa, flu_J_K_tot_d, flu_J_Na_tot_d, flu_J_Cl_tot_d, flu_J_O2_vascular, flu_J_O2_background, flu_J_O2_pump, flu_J_CBF_norm;

    // State Variables:
    state_r  	  = u[i_radius];

    state_ca_i    = u[ca_i];
    state_ca_sr_i = u[ca_sr_i];
    state_v_i     = u[v_i];
    state_w_i     = u[w_i];
    state_ip3_i   = u[ip3_i];
    state_K_i     = u[K_i];

    state_R_k     = u[R_k];
    state_N_Na_k  = u[N_Na_k];
    state_N_K_k   = u[N_K_k];
    state_N_HCO3_k= u[N_HCO3_k];
    state_N_Cl_k  = u[N_Cl_k];
    state_N_Na_s  = u[N_Na_s];
    state_N_K_s   = u[N_K_s];
    state_N_HCO3_s= u[N_HCO3_s];
    state_K_p     = u[K_p];
    state_w_k     = u[w_k];

    state_ca_j    = u[ca_j];
    state_ca_er_j = u[ca_er_j];
    state_v_j     = u[v_j];
    state_ip3_j   = u[ip3_j];

    state_Mp      = u[Mp];
    state_AMp     = u[AMp];
    state_AM      = u[AM];

    state_v_sa	  = u[v_sa];
	state_v_d	  = u[v_d];
	state_K_sa	  = u[K_sa];
	state_Na_sa	  = u[Na_sa];
	state_Cl_sa	  = u[Cl_sa];
	state_K_d	  = u[K_d];
	state_Na_d	  = u[Na_d];
	state_Cl_d	  = u[Cl_d];
	state_K_e	  = u[K_e];
	state_Na_e	  = u[Na_e];
	state_Cl_e	  = u[Cl_e];
	
	state_Buff_e  = u[Buff_e];
	state_Buff_s  = u[Buff_s];
	state_O2	  = u[O2];
	state_CBV	  = u[CBV];
	state_K_buff  = u[K_buff];
	
	state_m1	  = u[m1];
	state_m2	  = u[m2];
	state_m3	  = u[m3];
	state_m4	  = u[m4];
	state_m5	  = u[m5];
	state_m6	  = u[m6];
	state_m7	  = u[m7];
	state_m8	  = u[m8];
	state_h1	  = u[h1];
	state_h2	  = u[h2];
	state_h3	  = u[h3];
	state_h4	  = u[h4];
	state_h5	  = u[h5];
	state_h6	  = u[h6];
	

// Fluxes:
    AMp_AM = state_AMp + state_AM;
    g  = pow(state_r, 4) / nvu_w->l;

    // pressure
    pt = 0.5 * (p + nvu_w->pcap);
    e  = 1.0 + nvu_w->a5 * AMp_AM;
    r0 = nvu_w->a3 * (1.0 - nvu_w->a4 * AMp_AM);
    //q  = (p - nvu_w->pcap) * g;

    // AC fluxes
    R_s        	    	= R_tot - state_R_k;                            // state_R_k is AC volume-area ratio, R_s is SC
    flu_N_Cl_s         	= state_N_Na_s + state_N_K_s - state_N_HCO3_s;  //
    flu_Na_k           	= state_N_Na_k / state_R_k;                     //
    flu_K_k            	= state_N_K_k / state_R_k;                      //
    flu_HCO3_k         	= state_N_HCO3_k / state_R_k;                   //
    flu_Cl_k           	= state_N_Cl_k / state_R_k;                     //
    flu_Na_s           	= state_N_Na_s / R_s;                       //
    flu_K_s            	= state_N_K_s / R_s;                        //
    flu_HCO3_s         	= state_N_HCO3_s / R_s;                     //
    flu_Cl_s           	= flu_N_Cl_s / R_s;                         //

    flu_E_Na_k         	= (R_gas * Temp) / (z_Na * F) * log(flu_Na_s / flu_Na_k);    // V
    flu_E_K_k          	= (R_gas * Temp) / (z_K  * F) * log(flu_K_s / flu_K_k );     // V
    flu_E_Cl_k         	= (R_gas * Temp) / (z_Cl * F) * log(flu_Cl_s / flu_Cl_k);    // V
    flu_E_NBC_k        	= (R_gas * Temp) / (z_NBC* F) * log((flu_Na_s * pow(flu_HCO3_s,2))/(flu_Na_k * pow(flu_HCO3_k,2)));     // V
    flu_E_BK_k         	= (R_gas * Temp) / (z_K  * F) * log(state_K_p / flu_K_k);   // V
    flu_J_NaK_k        	= J_NaK_max * ( pow(flu_Na_k,1.5) / ( pow(flu_Na_k,1.5) + pow(K_Na_k,1.5) ) ) * ( flu_K_s / (flu_K_s + K_K_s) );    // uMm s-1
    flu_v_k            	= ( g_Na_k * flu_E_Na_k + g_K_k * flu_E_K_k + g_Cl_k  * flu_E_Cl_k + g_NBC_k * flu_E_NBC_k - flu_J_NaK_k * F/unitcon + g_BK_k * state_w_k * flu_E_BK_k) / (g_Na_k + g_K_k + g_Cl_k + g_NBC_k + g_BK_k * state_w_k);  // V
    flu_J_KCC1_k       	= 1 * (R_gas * Temp * g_KCC1_k) / (pow(F,2)) * log(((flu_K_s) * (flu_Cl_s))/((flu_K_k)*(flu_Cl_k))) * unitcon;   //uMm s-1
    flu_J_NBC_k        	= g_NBC_k / F * ((flu_v_k) - (flu_E_NBC_k))*unitcon;       //uMm s-1
    flu_J_NKCC1_k     	= 1 * (g_NKCC1_k * R_gas * Temp) / (pow(F,2))  * log(((flu_K_s) * (flu_Na_s) * pow(flu_Cl_s,2)) /((flu_K_k) * (flu_Na_k) * pow(flu_Cl_k,2)))*unitcon;        //uMm s-1
    flu_J_Na_k   		= g_Na_k / F * (flu_v_k - flu_E_Na_k) * unitcon;              //uMm s-1
    flu_J_K_k    		= g_K_k  / F * ((flu_v_k) - (flu_E_K_k )) * unitcon;          //uMm s-1
    flu_J_BK_k   		= g_BK_k / F * state_w_k * (flu_v_k - flu_E_BK_k) * unitcon;  //uMm s-1
    flu_w_inf    		= 0.5 * (1+tanh(((flu_v_k)+v_6)/v_4));                            //[-]
    flu_phi_w    		= psi_w * cosh(((flu_v_k)+v_6)/(2*v_4));                          //s-1

    // SMC fluxes
    flu_M               = 1 - state_Mp - state_AM - state_AMp;
    flu_E_K_i           = ( R_gas * Temp ) / ( z_K  * F ) * unitcon * log( state_K_p / state_K_i );
    flu_h_r             = 0.1 * state_r; 												//(non-dimensional!)
    flu_v_cpl_i		    = - g_hat * ( state_v_i - state_v_j );
    flu_c_cpl_i         = - p_hat * ( state_ca_i - state_ca_j );
    flu_I_cpl_i         = - p_hatIP3 * ( state_ip3_i - state_ip3_j );
    flu_rho_i		    = 1;
    flu_ip3_i		    = Fmax_i *  pow(state_ip3_i,2) / ( pow(Kr_i,2) + pow(state_ip3_i,2) );
    flu_SRuptake_i      = B_i * pow(state_ca_i,2) / ( pow(state_ca_i,2) + pow(cb_i,2) );
    flu_CICR_i		    = C_i * pow(state_ca_sr_i,2) / ( pow(sc_i,2) + pow(state_ca_sr_i,2) ) *  ( pow(state_ca_i,4) ) / ( pow(cc_i,4) + pow(state_ca_i,4) );
    flu_extrusion_i	    = D_i * state_ca_i * (1 + ( state_v_i - vd_i ) / Rd_i );
    flu_leak_i 		    = L_i * state_ca_sr_i;
    flu_VOCC_i		    = G_Ca * ( state_v_i - v_Ca1 ) / ( 1 + exp( - ( state_v_i - v_Ca2 ) / ( R_Ca ) ) );
    flu_NaCa_i		    = G_NaCa * state_ca_i * ( state_v_i - v_NaCa ) / ( state_ca_i + c_NaCa ) ;
    flu_NaK_i		    = F_NaK;
    flu_Cl_i		    = G_Cl * (state_v_i - v_Cl);
    flu_K_i			    = G_K * state_w_i * ( state_v_i - vK_i );
    flu_Kactivation_i   = pow((state_ca_i + c_w),2) / ( pow((state_ca_i + c_w),2) + bet*exp(-(state_v_i - v_Ca3)/R_K) );  // see NO pathway!
    flu_degrad_i	    = k_i * state_ip3_i;
	P_str 				= (p*P0 + PCAP) / 2.0 * PA2MMHG;
    flu_J_stretch_i     = G_stretch/(1+exp(-alpha1*(P_str*state_r / flu_h_r - sig0))) * (state_v_i - Esac); 

    flu_v_KIR_i    		= z_1 * state_K_p / unitcon + z_2;                                  // mV           state_K_p,
    flu_G_KIR_i    		= exp( z_5 * state_v_i + z_3 * state_K_p / unitcon + z_4 );        // pS pF-1 =s-1  state_v_i, state_K_p
    flu_J_KIR_i    		= F_il/gam * (flu_G_KIR_i) * (state_v_i-(flu_v_KIR_i));            // mV s-1 //     state_v_i, state_K_p
    
    // EC fluxes
    flu_v_cpl_j			= - g_hat * ( state_v_j - state_v_i );
    flu_c_cpl_j			= - p_hat * ( state_ca_j - state_ca_i );
    flu_I_cpl_j			= - p_hatIP3 * ( state_ip3_j - state_ip3_i );
    flu_rho_j 			= 1;
    flu_O_j 			= JO_j;
    flu_ip3_j			= Fmax_j * ( pow(state_ip3_j,2) ) / ( pow(Kr_j,2) + pow(state_ip3_j,2) );
    flu_ERuptake_j      = B_j * ( pow(state_ca_j,2) ) / ( pow(state_ca_j,2) + pow(cb_j,2) );
    flu_CICR_j			= C_j *  ( pow(state_ca_er_j,2) ) / ( pow(sc_j,2) + pow(state_ca_er_j,2) ) *  ( pow(state_ca_j,4) ) / ( pow(cc_j,4) + pow(state_ca_j,4) );
    flu_extrusion_j     = D_j * state_ca_j;
    flu_leak_j          = L_j * state_ca_er_j;
    flu_cation_j 		= G_cat * ( E_Ca - state_v_j) * 0.5 * ( 1 + tanh( ( log10( state_ca_j ) - m3cat )  /  m4cat  ) );
    flu_BKCa_j 			= 0.2 * ( 1 + tanh( ( (  log10(state_ca_j) - c) * ( state_v_j - b ) - a1 ) / ( m3b* pow(( state_v_j + a2 * ( log10( state_ca_j ) - c ) - b),2) + m4b ) ) );
    flu_SKCa_j 			= 0.3 * ( 1 + tanh( ( log10(state_ca_j) - m3s ) /  m4s ));
    flu_K_j 			= G_tot * ( state_v_j - vK_j ) * ( flu_BKCa_j + flu_SKCa_j );
    flu_R_j 			= G_R * ( state_v_j - v_rest);
    flu_degrad_j 		= k_j * state_ip3_j;
    flu_J_stretch_j     = G_stretch / (1 + exp(-alpha1*(P_str * state_r / flu_h_r - sig0))) * (state_v_j - Esac);
		
	// Neuron fluxes
	flu_J_tot_sa		= J_Na_tot_sa + J_K_tot_sa + J_Cl_tot_sa;
	flu_J_tot_d			= J_Na_tot_d + J_K_tot_d + J_Cl_tot_d;
	flu_J_K_tot_sa		= J_KDR_sa + J_KA_sa + J_Kleak_sa + J_Kpump_sa;
	flu_J_Na_tot_sa		= J_NaP_sa + J_Naleak_sa + J_Napump_sa + J_NaT_sa;
	flu_J_Cl_tot_sa		= p.gClleak_sa * (v_sa - p.E_Cl_sa);
	flu_J_K_tot_d		= J_KDR_d + J_KA_d + J_Kleak_d + J_Kpump_d + J_NMDA_K_d;
	flu_J_Na_tot_d		= J_NaP_d + J_Naleak_d + J_Napump_d + J_NMDA_Na_d;
	flu_J_Cl_tot_d		= p.gClleak_d * (v_d - p.E_Cl_d);
	flu_J_O2_vascular	= CBF .* ((p.O2_b - O2) ./ (p.O2_b - p.O2_0));
	flu_J_O2_background = p.CBF_init * P_02 * (1 - p.gamm);
	flu_J_O2_pump 		= p.CBF_init * P_02 * p.gamm .* ((J_pump1_sa + J_pump1_d) ./ (J_pump1init_sa + J_pump1init_d));
	flu_J_CBF_norm 		= CBF/0.03219;


// Mech fluxes
    flu_K1_c         	= gam_cross * pow(state_ca_i,3);
    flu_K6_c        	= flu_K1_c;

// Differential Equations:
    du[i_radius	] = -nvu_w->a1 * e * (state_r / r0 - 1.0) + nvu_w->a2 * state_r * pt; // Radius (non-dimensional!)

    //SC:
    du[ N_Na_s  ] = - k_C * K_input(t,x,y) - du[ N_Na_k];                           // uMm s-1
    du[ N_K_s   ] = k_C * K_input(t,x,y) - du[ N_K_k] - flu_J_BK_k + R_s * ( (state_K_e - flu_K_s) / tau2);                 // uMm s-1
    du[ N_HCO3_s] = - du[ N_HCO3_k];                                                // uMm s-1

    //AC:
    du[ R_k     ] = L_p * (flu_Na_k + flu_K_k + flu_Cl_k + flu_HCO3_k - flu_Na_s - flu_K_s - flu_Cl_s - flu_HCO3_s + X_k / state_R_k);  // m s-1
    du[ N_Na_k  ] = -flu_J_Na_k - 3 * flu_J_NaK_k + flu_J_NKCC1_k + flu_J_NBC_k;    // uMm s-1
    du[ N_K_k   ] = -flu_J_K_k + 2 * flu_J_NaK_k + flu_J_NKCC1_k + flu_J_KCC1_k -flu_J_BK_k; // uMm s-1
    du[ N_HCO3_k] = 2 * flu_J_NBC_k;                                                // uMm s-1
    du[ N_Cl_k  ] = du[ N_Na_k] + du[ N_K_k] - du[ N_HCO3_k];                       // uMm s-1, modified equation compared to the one of Ostby  //
    du[ w_k     ] = flu_phi_w * (flu_w_inf - state_w_k);                            // s-1

    //PVS:
    du[ K_p     ] = flu_J_BK_k / (VR_pa * state_R_k) + flu_J_KIR_i / VR_ps - R_decay * (state_K_p - K_p_min); // + ( (state_K_e - state_K_p) / tau);         // uM s-1

    //SMC:
    du[ ca_i    ] = flu_c_cpl_i + flu_rho_i * ( flu_ip3_i - flu_SRuptake_i + flu_CICR_i - flu_extrusion_i + flu_leak_i - flu_VOCC_i + flu_NaCa_i + 0.1* flu_J_stretch_i);
    du[ ca_sr_i ] = flu_SRuptake_i - flu_CICR_i - flu_leak_i ;
    du[ v_i     ] = flu_v_cpl_i + gam * ( - flu_NaK_i - flu_Cl_i - 2 * flu_VOCC_i - flu_NaCa_i - flu_K_i - flu_J_stretch_i - flu_J_KIR_i );
    du[ w_i     ] = lam * (flu_Kactivation_i - state_w_i ) ;
    du[ ip3_i   ] = flu_I_cpl_i - flu_degrad_i ;          // **
    du[ K_i     ] = - flu_J_KIR_i - flu_K_i + flu_NaK_i;                                            // uM s-1

    //EC:
    du[ca_j     ] = flu_c_cpl_j + flu_rho_j * ( flu_ip3_j - flu_ERuptake_j + flu_CICR_j - flu_extrusion_j + flu_leak_j + flu_cation_j + flu_O_j + flu_J_stretch_j ) ;
    du[ca_er_j  ] = flu_ERuptake_j - flu_CICR_j - flu_leak_j ;
    du[v_j      ] = flu_v_cpl_j - 1/C_m * ( flu_K_j + flu_R_j ) ;
    du[ip3_j    ] = flu_I_cpl_j + J_PLC - flu_degrad_j ;  // **

    // Mech:
    du[ Mp   	] = K4_c * state_AMp + flu_K1_c * flu_M - (K2_c + K3_c) * state_Mp;
    du[ AMp  	] = K3_c * state_Mp + flu_K6_c * state_AM - (K4_c + K5_c) * state_AMp;
    du[ AM   	] = K5_c * state_AMp - ( K7_c + flu_K6_c ) * state_AM;
	

    // Neuron - ions
	du[ v_sa	] = 1/p.Cm * ( -flu_J_tot_sa + 1 / (2 * p.Ra * p.dhod.^2) * (state_v_d - state_v_sa) + current_input(t,x,y) );
	du[ v_d		] = 1/p.Cm * (-flu_J_tot_d + 1 / (2 * p.Ra * p.dhod.^2) * (state_v_sa - state_v_d));
    du[ K_sa	] = -p.As / (p.Farad * p.Vs) * flu_J_K_tot_sa + p.D_K * (p.Vd + p.Vs) ./ (2 * p.dhod.^2 * p.Vs) * (state_K_d - state_K_sa);
    du[ Na_sa	] = -p.As / (p.Farad * p.Vs) * flu_J_Na_tot_sa + p.D_Na * (p.Vd + p.Vs) ./ (2 * p.dhod.^2 * p.Vs) * (state_Na_d - state_Na_sa);
	du[ Cl_sa	] = -p.As / (p.Farad * p.Vs) * flu_J_Cl_tot_sa + p.D_Cl * (p.Vd + p.Vs) ./ (2 * p.dhod.^2 * p.Vs) * (state_Cl_d - state_Cl_sa);
	du[ K_d		] = -p.Ad / (p.Farad * p.Vd) * flu_J_K_tot_d + p.D_K * (p.Vs + p.Vd) ./ (2 * p.dhod.^2 * p.Vd) * (state_K_sa - state_K_d);
	du[ Na_d	] = -p.Ad / (p.Farad * p.Vd) * flu_J_Na_tot_d + p.D_Na * (p.Vs + p.Vd) ./ (2 * p.dhod.^2 * p.Vd) * (state_Na_sa - state_Na_d);
	du[ Cl_d	] = -p.Ad / (p.Farad * p.Vd) * flu_J_Cl_tot_d + p.D_Cl * (p.Vs + p.Vd) ./ (2 * p.dhod.^2 * p.Vd) * (state_Cl_sa - state_Cl_d);
	
	du[ Na_e	] = 1/(p.Farad * p.fe) * (((p.As * flu_J_Na_tot_sa) / p.Vs) + ((p.Ad * flu_J_Na_tot_d) / p.Vd));
	du[ Cl_e	] = 1/(p.Farad * p.fe) * (((p.As * flu_J_Cl_tot_sa) / p.Vs) + ((p.Ad * flu_J_Cl_tot_d) / p.Vd));
	
	// Neuron - other
	du[ Buff_e	] = p.Mu * state_K_e .* (p.B0 - state_Buff_e) ./ (1 + exp(-((state_K_e - 5.5) ./ 1.09))) - (p.Mu * state_Buff_e);
	du[ Buff_s	] = p.Mu * state_K_buff .* (p.B0 - state_Buff_s) ./ (1 + exp(-((state_K_buff - 5.5) ./ 1.09))) - (p.Mu * state_Buff_s);
	du[ O2		] = flu_J_O2_vascular - flu_J_O2_background - flu_J_O2_pump;
	du[ CBV		] = 1/(p.tns + p.g) .* ( flu_J_CBF_norm - state_CBV.^2.5 ); 
	du[ K_buff	] = p.SC_coup ./ (p.Farad * p.fe) * ( (p.As .* flu_J_K_tot_sa) / p.Vs  + ( p.Ad .* flu_J_K_tot_d) / (p.Vd ) ) - du[ Buff_s];
	
	du[ K_e		] = 1/(p.Farad * p.fe) * (((p.As * flu_J_K_tot_sa) / p.Vs)  + ((p.Ad * flu_J_K_tot_d) / p.Vd)) - du[ Buff_e] + ECS_input(t,x,y);
	
	// Neuron Gating Variables
	du[ m1		] =	1000 * ((m1alpha .* (1 - state_m1)) - (m1beta .* state_m1));
	du[ m2		] =	1000 * ((m2alpha .* (1 - state_m2)) - (m2beta .* state_m2));
	du[ m3		] =	1000 * ((m3alpha .* (1 - state_m3)) - (m3beta .* state_m3));
	du[ m4		] = 1000 * ((m4alpha .* (1 - state_m4)) - (m4beta .* state_m4));	
	du[ m5		] =	1000 * ((m5alpha .* (1 - state_m5)) - (m5beta .* state_m5));
	du[ m6		] =	1000 * ((m6alpha .* (1 - state_m6)) - (m6beta .* state_m6));
	du[ m7		] =	1000 * ((m7alpha .* (1 - state_m7)) - (m7beta .* state_m7)); 
	du[ m8		] =	1000 * ((m8alpha .* (1 - state_m8)) - (m8beta .* state_m8));
	du[ h1		] = 1000 * ((h1alpha .* (1 - state_h1)) - (h1beta .* state_h1));	
	du[ h2		] =	1000 * ((h2alpha .* (1 - state_h2)) - (h2beta .* state_h2));
	du[ h3		] =	1000 * ((h3alpha .* (1 - state_h3)) - (h3beta .* state_h3));
	du[ h4		] =	1000 * ((h4alpha .* (1 - state_h4)) - (h4beta .* state_h4));
	du[ h5		] =	1000 * ((h5alpha .* (1 - state_h5)) - (h5beta .* state_h5));
	du[ h6		] =	1000 * ((h6alpha .* (1 - state_h6)) - (h6beta .* state_h6));

    // State variables so they can be plotted in Paraview, but only for one time (initial condition set in nvu_ics, use t for when the signal is turned on)
    du[PLC_i		] = 0;
    du[current_i	] = 0;
}

// Time-varying pressure at the root of the tree. 1 is nominal value. If
// you want to work in unscaled units, make sure you *multiply* by P0
// afterwards
double nvu_p0(double t)
{
    //double p0 = 1. * 8000 / P0; // 8000 Pa   original: 1.5 * 8000 / P0;
    //double p0 = (0.5 * sin(t) + 1) * 8000 / P0; //
    double p0 = 1.5 * 8000 / P0;	// no time dependence?

    return p0;
}

// Space- & time-varying K+ input signal (simulating neuronal activity)
double current_input(double t, double x, double y)
{
    double current_input_min 	= 0;
    double current_input_max 	= 0.02;
    double t_up   		= 1000;
    double t_down 		= 2000;
    double lengthpulse 	= t_down - t_up;
    double lengtht1 	= 10;
    double t0 			= t_up;
    double t1 			= t0 + lengtht1;
    double t2 			= t0 + lengthpulse;
    double t3 			= t1 + lengthpulse;
    int alpha 			= 2;
    int beta 			= 5;
    double ampl = 3;
    double ramp = 0.003;
    double x_centre = 0;//-0.0008;
    double y_centre = 0;//-0.0008;
    double deltat		= 10;
    double gab 			= factorial(alpha + beta - 1);
    double ga 			= factorial(alpha - 1);
    double gb 			= factorial(beta - 1);

    //double current_space = fmin(1.0,ampl*(exp(- ((pow((x-x_centre),2)+pow((y-y_centre),2)) / (2 * pow(ramp,2))))));
    //double current_space =((0.5 + 0.5 *(tanh(1e5 * (x-0.0004)+1))) *(0.5 + 0.5 *(tanh(1e5 *(y-0.0004)+1))));

    double current_space;
    // only in corner
    if (x <= 0 && y <= 0)
    {
        current_space = 1;
    }
    else
    {
        current_space = 0;
    }

    double current_time;
    if (t >= t0 && t <= t1)
    {
        //current_time = F_input * gab / (ga * gb) * pow((1-(t-t0) / deltat),(beta - 1)) * pow(((t - t0) / deltat),(alpha-1)); 
        current_time = 1 * gab / (ga * gb) * pow((1 - (t - t0) / deltat), (beta - 1)) * pow(((t - t0) / deltat), (alpha-1));		//TODO: rectpuls for input
    }
    else if (t >= t2 && t <= t3)
    {
        //current_time = - F_input;
        current_time = - 1;
    }   
    else
    {
        current_time = 0;
    }

    double current_out = current_input_min + (current_input_max - current_input_min) * current_space * current_time; // 0 if t3 < t or x,y <= 0
    return current_out;
}


// Space- & time-varying PLC input signal
double PLC_input(double t, double x, double y)
{
    double PLC_min = 0.18;
    double PLC_max = 0.4;
    double t_up   = 1000;
    double t_down = 9000;
    double ampl = 3;
    double ramp = 0.003;//0.002;
    double x_centre = 0; // 0.0008 -> n_bif = 7; python: ((((2**(n_bif-1))**0.5)/4)*0.0004)
    double y_centre = 0;
    double PLC_space = fmin(1.0, ampl * (exp(-((pow((x - x_centre), 2) + pow((y - y_centre), 2)) / (2 * pow(ramp, 2))))));
    double PLC_time = 0.5 * tanh((t - t_up) / 0.05) - 0.5 * tanh((t - t_down) / 0.05);

    double PLC_out = PLC_min + (PLC_max-PLC_min) * PLC_space * PLC_time;
    //double PLC_out = PLC_space; // no time-dependency
    return PLC_out;
}

// Space- & time-varying ECS K+ input signal
double ECS_input(double t, double x, double y)
{
    double ECS_max 		= 9e3;
    double t_up   		= 100;
    double t_down 		= 200;
    double lengthpulse 	= t_down - t_up;
    double lengtht1 	= 20;
    double t0 			= t_up;
    double t1 			= t0 + lengtht1;
    double t2 			= t0 + lengthpulse;
    double t3 			= t1 + lengthpulse;

    double ampl = 3;
    double ramp = 0.003;
    double x_centre = 0;
    double y_centre = 0;

    double ECS_space = fmin(1.0, ampl * (exp(-((pow((x - x_centre), 2) + pow((y - y_centre), 2)) / (2 * pow(ramp, 2))))));

    double ECS_time;
    if (t >= t0 && t <= t1)
    {
        ECS_time = 1;
    }
    else if (t >= t2 && t <= t3)
    {
    	ECS_time = - 1;
    }
    else
    {
    	ECS_time = 0;
    }

    double ECS_out = ECS_max * ECS_space * ECS_time;

    return ECS_out;
}

double factorial(int c)
{
    double result = 1;

    for (int n = 1; n < c; n++)
    {
        result = result * n;
    }

    return result;
}

// Initial conditions. If you want spatial inhomegeneity, make it a
// function of the coordinates x and y. u0 is already allocated, you just
// need to fill in the entries
void nvu_ics(double *u0, double x, double y, nvu_workspace *nvu_w)
{
	
	//TODO: get ICs from stable state in single NVU code

    u0[i_radius]  = 1;   						//0

    u0[R_k]       = 6.1e-8;                    //1
    u0[N_Na_k]    = 9.9796e-4;                 //2
    u0[N_K_k]     = 5.52782e-3;                //3
    u0[N_HCO3_k]  = 0.58804e-3;                //4
    u0[N_Cl_k]    = 0.32879e-3;                //5
    u0[w_k]       = 0.1815e-3;                 //10

    u0[N_Na_s]    = 4.301041e-3;               //6
    u0[N_K_s]     = 0.0807e-3;                 //7
    u0[N_HCO3_s]  = 0.432552e-3;               //8 

    u0[K_p]       = 3e3;                       //9

    u0[ca_i]      = 0.1;                       //11
    u0[ca_sr_i]   = 0.1;                       //12
    u0[v_i]       = -60;                       //13
    u0[w_i]       = 0.1;                       //14
    u0[ip3_i]     = 0.1;                       //15
    u0[K_i]       = 1e5;                       //16

    u0[ca_j]      = 0.1;                       //17
    u0[ca_er_j]   = 0.1;                       //18
    u0[v_j]       = -75;                       //19
    u0[ip3_j]     = 0.1;                       //20

    u0[Mp]        = 0.25;                      //21
    u0[AMp]       = 0.25;                      //22
    u0[AM]        = 0.25;                      //23

    u0[v_sa]	  = -70;
	u0[v_d]	      = -70;
	u0[K_sa]	  = 133.5;
	u0[Na_sa]	  = 9.9854;
	u0[Cl_sa]	  = 10.464;
	u0[K_d]	  	  = 133.5;
	u0[Na_d]	  = 9.9853;
	u0[Cl_d]	  = 10.464;
	u0[K_e]	  	  = 3.5;
	u0[Na_e]	  = 139.76;
	u0[Cl_e]	  = 144.1;
	
	u0[Buff_e]	  = 170;
	u0[Buff_s]	  = 170;
	u0[O2]		  = 0.0227;
	u0[CBV]	  	  = 1;
	u0[K_buff]	  = 3.5;
	
	u0[m1]	  = 0.0128;
	u0[m2]	  = 0.001;
	u0[m3]	  = 0.119;
	u0[m4]	  = 0.012;
	u0[m5]	  = 0.00087;
	u0[m6]	  = 0.0012;
	u0[m7]	  = 0.1193;
	u0[m8]	  = 0.005;
	u0[h1]	  = 0.9718;
	u0[h2]	  = 0.12;
	u0[h3]	  = 0.9718;
	u0[h4]	  = 0.99;
	u0[h5]	  = 0.12;
	u0[h6]	  = 0.996;

    // Only here so they can be shown in Paraview for some time when the signals are turned on (as a check), so choose t within t0 and t1
    u0[PLC_i]     	 = PLC_input(150,x,y);
    u0[current_i]    = current_input(150,x,y);

}
