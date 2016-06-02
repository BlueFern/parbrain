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
    int dfdp_pattern[] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 	// neq * 1 column vector

    // Column: variables the eq depends on
    // Row: variables that depend on that eq variable
    int dfdx_pattern[] = {1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0, 		//remember indexing starts at 0
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    		              0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
    		              0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    		              1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,0,0,0,1,1,1,0,0,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
    		              1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,
						  0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }; // column major, neq*neq

    // Initialise the workspace
    nvu_w = malloc(sizeof *nvu_w);
    nvu_w->neq = 28; //27;

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
	const double Farad       = 96500         ;// [C mol-1] Faradays constant.
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

    const double J_PLC 		= 0.4;	//0.18 *****************************

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

    // ECS:
    const double VR_pe       = 0.001; 	// [-]       The estimated volume ratio of PVS to ECS: Model Estimation
    // tau is dx^2 / 2D where dx is the length and D is diffusion rate
    const double tau 		 = 0.7; 	// (sec) characteristic time scale for ion to travel one cell length (SMC length is ~50 um)
    const double tau2 		 = 2.8; 	// (sec) characteristic time scale for ion to travel from PVS to SC (AC length is ~100 um, based on protoplasmic astrocyte process length of ~50 um)

    double pt, e, r0, q, g; // pressure stuff

    // Initialise state variables
    double state_r, AMp_AM;
    double state_R_k,   state_N_Na_k, state_N_K_k, state_N_HCO3_k, state_N_Cl_k, state_N_Na_s, state_N_K_s, state_N_HCO3_s, state_K_p, state_w_k; // AC state
    double state_ca_i, state_ca_sr_i, state_v_i, state_w_i, state_ip3_i, state_K_i; // SMC state
    double state_ca_j, state_ca_er_j, state_v_j, state_ip3_j; // EC state
    double state_Mp, state_AM, state_AMp; // Mech state
    double state_K_e; // ECS state
    double state_PLC_i, state_K_df_i, state_K_flux_i; // input

    // Fluxes
    double R_s, flu_N_Cl_s, flu_Na_k, flu_K_k, flu_HCO3_k, flu_Cl_k, flu_Na_s, flu_K_s, flu_HCO3_s, flu_Cl_s, flu_E_Na_k, flu_E_K_k, flu_E_Cl_k, flu_E_NBC_k, flu_E_BK_k, flu_J_NaK_k, flu_v_k, flu_J_KCC1_k, flu_J_NBC_k, flu_J_NKCC1_k, flu_J_Na_k, flu_J_K_k, flu_J_BK_k, flu_w_inf, flu_phi_w; // AC fluxes
    double flu_M, flu_E_K_i, flu_h_r, flu_v_cpl_i, flu_c_cpl_i, flu_I_cpl_i, flu_rho_i, flu_ip3_i, flu_SRuptake_i, flu_CICR_i, flu_extrusion_i, flu_leak_i, flu_VOCC_i, flu_NaCa_i, flu_NaK_i, flu_Cl_i, flu_K_i, flu_Kactivation_i, flu_degrad_i, flu_v_KIR_i, flu_G_KIR_i, flu_J_KIR_i, flu_J_stretch_i; // SMC fluxes
    double flu_v_cpl_j, flu_c_cpl_j, flu_I_cpl_j, flu_rho_j, flu_O_j, flu_ip3_j, flu_ERuptake_j, flu_CICR_j, flu_extrusion_j, flu_leak_j, flu_cation_j, flu_BKCa_j, flu_SKCa_j, flu_K_j, flu_R_j, flu_degrad_j, flu_J_stretch_j; // EC fluxes
    double flu_K1_c, flu_K6_c; // Mech fluxes   
    double P_str;

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

    state_K_e	  = u[K_e];

    state_PLC_i   = u[PLC_i];
    state_K_df_i  = u[K_df_i];
    state_K_flux_i= u[K_flux_i];

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

    flu_E_Na_k         	= (R_gas * Temp) / (z_Na * Farad) * log(flu_Na_s / flu_Na_k);    // V
    flu_E_K_k          	= (R_gas * Temp) / (z_K  * Farad) * log(flu_K_s / flu_K_k );     // V
    flu_E_Cl_k         	= (R_gas * Temp) / (z_Cl * Farad) * log(flu_Cl_s / flu_Cl_k);    // V
    flu_E_NBC_k        	= (R_gas * Temp) / (z_NBC* Farad) * log((flu_Na_s * pow(flu_HCO3_s,2))/(flu_Na_k * pow(flu_HCO3_k,2)));     // V
    flu_E_BK_k         	= (R_gas * Temp) / (z_K  * Farad) * log(state_K_p / flu_K_k);   // V
    flu_J_NaK_k        	= J_NaK_max * ( pow(flu_Na_k,1.5) / ( pow(flu_Na_k,1.5) + pow(K_Na_k,1.5) ) ) * ( flu_K_s / (flu_K_s + K_K_s) );    // uMm s-1
    flu_v_k            	= ( g_Na_k * flu_E_Na_k + g_K_k * flu_E_K_k + g_Cl_k  * flu_E_Cl_k + g_NBC_k * flu_E_NBC_k - flu_J_NaK_k * Farad/unitcon + g_BK_k * state_w_k * flu_E_BK_k) / (g_Na_k + g_K_k + g_Cl_k + g_NBC_k + g_BK_k * state_w_k);  // V
    flu_J_KCC1_k       	= flux_ft(t,x,y) * (R_gas * Temp * g_KCC1_k) / (pow(Farad,2)) * log(((flu_K_s) * (flu_Cl_s))/((flu_K_k)*(flu_Cl_k))) * unitcon;   //uMm s-1
    flu_J_NBC_k        	= g_NBC_k / Farad * ((flu_v_k) - (flu_E_NBC_k))*unitcon;       //uMm s-1
    flu_J_NKCC1_k     	= flux_ft(t,x,y) * (g_NKCC1_k * R_gas * Temp) / (pow(Farad,2))  * log(((flu_K_s) * (flu_Na_s) * pow(flu_Cl_s,2)) /((flu_K_k) * (flu_Na_k) * pow(flu_Cl_k,2)))*unitcon;        //uMm s-1
    flu_J_Na_k   		= g_Na_k / Farad * (flu_v_k - flu_E_Na_k) * unitcon;              //uMm s-1
    flu_J_K_k    		= g_K_k  / Farad * ((flu_v_k) - (flu_E_K_k )) * unitcon;          //uMm s-1
    flu_J_BK_k   		= g_BK_k / Farad * state_w_k * (flu_v_k - flu_E_BK_k) * unitcon;  //uMm s-1
    flu_w_inf    		= 0.5 * (1+tanh(((flu_v_k)+v_6)/v_4));                            //[-]
    flu_phi_w    		= psi_w * cosh(((flu_v_k)+v_6)/(2*v_4));                          //s-1

    // SMC fluxes
    flu_M               = 1 - state_Mp - state_AM - state_AMp;
    flu_E_K_i           = ( R_gas * Temp ) / ( z_K  * Farad ) * unitcon * log( state_K_p / state_K_i );
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
    du[ K_p     ] = flu_J_BK_k / (VR_pa * state_R_k) + flu_J_KIR_i / VR_ps - R_decay * (state_K_p - K_p_min) + ( (state_K_e - state_K_p) / tau);         // uM s-1

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

    //ECS:				smc efflux				SC flux					 				PVS flux								decay term
    du[ K_e		] = - flu_NaK_i + flu_K_i - ( (state_K_e - flu_K_s) / tau2) - VR_pe * ( (state_K_e - state_K_p) / tau); // - 0.05 * state_K_e;
//    du[ K_e		] = 0; // for only the diffusion eq

    // State variables so they can be plotted in Paraview, but only for one time (initial condition set in nvu_ics, use t for when the signal is turned on)
    du[PLC_i	] = 0;
    du[K_df_i	] = 0;
    du[K_flux_i	] = 0;
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
double K_input(double t, double x, double y)
{
    double K_input_min 	= 0;
    double K_input_max 	= 2.5;
    double t_up   		= 0;				// Signal starts at time = 5 ***
    double t_down 		= 350;
    double lengthpulse 	= t_down - t_up;
    double lengtht1 	= 15;
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

    //double K_space = fmin(1.0,ampl*(exp(- ((pow((x-x_centre),2)+pow((y-y_centre),2)) / (2 * pow(ramp,2))))));
    //double K_space =((0.5 + 0.5 *(tanh(1e5 * (x-0.0004)+1))) *(0.5 + 0.5 *(tanh(1e5 *(y-0.0004)+1))));

    double K_space;
    // only in corner
    if (x <= 0 && y <= 0)
    {
        K_space = 1;
    }
    else
    {
        K_space = 0;
    }

    double K_time;
    if (t >= t0 && t <= t1)
    {
        //K_time = F_input * gab / (ga * gb) * pow((1-(t-t0) / deltat),(beta - 1)) * pow(((t - t0) / deltat),(alpha-1)); 
        K_time = 1 * gab / (ga * gb) * pow((1 - (t - t0) / deltat), (beta - 1)) * pow(((t - t0) / deltat), (alpha-1));
    }
    else if (t >= t2 && t <= t3)
    {
        //K_time = - F_input;
        K_time = - 1;
    }   
    else
    {
        K_time = 0;
    }

    double K_out = K_input_min + (K_input_max - K_input_min) * K_space * K_time; // 0 if t3 < t or x,y <= 0
    return K_out;
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

// Block function to switch cotransporter channels on and off (K+ input)
double flux_ft(double t, double x, double y)
{
    double flux_min 	= 0;
    double flux_max 	= 1;
    double t_up   		= 0;					// Channels turn on at time = 5 ***
    double t_down 		= 350;
    double lengthpulse 	= t_down - t_up;
    double lengtht1 	= 15;
    double t0 			= t_up;
    double t1 			= t0 + lengtht1;
    double ampl = 3;
    double ramp = 0.003;
    double x_centre = 0;//-0.0008;
    double y_centre = 0;//-0.0008;

    double flux_time = 0.5 * tanh((t - t0) / 0.0005) - 0.5 * tanh((t - t1 - lengthpulse) / 0.0005);
    //double flux_time = 0.5 * tanh((t-t0)/0.005) - 0.5 * tanh((t-t1-lengthpulse)/0.005);

    //double flux_space = fmin(1.0,ampl*(exp(- ((pow((x-x_centre),2)+pow((y-y_centre),2)) / (2 * pow(ramp,2))))));
    //double flux_space =((0.5 + 0.5 *(tanh(1e5 * (x-0.0004)+1))) *(0.5 + 0.5 *(tanh(1e5 *(y-0.0004)+1))));  

    double flux_space;
    // only in corner
    if (x <= 0 && y <= 0)
    {
        flux_space = 1;
    }
    else
    {
        flux_space = 0;
    }

    double flux_out = flux_min + (flux_max-flux_min) * flux_time * flux_space;

    return flux_out;
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



// Initial conditions. If you want spatial inhomegeneity, make it a
// function of the coordinates x and y. u0 is already allocated, you just
// need to fill in the entries
void nvu_ics(double *u0, double x, double y, nvu_workspace *nvu_w)
{

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

    u0[K_e]		  = 3e3;

    // Only here so they can be shown in Paraview for some time when the signals are turned on (as a check), so choose t within t0 and t1
    u0[PLC_i]     = PLC_input(15,x,y);
    u0[K_df_i]    = K_input(15,x,y);
    u0[K_flux_i]  = flux_ft(15,x,y);

}
