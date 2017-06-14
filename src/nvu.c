// NVU 2.0

// TODO: ini file with the following: J_PLC, switches (TRPV, NO, K, Ca etc), ECS input size, K input size, stimulation time t_0 and t_2 for each input, type of spatial input (centre vs corner etc), end of simulation time, number of cores etc etc

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
static const double PA2MMHG    = 0.00750061683; // convert from Pa to mmHg

// nvu_init: this user-supplied function does any precomputation required
// for the model.
nvu_workspace *nvu_init(void)
{
    nvu_workspace *nvu_w;

    // Initialise the workspace
    nvu_w = malloc(sizeof *nvu_w);
    nvu_w->neq = 67;                 // Must specify here

    // Sparsity patterns (approximated by matrices full of 1s)
    int dfdp_pattern[nvu_w->neq],i;
    for(i = 0; i < nvu_w->neq ; dfdp_pattern[i++] = 1);

    int dfdx_pattern[nvu_w->neq*nvu_w->neq];
    for(i = 0; i < nvu_w->neq*nvu_w->neq ; dfdx_pattern[i++] = 1);

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

    // Allocate other nvu workspace parameters ** remove obsolete variables
//    double Rstar = R0;
//    double hstar = HRR * Rstar;
//    nvu_w->a1 = E0 * T0 * Rstar / (ETA * R0);
//    nvu_w->a2 = P0 * Rstar * T0 / (ETA * hstar);
//    nvu_w->a3 = Rstar / R0;
//    nvu_w->a4 = 1 - RSCALE;
//    nvu_w->a5 = EACTIVE / EPASSIVE - 1;
    nvu_w->pcap  = PCAP / P0; // pressure is nondimensional
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

/****** Model parameters ******/

// general constants:
	const double F           = 96500         ;// [C mol-1] Faradays constant.
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

// Perivascular Space constants:
	const double R_decay     = 0.15; 	// s^-1
	const double K_p_min 	 = 3e3; 	// uM

// BK channel constants:
    const double A_ef_k      = 3.7e-9        ; 				// m2       Area of an endfoot of an astrocyte, equal to Area astrocyte at synaptic cleft
    const double v_4         = 8e-3       ;				// V        A measure of the spread of the distribution
    const double psi_w       = 2.664         ; 				// s-1      A characteristic time
    const double G_BK_k      = 225         ; 				// !!!
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
    const double v_Ca3		= -27; 			// correct
    const double R_K		= 12;
    const double k_i		= 0.1;

// Stretch-activated channels
    const double G_stretch   = 0.0061;       // uM mV-1 s-1   (stretch activated channels)
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
    const double J_PLC 		= 0.11;	//0.11 or 0.3 for oscillations *****************************
    const double g_hat      = 0.5;
    const double p_hat      = 0.05;
    const double p_hatIP3   = 0.05;
    const double C_Hillmann = 1;
    const double K3_c        = 0.4 * C_Hillmann;
    const double K4_c        = 0.1 * C_Hillmann;
    const double K7_c        = 0.1 * C_Hillmann;
    const double gam_cross   = 17 * C_Hillmann;
    const double LArg_j		 = 100;

    // NO pathway

	const double LArg        = 100;
	const double V_spine     = 8e-8;
	const double k_ex        = 1600;
	const double Ca_rest     = 0.1;
	const double lambda      = 20;
	const double V_maxNOS    = 25e-3;
	const double V_max_NO_n  = 4.22;
	const double K_mO2_n 	 = 243;
	const double K_mArg_n	 = 1.5;
	const double K_actNOS    = 9.27e-2;
	const double D_NO 	     = 3300;
	const double k_O2        = 9.6e-6;
	const double On          = 200;
	const double v_n         = -0.04;
	const double Ok          = 200;
	const double G_M         = 46000;
	const double dist_nk     = 25;
	const double dist_ki     = 25;
	const double dist_ij     = 3.75;
	const double tau_nk      = pow(dist_nk,2)/(2*D_NO);
	const double tau_ki      = pow(dist_ki,2)/(2*D_NO);
	const double tau_ij      = pow(dist_ij,2)/(2*D_NO);
	const double P_Ca_P_M    = 3.6;
	const double Ca_ex       = 2e3;
	const double M           = 1.3e5;
	const double betA        = 650 ;
	const double betB        = 2800 ;
	const double Oj          = 200;
	const double K_dis       = 9e-2;
	const double K_eNOS      = 4.5e-1;
	const double mu2         = 0.0167;
	const double g_max       = 0.06;
	const double alp         = 2;
	const double W_0         = 1.4;
	const double delt_wss    = 2.86;
	const double k_dno       = 0.01;
	const double k1          = 2e3 ;
	const double k2          = 0.1;
	const double k3          = 3;
	const double k_1         = 100;
	const double V_max_sGC   = 0.8520;  //\muM s{-1}; (for m = 2)
	const double k_pde       = 0.0195;// s{-1} (for m = 2)
	const double C_4         = 0.011; // [s{-1} microM{-2}] (note: the changing units are correct!) (for m = 2)
	const double K_m_pde     = 2;           		// [microM]
	const double k_mlcp_b    = 0.0086;         // [s{-1}]
	const double k_mlcp_c    = 0.0327;          //[s{-1}]
	const double K_m_mlcp    = 5.5;        		// [microM]
	const double bet_i       = 0.13; // translation factor for membrane potential dependence of KCa channel activation sigmoidal [microM2]
	const double m 			 = 2;
	const double gam_eNOS    = 0.1; // [-]
	const double K_mO2_j     = 7.7;
	const double V_NOj_max   = 1.22;
	const double K_mArg_j    = 1.5;

	// AC Ca2+
	const double r_buff			= 0.05;
	const double G_TRPV_k		= 50;
	const double g_TRPV_k   	= G_TRPV_k * 1e-12 / A_ef_k;
	const double J_max			= 2880;
	const double K_act			= 0.17;
	const double K_I			= 0.03;
	const double P_L			= 0.0804;
	const double k_pump			= 0.24;
	const double V_max			= 20;
	const double C_astr_k		= 40;
	const double gamma_k		= 834.3;
	const double B_ex 			= 11.35;
	const double BK_end			= 40;
	const double K_ex			= 0.26;
	const double delta			= 1.235e-2;
	const double K_G			= 8.82;
	const double Ca_3			= 0.4;
	const double Ca_4			= 0.35;
	const double v_5			= 15e-3;
	const double v_7			= -55e-3;
	const double eet_shift		= 2e-3;
	const double gam_cae_k		= 200;
	const double gam_cai_k		= 0.01;
	const double R_0_passive_k	= 20e-6;
	const double epshalf_k		= 0.1;
	const double kappa_k		= 0.1;
	const double v1_TRPV_k		= 0.12;
	const double v2_TRPV_k		= 0.013;
	const double t_TRPV_k		= 0.9;
	const double VR_ER_cyt		= 0.185;
	const double K_inh			= 0.1;
	const double k_on			= 2;
	const double k_deg			= 1.25;
	const double r_h			= 4.8;
	const double Ca_k_min		= 0.1;
	const double k_eet			= 7.2;
	const double V_eet			= 72;
	const double Ca_decay_k		= 0.5;
	const double Capmin_k		= 2000;
	const double reverseBK		= 0;
	const double switchBK		= 1;
	const double trpv_switch	= 1;
	const double z_Ca			= 2;
	const double m_c			= 4;

	// Glutamate constants
	const double Glu_max		= 1846;
	const double Glu_slope 		= 0.1;
	const double Ke_switch		= 5.5;

	// Neuron constants
	const double SC_coup	 = 11.5;
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
	const double O2_0		 = 0.02;
	const double alpha_O2     = 0.05;
	const double D_Na 		 = 1.33e-5;
	const double D_K 		 = 1.96e-5;
	const double K_init_e 	 = 2.9;
	const double Na_init_sa  = 10;
	const double Na_init_d 	 = 10;
	const double R_init	     = 1.9341e-5;
	const double CBF_init 	 = 0.032;
	const double O2_b 		 = 0.04;
	const double gamma_O2    = 0.1;
	const double Mg			 = 1.2;

	const double gNaleak_sa	 = 6.2378e-5;
	const double gKleak_sa	 = 2.1989e-4;
	const double gleak_sa	 = 10*6.2378e-5;
	const double gNaleak_d	 = 6.2961e-5;
	const double gKleak_d	 = 2.1987e-4;
	const double gleak_d	 = 10*6.2961e-5;
	const double Imax		 = 0.013*6;

	// BOLD constants
	const double tau_MTT	= 3;
	const double tau_TAT	= 20;
	const double d			= 4;
	const double a_1		= 3.4;
	const double a_2		= 1;
	const double V_0		= 0.03;
	const double E_0		= 0.4;

	// Steady state values used for normalisation of BOLD signal and change in CBF
		// J_PLC = 0.11
		const double DHG_0		= 0.6662;
		const double CBV_0		= 1.317;
		const double CBF_0		= 0.0637;

		// J_PLC = 0.3
//		const double DHG_0		= 1.0753;
//		const double CBV_0		= 0.9703;
//		const double CBF_0		= 0.0295;

	/****** Model fluxes/algebraic variables and state variables ******/

    double trans_p, trans_P_mmHg, delta_p; // pressure stuff

    // Initialise state variables
    double state_r;
    double state_R_k, state_N_Na_k, state_N_K_k, state_N_HCO3_k, state_N_Cl_k, state_N_Na_s, state_N_K_s, state_N_HCO3_s, state_K_p, state_w_k; // AC state
    double state_ca_i, state_ca_sr_i, state_v_i, state_w_i, state_ip3_i; // SMC state
    double state_ca_j, state_ca_er_j, state_v_j, state_ip3_j; // EC state
    double state_Mp, state_AM, state_AMp; // Mech state
    double state_ca_n, state_nNOS, state_NOn, state_NOi, state_E_b, state_E_6c, state_cGMP, state_eNOS, state_NOj, state_NOk; // NO pathway state
    double state_ca_k, state_s_k, state_h_k, state_ip3_k, state_eet_k, state_m_k, state_ca_p; // AC Ca2+ state
    double state_v_sa, state_v_d, state_K_sa, state_Na_sa, state_K_d, state_Na_d, state_K_e, state_Na_e; // NE ions state
	double state_Buff_e, state_O2, state_CBV, state_DHG; // NE other state
	double state_m1, state_m2, state_m3, state_m4, state_m5, state_m6, state_m7, state_m8, state_h1, state_h2, state_h3, state_h4, state_h5, state_h6; // NE gating state

    // Fluxes
    double R_s, flu_N_Cl_s, flu_Na_k, flu_K_k, flu_HCO3_k, flu_Cl_k, flu_Na_s, flu_K_s, flu_HCO3_s, flu_Cl_s, flu_E_Na_k, flu_E_K_k, flu_E_Cl_k, flu_E_NBC_k, flu_E_BK_k, flu_J_NaK_k, flu_v_k, flu_J_KCC1_k, flu_J_NBC_k, flu_J_NKCC1_k, flu_J_Na_k, flu_J_K_k, flu_J_BK_k, flu_w_inf, flu_phi_w; // AC fluxes
    double flu_M, flu_h_r, flu_v_cpl_i, flu_c_cpl_i, flu_I_cpl_i, flu_rho_i, flu_ip3_i, flu_SRuptake_i, flu_CICR_i, flu_extrusion_i, flu_leak_i, flu_VOCC_i, flu_NaCa_i, flu_NaK_i, flu_Cl_i, flu_K_i, flu_Kactivation_i, flu_degrad_i, flu_v_KIR_i, flu_G_KIR_i, flu_J_KIR_i, flu_J_stretch_i; // SMC fluxes
    double flu_v_cpl_j, flu_c_cpl_j, flu_I_cpl_j, flu_rho_j, flu_O_j, flu_ip3_j, flu_ERuptake_j, flu_CICR_j, flu_extrusion_j, flu_leak_j, flu_cation_j, flu_BKCa_j, flu_SKCa_j, flu_K_j, flu_R_j, flu_degrad_j, flu_J_stretch_j; // EC fluxes
    double flu_K1_c, flu_K6_c; // Mech fluxes
    double flu_c_w, flu_P_NR2AO, flu_P_NR2BO, flu_I_Ca, flu_CaM, flu_W_tau_w, flu_F_tau_w, flu_k4, flu_R_cGMP2, flu_K2_c, flu_K5_c, flu_tau_w, E_5c, V_max_pde;    // NO pathway fluxes
    double rho, flu_ip3_k, flu_er_leak, flu_pump_k, flu_I_TRPV_k, flu_TRPV_k, flu_E_TRPV_k, B_cyt, G, v_3, H_Ca_k, eta, minf_k, t_Ca_k, flu_VOCC_k;  // AC Ca2+ fluxes
    double F_r, E, R_0, state_R_dim;
    double E_Na_sa, E_K_sa, E_Na_d, E_K_d;
    double J_Naleak_sa, J_Kleak_sa, J_Naleak_d, J_Kleak_d;
    double m1alpha, m1beta, h1alpha, h1beta, J_NaP_sa;
    double m8alpha, m8beta, h6alpha, h6beta, J_NaT_sa;
    double m2alpha, m2beta, J_KDR_sa;
    double m3alpha, m3beta, h2alpha, h2beta, J_KA_sa;
    double m4alpha, m4beta, h3alpha, h3beta, J_NaP_d;
    double m5alpha, m5beta, h4alpha, h4beta, J_NMDA_K_d, J_NMDA_Na_d;
    double m6alpha, m6beta, J_KDR_d;
    double m7alpha, m7beta, h5alpha, h5beta, J_KA_d;
    double J_pump1_sa, J_pump1init_sa, J_pump1_d, J_pump1init_d, J_pump2, J_pump_sa, J_pump_d, J_Napump_sa, J_Kpump_sa, J_Napump_d, J_Kpump_d, J_pump2_0, J_pump2_O2_0;
    double J_Na_tot_sa, J_K_tot_sa, J_leak_tot_sa, J_Na_tot_d, J_K_tot_d, J_leak_tot_d, J_tot_sa, J_tot_d;
    double P_02, CBF, J_O2_vascular, J_O2_background, J_O2_pump;
    double f_out, CMRO2, CMRO2_init, OEF, BOLD, Delta_CBF;
    double Glu;


    // State Variables:
    state_r  	  = u[i_radius];

    state_ca_i    = u[ca_i];
    state_ca_sr_i = u[ca_sr_i];
    state_v_i     = u[v_i];
    state_w_i     = u[w_i];
    state_ip3_i   = u[ip3_i];

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

    state_ca_n    = u[ca_n];
    state_nNOS    = u[nNOS];
    state_NOn     = u[NOn];
    state_NOi     = u[NOi];
    state_E_b     = u[E_b];
    state_E_6c    = u[E_6c];
    state_cGMP    = u[cGMP];
    state_eNOS    = u[eNOS];
    state_NOj     = u[NOj];
    state_NOk     = u[NOk];

    state_ca_k     = u[ca_k];
    state_s_k      = u[s_k];
    state_h_k      = u[h_k];
    state_ip3_k    = u[ip3_k];
    state_eet_k    = u[eet_k];
    state_m_k      = u[m_k];
    state_ca_p     = u[ca_p];

    state_v_sa		= u[v_sa];
    state_v_d		= u[v_d];
	state_K_sa		= u[K_sa];
	state_Na_sa		= u[Na_sa];
	state_K_d		= u[K_d];
	state_Na_d		= u[Na_d];
	state_K_e		= u[K_e];
	state_Na_e		= u[Na_e];
	state_Buff_e	= u[Buff_e];
	state_O2		= u[O2];
	state_CBV		= u[CBV];
	state_DHG		= u[DHG];
	state_m1		= u[m1];
	state_m2		= u[m2];
	state_m3		= u[m3];
	state_m4		= u[m4];
	state_m5		= u[m5];
	state_m6		= u[m6];
	state_m7		= u[m7];
	state_m8		= u[m8];
	state_h1		= u[h1];
	state_h2		= u[h2];
	state_h3		= u[h3];
	state_h4		= u[h4];
	state_h5		= u[h5];
	state_h6		= u[h6];

// Fluxes:

	// Dimensional radius
	state_R_dim			= state_r * R_0_passive_k;

    // pressure
    trans_p 			= P0 / 2 * (p + nvu_w->pcap); // x P0 to make dimensional, transmural pressure in Pa
	trans_P_mmHg 		= trans_p * PA2MMHG; // transmural pressure in mmHg
	delta_p				= P0 * (p - nvu_w->pcap); // dimensional pressure drop over leaf vessel

    // SC fluxes
    R_s        	    	= R_tot - state_R_k;                            // state_R_k is AC volume-area ratio, R_s is SC
    flu_N_Cl_s         	= state_N_Na_s + state_N_K_s - state_N_HCO3_s;  //
    flu_Na_s           	= state_N_Na_s / R_s;                       //
    flu_K_s            	= state_N_K_s / R_s;                        //
    flu_HCO3_s         	= state_N_HCO3_s / R_s;                     //
    flu_Cl_s           	= flu_N_Cl_s / R_s;                         //

    // AC fluxes
	flu_E_TRPV_k	= R_gas * Temp / (z_Ca * F) * log(state_ca_p / state_ca_k); // TRPV4 channel Nernst Potential

    flu_Na_k           	= state_N_Na_k / state_R_k;                     //
    flu_K_k            	= state_N_K_k / state_R_k;                      //
    flu_HCO3_k         	= state_N_HCO3_k / state_R_k;                   //
    flu_Cl_k           	= state_N_Cl_k / state_R_k;                     //
    flu_E_Na_k         	= (R_gas * Temp) / (z_Na * F) * log(flu_Na_s / flu_Na_k);    // V
    flu_E_K_k          	= (R_gas * Temp) / (z_K  * F) * log(flu_K_s / flu_K_k );     // V
    flu_E_Cl_k         	= (R_gas * Temp) / (z_Cl * F) * log(flu_Cl_s / flu_Cl_k);    // V
    flu_E_NBC_k        	= (R_gas * Temp) / (z_NBC* F) * log((flu_Na_s * pow(flu_HCO3_s,2))/(flu_Na_k * pow(flu_HCO3_k,2)));     // V
    flu_E_BK_k         	= reverseBK + switchBK * ((R_gas * Temp) / (z_K  * F) * log(state_K_p / flu_K_k));   // V
    flu_J_NaK_k        	= J_NaK_max * ( pow(flu_Na_k,1.5) / ( pow(flu_Na_k,1.5) + pow(K_Na_k,1.5) ) ) * ( flu_K_s / (flu_K_s + K_K_s) );    // uMm s-1
    flu_v_k         	= ( g_Na_k * flu_E_Na_k + g_K_k * flu_E_K_k + g_TRPV_k * state_m_k * flu_E_TRPV_k + g_Cl_k * flu_E_Cl_k + g_NBC_k * flu_E_NBC_k + g_BK_k * state_w_k * flu_E_BK_k - flu_J_NaK_k * F / unitcon ) / ( g_Na_k + g_K_k + g_Cl_k + g_NBC_k + g_TRPV_k * state_m_k + g_BK_k * state_w_k );
    flu_J_KCC1_k       	= (R_gas * Temp * g_KCC1_k) / (pow(F,2)) * log(((flu_K_s) * (flu_Cl_s))/((flu_K_k)*(flu_Cl_k))) * unitcon;   //uMm s-1
    flu_J_NBC_k        	= g_NBC_k / F * (flu_v_k - flu_E_NBC_k) * unitcon;       //uMm s-1
    flu_J_NKCC1_k     	= (g_NKCC1_k * R_gas * Temp) / (pow(F,2))  * log(((flu_K_s) * (flu_Na_s) * pow(flu_Cl_s,2)) /((flu_K_k) * (flu_Na_k) * pow(flu_Cl_k,2)))*unitcon;        //uMm s-1
    flu_J_Na_k   		= g_Na_k / F * (flu_v_k - flu_E_Na_k) * unitcon;              //uMm s-1
    flu_J_K_k    		= g_K_k  / F * ((flu_v_k) - (flu_E_K_k )) * unitcon;          //uMm s-1
    flu_J_BK_k   		= g_BK_k / F * state_w_k * (flu_v_k - flu_E_BK_k) * unitcon;  //uMm s-1

    // SMC fluxes
    flu_M               = 1 - state_Mp - state_AM - state_AMp;
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
    flu_degrad_i	    = k_i * state_ip3_i;
    flu_J_stretch_i     = G_stretch/(1 + exp( -alpha1 * (trans_P_mmHg * state_r / flu_h_r - sig0))) * (state_v_i - Esac);
    flu_v_KIR_i    		= z_1 * state_K_p / unitcon + z_2;                                  // mV           state_K_p,
    flu_G_KIR_i    		= exp( z_5 * state_v_i + z_3 * state_K_p / unitcon + z_4 );        // pS pF-1 =s-1  state_v_i, state_K_p
    flu_J_KIR_i    		= F_il/gam * (flu_G_KIR_i) * (state_v_i - (flu_v_KIR_i));            // mV s-1 //     state_v_i, state_K_p
    
    // EC fluxes
    flu_v_cpl_j			= - g_hat * ( state_v_j - state_v_i );
    flu_c_cpl_j			= - p_hat * ( state_ca_j - state_ca_i );
    flu_I_cpl_j			= - p_hatIP3 * ( state_ip3_j - state_ip3_i );
    flu_rho_j 			= 1;
    flu_O_j 			= JO_j;
    flu_ip3_j			= Fmax_j * ( pow(state_ip3_j, 2) ) / ( pow(Kr_j, 2) + pow(state_ip3_j, 2) );
    flu_ERuptake_j      = B_j * ( pow(state_ca_j, 2) ) / ( pow(state_ca_j, 2) + pow(cb_j, 2) );
    flu_CICR_j			= C_j *  ( pow(state_ca_er_j, 2) ) / ( pow(sc_j, 2) + pow(state_ca_er_j, 2) ) *  ( pow(state_ca_j, 4) ) / ( pow(cc_j,4) + pow(state_ca_j,4) );
    flu_extrusion_j     = D_j * state_ca_j;
    flu_leak_j          = L_j * state_ca_er_j;
    flu_cation_j 		= G_cat * ( E_Ca - state_v_j) * 0.5 * ( 1 + tanh( ( log10( state_ca_j ) - m3cat )  /  m4cat  ) );
    flu_BKCa_j 			= 0.2 * ( 1 + tanh( ( (  log10(state_ca_j) - c) * ( state_v_j - b ) - a1 ) / ( m3b* pow(( state_v_j + a2 * ( log10( state_ca_j ) - c ) - b),2) + m4b ) ) );
    flu_SKCa_j 			= 0.3 * ( 1 + tanh( ( log10(state_ca_j) - m3s ) /  m4s ));
    flu_K_j 			= G_tot * ( state_v_j - vK_j ) * ( flu_BKCa_j + flu_SKCa_j );
    flu_R_j 			= G_R * ( state_v_j - v_rest);
    flu_degrad_j 		= k_j * state_ip3_j;
    flu_J_stretch_j     = G_stretch / (1 + exp(-alpha1*(trans_P_mmHg * state_r / flu_h_r - sig0))) * (state_v_j - Esac);

// Mech fluxes
    flu_K1_c         	= gam_cross * pow(state_ca_i,3);
    flu_K6_c        	= flu_K1_c;
    F_r			        = state_AMp + state_AM;
    E 			        = EPASSIVE + F_r * (EACTIVE - EPASSIVE);
    R_0 			    = R_0_passive_k + F_r * (RSCALE - 1) * R_0_passive_k;

// NO pathway fluxes

	Glu 			= 0.5 * Glu_max * ( 1 + tanh( (state_K_e - Ke_switch) / Glu_slope) );

	flu_P_NR2AO         = Glu / (betA + Glu);
	flu_P_NR2BO         = Glu / (betB + Glu);
	flu_I_Ca            = (-4 * v_n * G_M * P_Ca_P_M * (Ca_ex/M)) / (1+exp(-80*(v_n+0.02))) * (exp(2 * v_n * F / (R_gas * Temp))) / (1 - exp(2 * v_n * F / (R_gas * Temp))) * (0.63 * flu_P_NR2AO + 11 * flu_P_NR2BO);
	flu_CaM             = state_ca_n / m_c;                                      // concentration of calmodulin / calcium complexes ; (100)
    flu_tau_w			= state_R_dim * delta_p / (2*200e-6); // WSS using pressure from the H tree. L_0 = 200 um

	flu_W_tau_w         = W_0 * pow((flu_tau_w + sqrt(16 * pow(delt_wss,2) + pow(flu_tau_w,2)) - 4 * delt_wss),2) / (flu_tau_w + sqrt(16 * pow(delt_wss,2) + pow(flu_tau_w,2))) ;  // - tick
	flu_F_tau_w         = 1 / (1 + alp * exp(-flu_W_tau_w)) - 1 / (1 + alp); // -(1/(1+alp)) was added to get no NO at 0 wss (!) - tick
	flu_k4              = C_4 * pow(state_cGMP, m);
	flu_R_cGMP2         = (pow(state_cGMP, 2)) / (pow(state_cGMP, 2) + pow(K_m_mlcp, 2)); // - tick
	flu_K2_c            = 58.1395 * k_mlcp_b + 58.1395 * k_mlcp_c * flu_R_cGMP2;  // Factor is chosen to relate two-state model of Yang2005 to Hai&Murphy model
	flu_K5_c            = flu_K2_c;
	flu_c_w             = 1/2 * ( 1 + tanh( (state_cGMP - 10.75) / 0.668 ) );

	flu_Kactivation_i   = pow((state_ca_i + flu_c_w),2) / (pow((state_ca_i + flu_c_w),2) + bet_i * exp(-(state_v_i - v_Ca3) / R_K));
	E_5c				= 1 - state_E_b - state_E_6c;
	V_max_pde			= k_pde * state_cGMP;

// AC Ca2+ fluxes

	rho				= 0.1 + 0.6/1846 * Glu;

    flu_ip3_k 		= J_max * pow(( state_ip3_k / ( state_ip3_k + K_I ) * state_ca_k / ( state_ca_k + K_act ) * state_h_k ) , 3) * (1.0 - state_ca_k / state_s_k);
	flu_er_leak 	= P_L * ( 1.0 - state_ca_k / state_s_k );
	flu_pump_k  	= V_max * pow(state_ca_k, 2) / ( pow(state_ca_k, 2) + pow(k_pump, 2) );
	flu_I_TRPV_k	= G_TRPV_k * state_m_k * (flu_v_k - flu_E_TRPV_k) * unitcon;
	flu_TRPV_k		= -0.5 * flu_I_TRPV_k / ( C_astr_k * gamma_k );
	B_cyt 			= 1.0 / (1.0 + BK_end + K_ex * B_ex / pow((K_ex + state_ca_k), 2) );
	G				= ( rho + delta ) / ( K_G + rho + delta );
	v_3				= -v_5 / 2.0 * tanh((state_ca_k - Ca_3) / Ca_4) + v_7;
	flu_w_inf    	= 0.5 * ( 1 + tanh( ( flu_v_k + eet_shift * state_eet_k - v_3 ) / v_4 ) );
    flu_phi_w    	= psi_w * cosh( (flu_v_k - v_3) / (2*v_4) );
    H_Ca_k			= state_ca_k / gam_cai_k + state_ca_p / gam_cae_k;
    eta 			= (state_R_dim - R_0_passive_k) / R_0_passive_k;
    minf_k 			= ( 1 / ( 1 + exp( - (eta - epshalf_k) / kappa_k ) ) ) * ( ( 1 / (1 + H_Ca_k) ) * (H_Ca_k + tanh(( flu_v_k - v1_TRPV_k) / v2_TRPV_k )));
    t_Ca_k 			= t_TRPV_k / state_ca_p;
    flu_VOCC_k		= flu_VOCC_i;

// Neuron fluxes
    // Nernst potentials //
    E_Na_sa 		= ph * log(state_Na_e / state_Na_sa);
	E_K_sa 			= ph * log(state_K_e / state_K_sa);
	E_Na_d 			= ph * log(state_Na_e / state_Na_d);
	E_K_d 			= ph * log(state_K_e / state_K_d);

	// Ion fluxes //
	// Leak fluxes
	J_Naleak_sa 	= gNaleak_sa * (state_v_sa - E_Na_sa);
	J_Kleak_sa 		= gKleak_sa * (state_v_sa - E_K_sa);
	J_Naleak_d 		= gNaleak_d * (state_v_d - E_Na_d);
	J_Kleak_d 		= gKleak_d * (state_v_d - E_K_d);

	// Na+ flux through NaP in soma/axon
	m1alpha 		= 1 / (6 * (1 + exp(-((0.143 * state_v_sa) + 5.67))));
	m1beta 			= exp(-((0.143 * state_v_sa) + 5.67)) / (6 * (1 + exp(-((0.143 * state_v_sa) + 5.67))));
	h1alpha 		= 5.12e-8 * exp(-((0.056 * state_v_sa) + 2.94));
	h1beta 			= 1.6e-6 / (1 + exp(-(((0.2 * state_v_sa)) + 8)));
	J_NaP_sa 		= (pow(m1,2) * h1 * gNaP_GHk * Farad * state_v_sa * (state_Na_sa - (exp(-state_v_sa / ph) * state_Na_e))) / (ph * (1 - exp(-state_v_sa / ph)));

	// Na+ flux through NaT in soma/axon
	m8alpha 		= 0.32 * ((-state_v_sa - 51.9) / (exp(-(0.25 * state_v_sa + 12.975)) - 1));
	m8beta 			= 0.28 * ((state_v_sa + 24.89) / (exp(0.2 * state_v_sa + 4.978) - 1));
	h6alpha 		= 0.128 * exp(-(0.056 * state_v_sa + 2.94));
	h6beta 			= 4 / (1 + exp(-(0.2 * state_v_sa + 6)));
	J_NaT_sa 		= (pow(m8,3) * h6 * gNaT_GHk * Farad * state_v_sa * (state_Na_sa - (exp(-state_v_sa / ph) * state_Na_e))) / (ph * (1 - exp(-state_v_sa / ph)));

	// K+ flux through KDR in soma/axon
	m2alpha 		= 0.016 * ((state_v_sa + 34.9) / (1 - exp(-((0.2 * state_v_sa) + 6.98))));
	m2beta 			= 0.25 * exp(-((0.025 * state_v_sa) + 1.25));
	J_KDR_sa 		= (pow(m2,2) * gKDR_GHk * Farad * state_v_sa * (state_K_sa - (exp(-state_v_sa / ph) * state_K_e))) / (ph * (1 - exp(-state_v_sa / ph)));

	// K+ flux through KA in soma/axon
	m3alpha 		= 0.02 * ((state_v_sa + 56.9) / (1 - exp(-((0.1 * state_v_sa) + 5.69))));
	m3beta 			= 0.0175 * ((state_v_sa + 29.9) / (exp(((0.1 * state_v_sa) + 2.99)) - 1));
	h2alpha 		= 0.016 * exp(-((0.056 * state_v_sa) + 4.61));
	h2beta 			= 0.5 / (1 + exp(-((0.2 * state_v_sa) + 11.98)));
	J_KA_sa 		= (pow(m3,2) * h2 * gKA_GHk * Farad * state_v_sa * (state_K_sa - (exp(-state_v_sa / ph) * state_K_e))) / (ph * (1 - exp(-state_v_sa / ph)));

	// Na+ flux through NaP in dendrite
	m4alpha 		= 1 / (6 * (1 + exp(-((0.143 * state_v_d) + 5.67))));
	m4beta 			= exp(-((0.143 * state_v_d) + 5.67)) / (6 * (1 + exp(-((0.143 * state_v_d) + 5.67))));
	h3alpha 		= 5.12e-8 * exp(-((0.056 * state_v_d) + 2.94));
	h3beta 			= 1.6e-6 / (1 + exp(-(((0.2 * state_v_d)) + 8)));
	J_NaP_d 		= (pow(m4,2) * h3 * gNaP_GHk * Farad * state_v_d * (state_Na_d - (exp(-state_v_d / ph) * state_Na_e))) / (ph * (1 - exp(-state_v_d / ph)));

	// Na+ and K+ flux through NMDA in dendrite
	m5alpha 		= 0.5 / (1 + exp((13.5 - state_K_e) / 1.42));
	m5beta 			= 0.5 - m5alpha;
	h4alpha 		= 1 / (2000 * (1 + exp((state_K_e - 6.75) / 0.71)));
	h4beta 			= 5e-4 - h4alpha;
	J_NMDA_K_d 		= ( (m5 * h4 * gNMDA_GHk * Farad * state_v_d * (state_K_d - (exp(-state_v_d / ph) * state_K_e))) / (ph * (1 - exp(-state_v_d / ph))) ) / (1 + 0.33 * Mg * exp(-(0.07 * state_v_d + 0.7)));
	J_NMDA_Na_d 	= ( (m5 * h4 * gNMDA_GHk * Farad * state_v_d * (state_Na_d - (exp(-state_v_d / ph) * state_Na_e))) / (ph * (1 - exp(-state_v_d / ph))) ) / (1 + 0.33 * Mg * exp(-(0.07 * state_v_d + 0.7)));

	// K+ flux through KDR in dendrite
	m6alpha 		= 0.016 * ((state_v_d + 34.9) / (1 - exp(-((0.2 * state_v_d) + 6.98))));
	m6beta 			= 0.25 * exp(-((0.025 * state_v_d) + 1.25));
	J_KDR_d 		=(pow(m6,2) * gKDR_GHk * Farad * state_v_d * (state_K_d - (exp(-state_v_d / ph) * state_K_e))) / (ph * (1 - exp(-state_v_d / ph)));

	// K+ flux through KA in dendrite
	m7alpha 		= 0.02 * ((state_v_d + 56.9) / (1 - exp(-((0.1 * state_v_d) + 5.69))));
	m7beta 			= 0.0175 * ((state_v_d + 29.9) / (exp(((0.1 * state_v_d) + 2.99)) - 1));
	h5alpha 		= 0.016 * exp(-((0.056 * state_v_d) + 4.61));
	h5beta 			= 0.5 / (1 + exp(-((0.2 * state_v_d) + 11.98)));
	J_KA_d 			= (pow(m7,2) * h5 * gKA_GHk * Farad * state_v_d * (state_K_d - (exp(-state_v_d / ph) * state_K_e))) / (ph * (1 - exp(-state_v_d / ph)));

	// ATPase pump
	J_pump1_sa 		= pow((1 + (K_init_e / state_K_e)),-2) * pow((1 + (Na_init_sa / state_Na_sa)),-3);
	J_pump1init_sa 	= 0.0312;				// (1 + (K_init_e / K_init_e))^(-2) * (1 + (Na_init_sa / Na_init_sa)) ^ (-3);
	J_pump1_d 		= pow((1 + (K_init_e / state_K_e)),-2) * pow((1 + (Na_init_d / state_Na_d)),-3);
	J_pump1init_d 	= 0.0312;				// (1 + (K_init_e / K_init_e))^(-2) * (1 + (Na_init_d / Na_init_d)) ^ (-3);
	J_pump2 		= 2 * pow((1 + O2_0 / (((1 - alpha_O2) * state_O2) + alpha_O2 * O2_0)),-1);
	J_pump_sa 		= Imax * J_pump1_sa * J_pump2;
	J_pump_d 		= Imax * J_pump1_d * J_pump2;
	J_Napump_sa 	= 3 * J_pump_sa;
	J_Kpump_sa 		= -2 * J_pump_sa;
	J_Napump_d 		= 3 * J_pump_d;
	J_Kpump_d 		= -2 * J_pump_d;

	// Total fluxes
	J_Na_tot_sa 	= J_NaP_sa + J_Naleak_sa + J_Napump_sa + J_NaT_sa;
	J_K_tot_sa 		= J_KDR_sa + J_KA_sa + J_Kleak_sa + J_Kpump_sa;
	J_leak_tot_sa 	= gleak_sa * (state_v_sa - E_Cl_sa);
	J_Na_tot_d 		= J_NaP_d + J_Naleak_d + J_Napump_d + J_NMDA_Na_d;
	J_K_tot_d 		= J_KDR_d + J_KA_d + J_Kleak_d + J_Kpump_d + J_NMDA_K_d;
	J_leak_tot_d 	= gleak_d * (state_v_d - E_Cl_d);
	J_tot_sa 		= J_Na_tot_sa + J_K_tot_sa + J_leak_tot_sa;
	J_tot_d 		= J_Na_tot_d + J_K_tot_d + J_leak_tot_d;

	// Oxygen //
	J_pump2_0		= 0.0952; 	// 2 * (1 + O2_0 / (((1 - alpha_O2) * 0) + alpha_O2 * O2_0))^(-1)
	J_pump2_O2_0	= 1; 		// 2 * (1 + O2_0 / (((1 - alpha_O2) * O2_0) + alpha_O2 * O2_0))^(-1)
	P_02 			= (J_pump2 - J_pump2_0 ) / ( J_pump2_O2_0 - J_pump2_0);
	CBF 			= CBF_init * (pow(state_R_dim,4) / pow(R_init,4));
	J_O2_vascular 	= CBF * ((O2_b - state_O2) / (O2_b - O2_0));
	J_O2_background = CBF_init * P_02 * (1 - gamma_O2);
	J_O2_pump 		= CBF_init * P_02 * gamma_O2 * ((J_pump1_sa + J_pump1_d) / (J_pump1init_sa + J_pump1init_d));

	// BOLD //
	f_out 			= pow(state_CBV,(1/d)) + tau_TAT * (1/(tau_MTT + tau_TAT) * ( CBF/CBF_init  - pow(state_CBV,(1/d)) ));
	CMRO2 			= J_O2_background + J_O2_pump;
	CMRO2_init 		= CBF_init * P_02;
	OEF 			= CMRO2 * E_0 / CBF;

	BOLD 			= 100 * V_0 * ( a_1 * (1 - state_DHG/DHG_0) - a_2 * (1 - state_CBV/CBV_0) ); // divides DHG, CBV and CBF by initial values, will depend on J_PLC
	Delta_CBF		= (CBF - CBF_0)/CBF_0;


// Differential Equations:

    /***********Neuron dynamics*********/

	// Neuron - other
    du[Buff_e]   = Mu * state_K_e * (B0 - state_Buff_e) / (1 + exp(-((state_K_e - 5.5) / 1.09))) - (Mu * state_Buff_e); // [mM]?
    du[O2]       = J_O2_vascular - J_O2_background - J_O2_pump; 						// [mM]
    du[CBV]      = 1/(tau_MTT + tau_TAT) * ( CBF/CBF_init  - pow(state_CBV,(1/d)) );	// [-]
    du[DHG]      = 1/tau_MTT * ( CMRO2/CMRO2_init - state_DHG/state_CBV * f_out );		// [-]

    // Neuron - ions (mV or mM)
    du[v_sa]	 = 1/Cm * ( -J_tot_sa + 1 / (2 * Ra * pow(dhod,2)) * (state_v_d - state_v_sa) + current_input(t,x,y) );
	du[v_d]      = 1/Cm * (-J_tot_d + 1 / (2 * Ra * pow(dhod,2)) * (state_v_sa - state_v_d));
    du[K_sa]     = -As / (Farad * Vs) * J_K_tot_sa + D_K * (Vd + Vs) / (2 * pow(dhod,2) * Vs) * (state_K_d - state_K_sa);
    du[Na_sa]    = -As / (Farad * Vs) * J_Na_tot_sa + D_Na * (Vd + Vs) / (2 * pow(dhod,2) * Vs) * (state_Na_d - state_Na_sa);
    du[K_d]      = -Ad / (Farad * Vd) * J_K_tot_d + D_K * (Vs + Vd) / (2 * pow(dhod,2) * Vd) * (state_K_sa - state_K_d);
    du[Na_d]     = -Ad / (Farad * Vd) * J_Na_tot_d + D_Na * (Vs + Vd) / (2 * pow(dhod,2) * Vd) * (state_Na_sa - state_Na_d);
    du[K_e]      = 1/(Farad * fe) * (((As * J_K_tot_sa) / Vs)  + ((Ad * J_K_tot_d) / Vd)) - du[Buff_e];
    du[Na_e]     = 1/(Farad * fe) * (((As * J_Na_tot_sa) / Vs) + ((Ad * J_Na_tot_d) / Vd));

	// Neuron - gating variables ([-])
    du[m1]     	= 1000 * ((m1alpha * (1 - state_m1)) - (m1beta * state_m1));
    du[m2]     	= 1000 * ((m2alpha * (1 - state_m2)) - (m2beta * state_m2));
    du[m3]     	= 1000 * ((m3alpha * (1 - state_m3)) - (m3beta * state_m3));
    du[m4]     	= 1000 * ((m4alpha * (1 - state_m4)) - (m4beta * state_m4));
    du[m5]     	= 1000 * ((m5alpha * (1 - state_m5)) - (m5beta * state_m5));
    du[m6]     	= 1000 * ((m6alpha * (1 - state_m6)) - (m6beta * state_m6));
    du[m7]     	= 1000 * ((m7alpha * (1 - state_m7)) - (m7beta * state_m7));
    du[m8]     	= 1000 * ((m8alpha * (1 - state_m8)) - (m8beta * state_m8));
    du[h1]     	= 1000 * ((h1alpha * (1 - state_h1)) - (h1beta * state_h1));
    du[h2]     	= 1000 * ((h2alpha * (1 - state_h2)) - (h2beta * state_h2));
    du[h3]     	= 1000 * ((h3alpha * (1 - state_h3)) - (h3beta * state_h3));
    du[h4]     	= 1000 * ((h4alpha * (1 - state_h4)) - (h4beta * state_h4));
    du[h5]     	= 1000 * ((h5alpha * (1 - state_h5)) - (h5beta * state_h5));
    du[h6]     	= 1000 * ((h6alpha * (1 - state_h6)) - (h6beta * state_h6));

	/***********General dynamics**********/

    du[i_radius		] = 1 / ETA * (state_r * trans_p / flu_h_r - E * (state_R_dim - R_0)/R_0); // Radius - nondimensional (state_R_dim: dimensional)

    //AC:
    du[ R_k     ] = 0; //L_p * (flu_Na_k + flu_K_k + flu_Cl_k + flu_HCO3_k - flu_Na_s - flu_K_s - flu_Cl_s - flu_HCO3_s + X_k / state_R_k);  // m s-1
    du[ N_Na_k  ] = -flu_J_Na_k - 3 * flu_J_NaK_k + flu_J_NKCC1_k + flu_J_NBC_k;    // uMm s-1
    du[ N_K_k   ] = -flu_J_K_k + 2 * flu_J_NaK_k + flu_J_NKCC1_k + flu_J_KCC1_k -flu_J_BK_k; // uMm s-1
    du[ N_HCO3_k] = 2 * flu_J_NBC_k;                                                // uMm s-1
    du[ N_Cl_k  ] = du[ N_Na_k] + du[ N_K_k] - du[ N_HCO3_k];                       // uMm s-1, modified equation compared to the one of Ostby  //
    du[ w_k     ] = flu_phi_w * (flu_w_inf - state_w_k);                            // s-1

    //SC:
    du[ N_Na_s  ] = - du[ N_Na_k] - SC_coup * du[K_e] * 1000 * R_s;                           // uMm s-1
    du[ N_K_s   ] = flu_J_K_k - 2 * flu_J_NaK_k - flu_J_NKCC1_k - flu_J_KCC1_k + SC_coup * du[K_e] * 1000 * R_s;                 // uMm s-1
    du[ N_HCO3_s] = - du[ N_HCO3_k];                                                // uMm s-1

    //PVS:
    du[ K_p     ] = flu_J_BK_k / (VR_pa * state_R_k) + flu_J_KIR_i / VR_ps - R_decay * (state_K_p - K_p_min);         // uM s-1

    //SMC:
    du[ ca_i    ] = flu_c_cpl_i + flu_rho_i * ( flu_ip3_i - flu_SRuptake_i + flu_CICR_i - flu_extrusion_i + flu_leak_i - flu_VOCC_i + flu_NaCa_i - 0.1* flu_J_stretch_i);
    du[ ca_sr_i ] = flu_SRuptake_i - flu_CICR_i - flu_leak_i ;
    du[ v_i     ] = flu_v_cpl_i + gam * ( - flu_NaK_i - flu_Cl_i - 2 * flu_VOCC_i - flu_NaCa_i - flu_K_i - flu_J_stretch_i - flu_J_KIR_i );
    du[ w_i     ] = lam * (flu_Kactivation_i - state_w_i ) ;
    du[ ip3_i   ] = flu_I_cpl_i - flu_degrad_i ;          // **
    du[ K_i     ] = - flu_J_KIR_i - flu_K_i + flu_NaK_i;                                            // uM s-1

    //EC:
    du[ca_j     ] = flu_c_cpl_j + flu_rho_j * ( flu_ip3_j - flu_ERuptake_j + flu_CICR_j - flu_extrusion_j + flu_leak_j + flu_cation_j + flu_O_j - flu_J_stretch_j ) ;
    du[ca_er_j  ] = flu_ERuptake_j - flu_CICR_j - flu_leak_j ;
    du[v_j      ] = flu_v_cpl_j - 1/C_m * ( flu_K_j + flu_R_j ) ;
    du[ip3_j    ] = flu_I_cpl_j + J_PLC - flu_degrad_j ;  // **

    // Mech:
    du[ Mp   	] = K4_c * state_AMp + flu_K1_c * flu_M - (flu_K2_c + K3_c) * state_Mp;
    du[ AMp  	] = K3_c * state_Mp + flu_K6_c * state_AM - (K4_c + flu_K5_c) * state_AMp;
    du[ AM   	] = flu_K5_c * state_AMp - ( K7_c + flu_K6_c ) * state_AM;

    /***********NO pathway***********/

    // NE:
    du[ca_n]       = (flu_I_Ca / (2*F * V_spine) - (k_ex * (state_ca_n - Ca_rest))) / (1 + lambda);     //\muM
    du[nNOS]       = V_maxNOS * flu_CaM / (K_actNOS + flu_CaM) - mu2 * state_nNOS ;                  //\muM
    du[NOn]       = state_nNOS * V_max_NO_n * On / (K_mO2_n + On) * LArg / (K_mArg_n + LArg) - ((state_NOn - state_NOk) / tau_nk) - (k_O2* pow(state_NOn,2) * On);

    // AC:
    du[NOk]       = (state_NOn - state_NOk) / tau_nk + (state_NOi - state_NOk) / tau_ki - k_O2 * pow(state_NOk,2) * Ok;

    // SMC:
    du[NOi]       = (state_NOk - state_NOi) / tau_ki + (state_NOj - state_NOi) / tau_ij - k_dno * state_NOi ;
    du[E_b]        = -k1 * state_E_b * state_NOi + k_1 * state_E_6c + flu_k4 * E_5c;
    du[E_6c]       = k1 * state_E_b * state_NOi - k_1 * state_E_6c - k2 * state_E_6c - k3 * state_E_6c * state_NOi ;
    du[cGMP]       = V_max_sGC * E_5c - V_max_pde * state_cGMP / (K_m_pde + state_cGMP);

    // EC:
    du[eNOS]       = gam_eNOS * (K_dis * state_ca_j / (K_eNOS + state_ca_j))  // Ca-dependent activation - tick
					+ (1 - gam_eNOS) * (g_max * flu_F_tau_w) // wss-dependent activation - tick
					- mu2 * state_eNOS;      // deactivation - tick
    du[NOj]       = V_NOj_max * state_eNOS * Oj / (K_mO2_j + Oj) * LArg_j / (K_mArg_j + LArg_j) // production - tick
					- k_O2 * pow(state_NOj,2) * Oj // consumption
					+ (state_NOi - state_NOj) / tau_ij - state_NOj * 4 * D_NO / (pow(25,2));

    /**********Astrocytic Calcium*******/

    // AC:
    du[ca_k]	= B_cyt * (flu_ip3_k - flu_pump_k + flu_er_leak + flu_TRPV_k/r_buff);
    du[s_k]		= -(B_cyt * (flu_ip3_k - flu_pump_k + flu_er_leak)) / (VR_ER_cyt);
    du[h_k]		= k_on * (K_inh - (state_ca_k + K_inh) * state_h_k);
    du[ip3_k]	= r_h * G - k_deg * state_ip3_k;
    du[eet_k]	= V_eet * fmax((state_ca_k - Ca_k_min), 0) - k_eet * state_eet_k;
    du[m_k]		= trpv_switch * (minf_k - state_m_k) / (t_Ca_k * state_ca_p);

	//PVS:
    du[ca_p]	= (-flu_TRPV_k / VR_pa) + (flu_VOCC_k / VR_ps) - Ca_decay_k * (state_ca_p - Capmin_k);


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

double factorial(int c)
{
    double result = 1;

    for (int n = 1; n < c; n++)
    {
        result = result * n;
    }

    return result;
}

// Space- & time-varying K+ input signal (simulating neuronal activity)
double current_input(double t, double x, double y)
{
    double I_strength 	= 0.025;
    double t_up   		= 1000;
    double t_down 		= 2000;

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
    if (t >= t_up && t <= t_down)
    {
         current_time = 1;
    }
    else
    {
        current_time = 0;
    }

    double current_out = I_strength * current_space * current_time; // 0 if t3 < t or x,y <= 0
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

double ECS_input(double t, double x, double y)
{
    double ECS_max 		= 9e3;
    double t_up   		= 1000;
    double t_down 		= 2000;
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

// Initial conditions. If you want spatial inhomegeneity, make it a
// function of the coordinates x and y. u0 is already allocated, you just
// need to fill in the entries
void nvu_ics(double *u0, double x, double y, nvu_workspace *nvu_w)
{
		// Current ICs optimised for J_PLC = 0.11
		u0[i_radius]  = 1.49986;

	    u0[R_k]       = 0.06e-6;
	    u0[N_Na_k]    = 0.00115629;
	    u0[N_K_k]     = 0.00554052;
	    u0[N_HCO3_k]  = 0.000582264;
	    u0[N_Cl_k]    = 0.000505576;
	    u0[w_k]       = 3.61562e-5;

	    u0[N_Na_s]    = 0.00414271;
	    u0[N_K_s]     = 7.27797e-5;
	    u0[N_HCO3_s]  = 0.000438336;

	    u0[K_p]       = 3246.44;

	    u0[ca_i]      = 0.137077;
	    u0[ca_sr_i]   = 1.20037;
	    u0[v_i]       = -58.5812;
	    u0[w_i]       = 0.38778;
	    u0[ip3_i]     = 0.45;
	    u0[K_i]       = 99994.8;

	    u0[ca_j]      = 0.537991;
	    u0[ca_er_j]   = 0.872007;
	    u0[v_j]       = -64.8638;
	    u0[ip3_j]     = 1.35;

	    u0[Mp]        = 0.0165439;
	    u0[AMp]       = 0.00434288;
	    u0[AM]        = 0.0623458;

	    // NO pathway
		u0[NOn]       = 0.273264;
		u0[NOk]       = 0.21676;
		u0[NOi]       = 0.160269;
		u0[NOj]       = 0.159001;
		u0[cGMP]      = 11.6217;
		u0[eNOS]      = 2.38751;
		u0[nNOS]      = 0.317995;
		u0[ca_n]      = 0.1;
		u0[E_b]       = 0.184071;
		u0[E_6c]      = 0.586609;

		// AC Ca2+ pathway
		u0[ca_k]     = 0.133719;
		u0[s_k]      = 502.461;
		u0[h_k]      = 0.427865;
		u0[ip3_k]    = 0.048299;
		u0[eet_k]    = 0.337187;
		u0[m_k]      = 0.896358;
		u0[ca_p]     = 1713.39;

		// Neuron - ions
		u0[v_sa]     = -70;
		u0[v_d]      = -70;
		u0[K_sa]     = 133.5;
		u0[Na_sa]    = 9.9854;
		u0[K_d]      = 133.5;
		u0[Na_d]     = 9.9853;
		u0[K_e]      = 3.5006;
		u0[Na_e]     = 139.76;

		// Neuron - other
		u0[Buff_e]   = 170;
		u0[O2]       = 0.022715;
		u0[CBV]      = 1;
		u0[DHG]      = 1;

		// Neuron - gating variables
		u0[m1]     	= 0.012869;
		u0[m2]     	= 0.0012175;
		u0[m3]     	= 0.1193;
		u0[m4]     	= 0.012869;
		u0[m5]     	= 0.00087377;
		u0[m6]     	= 0.0012175;
		u0[m7]     	= 0.1193;
		u0[m8]     	= 0.005;
		u0[h1]     	= 0.9718;
		u0[h2]     	= 0.12053;
		u0[h3]     	= 0.9718;
		u0[h4]     	= 0.99005;
		u0[h5]     	= 0.12053;
		u0[h6]     	= 0.9961;


}
