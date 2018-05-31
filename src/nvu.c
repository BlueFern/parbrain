// NVU 2.0


#include "nvu.h"
#include <math.h>

// nvu_init: this user-supplied function does any precomputation required
// for the model.
nvu_workspace *nvu_init(void)
{
    nvu_workspace *nvu_w;

    // Initialise the workspace
    nvu_w = malloc(sizeof *nvu_w);
    nvu_w->neq = NEQ;

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

    // Allocate other nvu workspace parameters
    nvu_w->pcap  = PCAP / P0; // pressure is nondimensional
    nvu_w->l  = 1; // normalised away
    
    // Obtain theta map from csv file and input into data_theta array for use in theta(x,y) function later
    
    
	FILE *inFile = fopen("thetamap64.csv", "r");
	
	// How many NVUs on each side of tissue slice
		int row     = pow(pow(2,NTREE-1),0.5);
		int col     = pow(pow(2,NTREE-1),0.5);

	double **data_theta;
	data_theta = (double **)malloc(row * sizeof(double *));
	for (int i = 0; i < row; ++i){
		data_theta[i] = (double *)malloc(col * sizeof(double));
	}
	
	int i_th,j_th;
	double temp;
	
	for(i_th=0;i_th<col;i_th++)
	{
		for(j_th=0;j_th<row;j_th++)
		{
		int warning = fscanf(inFile,"%lf%*[, \t\n]",&temp);  		// reads number from csv/txt file and puts it in temp, stopping at delimiters (comma, space, tab or \n)
		data_theta[i_th][j_th] = temp;					// puts number into data_theta array
//		printf("%f ",temp);
		}
//	printf("\n");
	}
	
	nvu_w->data_theta = data_theta;

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
	/*** Model fluxes/algebraic variables and state variables ***/

    double trans_p, trans_P_mmHg, delta_p_L; // pressure stuff

    // Initialise state variables
    double state_r;
    double state_v_k, state_Na_k, state_K_k, state_HCO3_k, state_Cl_k, state_Na_s, state_K_s, state_HCO3_s, state_K_p, state_w_k; // AC state
    double state_ca_i, state_ca_sr_i, state_v_i, state_w_i, state_ip3_i; // SMC state
    double state_ca_j, state_ca_er_j, state_v_j, state_ip3_j; // EC state
    double state_Mp, state_AM, state_AMp; // Mech state
    double state_ca_n, state_nNOS, state_NOn, state_NOi, state_E_b, state_E_6c, state_cGMP, state_eNOS, state_NOj, state_NOk; // NO pathway state
    double state_ca_k, state_s_k, state_h_k, state_ip3_k, state_eet_k, state_m_k, state_ca_p; // AC Ca2+ state
    double state_v_sa, state_v_d, state_K_sa, state_Na_sa, state_K_d, state_Na_d, state_K_e, state_Na_e; // NE ions state
	double state_Buff_e, state_O2, state_CBV, state_HbR; // NE other state
	double state_m1, state_m2, state_m3, state_m4, state_m5, state_m6, state_m7, state_m8, state_h1, state_h2, state_h3, state_h4, state_h5, state_h6; // NE gating state

    // Fluxes
    double Cl_s, flu_E_Na_k, flu_E_K_k, flu_E_Cl_k, flu_E_NBC_k, flu_E_BK_k, flu_J_Cl_k, flu_J_NaK_k,  flu_J_KCC1_k, flu_J_NBC_k, flu_J_NKCC1_k, flu_J_Na_k, flu_J_K_k, flu_J_BK_k, flu_w_inf, flu_phi_w; // AC fluxes
    double flu_M, flu_h_r, flu_v_cpl_i, flu_c_cpl_i, flu_I_cpl_i, flu_rho_i, flu_ip3_i, flu_SRuptake_i, flu_CICR_i, flu_extrusion_i, flu_leak_i, flu_VOCC_i, flu_NaCa_i, flu_NaK_i, flu_Cl_i, flu_K_i, flu_Kactivation_i, flu_degrad_i, flu_v_KIR_i, flu_G_KIR_i, flu_J_KIR_i, flu_J_stretch_i; // SMC fluxes
    double flu_v_cpl_j, flu_c_cpl_j, flu_I_cpl_j, flu_rho_j, flu_O_j, flu_ip3_j, flu_ERuptake_j, flu_CICR_j, flu_extrusion_j, flu_leak_j, flu_cation_j, flu_BKCa_j, flu_SKCa_j, flu_K_j, flu_R_j, flu_degrad_j, flu_J_stretch_j; // EC fluxes
    double flu_K1_c, flu_K6_c; // Mech fluxes
    double flu_c_w, flu_P_NR2AO, flu_P_NR2BO, flu_I_Ca, flu_CaM, flu_W_tau_w, flu_F_tau_w, flu_k4, flu_R_cGMP2, flu_K2_c, flu_K5_c, flu_tau_w, E_5c, V_max_pde;    // NO pathway fluxes
    double flu_p_NO_n, flu_c_NO_n, flu_d_NO_n, flu_p_NO_j, flu_c_NO_j, flu_d_NO_j;
    double rho, flu_ip3_k, flu_er_leak, flu_pump_k, flu_I_TRPV_k, flu_J_TRPV_k, flu_E_TRPV_k, B_cyt, G, v_3, H_Ca_k, eta, minf_k, t_Ca_k;  // AC Ca2+ fluxes
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
    double J_pump1_sa, J_pump1init_sa, J_pump1_d, J_pump1init_d, J_pump2, O2_p, O2_j, J_pump_sa, J_pump_d, J_Napump_sa, J_Kpump_sa, J_Napump_d, J_Kpump_d, J_pump2_0, J_pump2_O2_0;
    double J_Na_tot_sa, J_K_tot_sa, J_leak_tot_sa, J_Na_tot_d, J_K_tot_d, J_leak_tot_d, J_tot_sa, J_tot_d;
    double P_02, CBF, J_O2_vascular, J_O2_background, J_O2_pump;
    double f_out, CMRO2, CMRO2_init;
    double Glu;


    // State Variables:
    state_r  	  = u[i_radius];

    state_ca_i    = u[i_ca_i];
    state_ca_sr_i = u[i_ca_sr_i];
    state_v_i     = u[i_v_i];
    state_w_i     = u[i_w_i];
    state_ip3_i   = u[i_ip3_i];

    state_v_k     = u[i_v_k];
    state_Na_k    = u[i_Na_k];
    state_K_k     = u[i_K_k];
    state_HCO3_k  = u[i_HCO3_k];
    state_Cl_k    = u[i_Cl_k];
    state_Na_s    = u[i_Na_s];
    state_K_s     = u[i_K_s];
    state_HCO3_s  = u[i_HCO3_s];
    state_K_p     = u[i_K_p];
    state_w_k     = u[i_w_k];

    state_ca_j    = u[i_ca_j];
    state_ca_er_j = u[i_ca_er_j];
    state_v_j     = u[i_v_j];
    state_ip3_j   = u[i_ip3_j];

    state_Mp      = u[i_Mp];
    state_AMp     = u[i_AMp];
    state_AM      = u[i_AM];

    state_ca_n    = u[i_ca_n];
    state_nNOS    = u[i_nNOS];
    state_NOn     = u[i_NOn];
    state_NOi     = u[i_NOi];
    state_E_b     = u[i_E_b];
    state_E_6c    = u[i_E_6c];
    state_cGMP    = u[i_cGMP];
    state_eNOS    = u[i_eNOS];
    state_NOj     = u[i_NOj];
    state_NOk     = u[i_NOk];

    state_ca_k     = u[i_ca_k];
    state_s_k      = u[i_s_k];
    state_h_k      = u[i_h_k];
    state_ip3_k    = u[i_ip3_k];
    state_eet_k    = u[i_eet_k];
    state_m_k      = u[i_m_k];
    state_ca_p     = u[i_ca_p];

    state_v_sa		= u[i_v_sa];
    state_v_d		= u[i_v_d];
	state_K_sa		= u[i_K_sa];
	state_Na_sa		= u[i_Na_sa];
	state_K_d		= u[i_K_d];
	state_Na_d		= u[i_Na_d];
	state_K_e		= u[i_K_e];
	state_Na_e		= u[i_Na_e];
	state_Buff_e	= u[i_Buff_e];
	state_O2		= u[i_O2];
	state_CBV		= u[i_CBV];
	state_HbR		= u[i_HbR];
	state_m1		= u[i_m1];
	state_m2		= u[i_m2];
	state_m3		= u[i_m3];
	state_m4		= u[i_m4];
	state_m5		= u[i_m5];
	state_m6		= u[i_m6];
	state_m7		= u[i_m7];
	state_m8		= u[i_m8];
	state_h1		= u[i_h1];
	state_h2		= u[i_h2];
	state_h3		= u[i_h3];
	state_h4		= u[i_h4];
	state_h5		= u[i_h5];
	state_h6		= u[i_h6];

// Fluxes:

	// Dimensional radius
	state_R_dim			= state_r * R_0_passive_k;

    // pressure
    trans_p 			= P0 / 2 * (p + nvu_w->pcap); 				// x P0 to make dimensional, transmural pressure in Pa.  4000 in matlab code
	trans_P_mmHg 		= trans_p * PA2MMHG; 						// transmural pressure in mmHg.   30 in matlab code
	delta_p_L			= P0 * (p - nvu_w->pcap) / L0; 			// dimensional pressure drop over leaf vessel divided by length of the vessel (200um). 18.2Pa / 200e-6m = 91000 in matlab code

    // SC fluxes
    Cl_s         		= state_Na_s + state_K_s - state_HCO3_s;  //

    // AC fluxes
	flu_E_TRPV_k		= ph / z_Ca * log(state_ca_p / state_ca_k); // TRPV4 channel Nernst Potential
    flu_E_Na_k         	= ph / z_Na * log(state_Na_s / state_Na_k);    // mV
    flu_E_K_k          	= ph / z_K * log(state_K_s / state_K_k );
    flu_E_Cl_k         	= ph / z_Cl * log(Cl_s / state_Cl_k);
    flu_E_NBC_k        	= ph / z_NBC * log((state_Na_s * pow(state_HCO3_s,2))/(state_Na_k * pow(state_HCO3_k,2)));
    flu_E_BK_k         	= ph / z_K * log(state_K_p / state_K_k);

    flu_J_NaK_k        	= J_NaK_max * ( pow(state_Na_k,1.5) / ( pow(state_Na_k,1.5) + pow(K_Na_k,1.5) ) ) * ( state_K_s / (state_K_s + K_K_s) );    // uM s-1
    flu_J_KCC1_k       	= ph * G_KCC1_k * log((state_K_s * Cl_s)/(state_K_k*state_Cl_k)) ;   //uM s-1
    flu_J_NBC_k        	= G_NBC_k * (state_v_k - flu_E_NBC_k);       //uM s-1
    flu_J_NKCC1_k     	= ph * G_NKCC1_k * log((state_K_s * state_Na_s * pow(Cl_s,2)) / (state_K_k * state_Na_k * pow(state_Cl_k,2)) );        //uM s-1
    flu_J_Na_k   		= G_Na_k * (state_v_k - flu_E_Na_k);              //uM s-1
    flu_J_K_k    		= G_K_k * (state_v_k - flu_E_K_k);          //uM s-1
    flu_J_BK_k   		= G_BK_k * state_w_k * (state_v_k - flu_E_BK_k);  //uM s-1
    flu_J_Cl_k			= G_Cl_k * (state_v_k - flu_E_Cl_k);

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
    flu_v_KIR_i    		= z_1 * state_K_p - z_2;
    flu_G_KIR_i    		= exp( z_5 * state_v_i + z_3 * state_K_p - z_4 );
    flu_J_KIR_i    		= F_il/gam * flu_G_KIR_i * (state_v_i - flu_v_KIR_i);

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
    flu_BKCa_j 			= 0.2 * ( 1 + tanh( ( (  log10(state_ca_j) - c_j) * ( state_v_j - b_j ) - a1_j ) / ( m3b* pow(( state_v_j + a2_j * ( log10( state_ca_j ) - c_j ) - b_j),2) + m4b ) ) );
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

	Glu 				= GluSwitch * ( 0.5 * Glu_max * ( 1 + tanh( (state_K_e - Ke_switch) / Glu_slope) ) );

	flu_P_NR2AO         = Glu / (betA + Glu);
	flu_P_NR2BO         = Glu / (betB + Glu);
	flu_I_Ca            = (-4 * v_n * G_M * P_Ca_P_M * (Ca_ex / M_mono)) / (1 + exp(-80 * (v_n + 0.02))) * (exp(2 * v_n * F / (R_gas * Temp))) / (1 - exp(2 * v_n * F / (R_gas * Temp))) * (0.63 * flu_P_NR2AO + 11 * flu_P_NR2BO);
	flu_CaM             = state_ca_n / m_c;                                      // concentration of calmodulin / calcium complexes ; (100)

	flu_p_NO_n			= NOswitch * ( state_nNOS * V_max_NO_n * On / (K_mO2_n + On) * LArg / (K_mArg_n + LArg) );
    flu_c_NO_n			= (k_O2* pow(state_NOn,2) * On);
	flu_d_NO_n			= - ((state_NOn - state_NOk) / tau_nk);

	flu_tau_w			= state_R_dim/2 * delta_p_L; // WSS using pressure from the H tree. L_0 = 200 um
	flu_W_tau_w         = W_0 * pow((flu_tau_w + sqrt(16 * pow(delt_wss,2) + pow(flu_tau_w,2)) - 4 * delt_wss),2) / (flu_tau_w + sqrt(16 * pow(delt_wss,2) + pow(flu_tau_w,2))) ;  // - tick
	flu_F_tau_w         = 1 / (1 + alp * exp(-flu_W_tau_w)) - 1 / (1 + alp); // -(1/(1+alp)) was added to get no NO at 0 wss (!) - tick
	flu_k4              = C_4 * pow(state_cGMP, m_4);
	flu_R_cGMP2         = (pow(state_cGMP, 2)) / (pow(state_cGMP, 2) + pow(K_m_mlcp, 2)); // - tick
	flu_K2_c            = 58.1395 * k_mlcp_b + 58.1395 * k_mlcp_c * flu_R_cGMP2;  // Factor is chosen to relate two-state model of Yang2005 to Hai&Murphy model
	flu_K5_c            = flu_K2_c;
	flu_c_w             = 1/2 * ( 1 + tanh( (state_cGMP - 10.75) / 0.668 ) );

	O2_j 				= state_O2 * 1e3;
	flu_p_NO_j			= NOswitch * ( V_NOj_max * state_eNOS * O2_j / (K_mO2_j + O2_j) * LArg_j / (K_mArg_j + LArg_j) );
	flu_c_NO_j			= k_O2 * pow(state_NOj,2) * O2_j;
	flu_d_NO_j			= (state_NOi - state_NOj) / tau_ij - state_NOj * 4 * D_NO / (pow(25,2));

	flu_Kactivation_i   = pow((state_ca_i + flu_c_w),2) / (pow((state_ca_i + flu_c_w),2) + bet_i * exp(-(state_v_i - v_Ca3) / R_K));
	E_5c				= 1 - state_E_b - state_E_6c;
	V_max_pde			= k_pde * state_cGMP;

// AC Ca2+ fluxes

	rho				= rho_min + (rho_max - rho_min)/Glu_max * Glu;

    flu_ip3_k 		= J_max * pow(( state_ip3_k / ( state_ip3_k + K_I ) * state_ca_k / ( state_ca_k + K_act ) * state_h_k ) , 3) * (1.0 - state_ca_k / state_s_k);
	flu_er_leak 	= P_L * ( 1.0 - state_ca_k / state_s_k );
	flu_pump_k  	= V_max * pow(state_ca_k, 2) / ( pow(state_ca_k, 2) + pow(k_pump, 2) );
	flu_I_TRPV_k	= G_TRPV_k * state_m_k * (state_v_k - flu_E_TRPV_k);
	flu_J_TRPV_k	= - flu_I_TRPV_k / ( z_Ca * C_astr_k * gam );
	B_cyt 			= 1.0 / (1.0 + BK_end + K_ex * B_ex / pow((K_ex + state_ca_k), 2) );
	G				= ( rho + delta ) / ( K_G + rho + delta );
	v_3				= -v_5 / 2.0 * tanh((state_ca_k - Ca_3) / Ca_4) + v_6;
	flu_w_inf    	= 0.5 * ( 1 + tanh( ( state_v_k + eet_shift * state_eet_k - v_3 ) / v_4 ) );
    flu_phi_w    	= psi_w * cosh( (state_v_k - v_3) / (2*v_4) );
    H_Ca_k			= state_ca_k / gam_cai_k + state_ca_p / gam_cae_k;
    eta 			= (state_R_dim - R_0_passive_k) / R_0_passive_k;
    minf_k 			= ( 1 / ( 1 + exp( - (eta - epshalf_k) / kappa_k ) ) ) * ( ( 1 / (1 + H_Ca_k) ) * (H_Ca_k + tanh(( state_v_k - v1_TRPV_k) / v2_TRPV_k )));
    t_Ca_k 			= t_TRPV_k / state_ca_p;

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
	J_NaP_sa 		= (pow(state_m1,2) * state_h1 * gNaP_GHk * Farad * state_v_sa * (state_Na_sa - (exp(-state_v_sa / ph) * state_Na_e))) / (ph * (1 - exp(-state_v_sa / ph)));

	// Na+ flux through NaT in soma/axon
	m8alpha 		= 0.32 * ((-state_v_sa - 51.9) / (exp(-(0.25 * state_v_sa + 12.975)) - 1));
	m8beta 			= 0.28 * ((state_v_sa + 24.89) / (exp(0.2 * state_v_sa + 4.978) - 1));
	h6alpha 		= 0.128 * exp(-(0.056 * state_v_sa + 2.94));
	h6beta 			= 4 / (1 + exp(-(0.2 * state_v_sa + 6)));
	J_NaT_sa 		= (pow(state_m8,3) * state_h6 * gNaT_GHk * Farad * state_v_sa * (state_Na_sa - (exp(-state_v_sa / ph) * state_Na_e))) / (ph * (1 - exp(-state_v_sa / ph)));

	// K+ flux through KDR in soma/axon
	m2alpha 		= 0.016 * ((state_v_sa + 34.9) / (1 - exp(-((0.2 * state_v_sa) + 6.98))));
	m2beta 			= 0.25 * exp(-((0.025 * state_v_sa) + 1.25));
	J_KDR_sa 		= (pow(state_m2,2) * gKDR_GHk * Farad * state_v_sa * (state_K_sa - (exp(-state_v_sa / ph) * state_K_e))) / (ph * (1 - exp(-state_v_sa / ph)));

	// K+ flux through KA in soma/axon
	m3alpha 		= 0.02 * ((state_v_sa + 56.9) / (1 - exp(-((0.1 * state_v_sa) + 5.69))));
	m3beta 			= 0.0175 * ((state_v_sa + 29.9) / (exp(((0.1 * state_v_sa) + 2.99)) - 1));
	h2alpha 		= 0.016 * exp(-((0.056 * state_v_sa) + 4.61));
	h2beta 			= 0.5 / (1 + exp(-((0.2 * state_v_sa) + 11.98)));
	J_KA_sa 		= (pow(state_m3,2) * state_h2 * gKA_GHk * Farad * state_v_sa * (state_K_sa - (exp(-state_v_sa / ph) * state_K_e))) / (ph * (1 - exp(-state_v_sa / ph)));

	// Na+ flux through NaP in dendrite
	m4alpha 		= 1 / (6 * (1 + exp(-((0.143 * state_v_d) + 5.67))));
	m4beta 			= exp(-((0.143 * state_v_d) + 5.67)) / (6 * (1 + exp(-((0.143 * state_v_d) + 5.67))));
	h3alpha 		= 5.12e-8 * exp(-((0.056 * state_v_d) + 2.94));
	h3beta 			= 1.6e-6 / (1 + exp(-(((0.2 * state_v_d)) + 8)));
	J_NaP_d 		= (pow(state_m4,2) * state_h3 * gNaP_GHk * Farad * state_v_d * (state_Na_d - (exp(-state_v_d / ph) * state_Na_e))) / (ph * (1 - exp(-state_v_d / ph)));

	// Na+ and K+ flux through NMDA in dendrite
	m5alpha 		= 0.5 / (1 + exp((13.5 - state_K_e) / 1.42));
	m5beta 			= 0.5 - m5alpha;
	h4alpha 		= 1 / (2000 * (1 + exp((state_K_e - 6.75) / 0.71)));
	h4beta 			= 5e-4 - h4alpha;
	J_NMDA_K_d 		= ( (state_m5 * state_h4 * gNMDA_GHk * Farad * state_v_d * (state_K_d - (exp(-state_v_d / ph) * state_K_e))) / (ph * (1 - exp(-state_v_d / ph))) ) / (1 + 0.33 * Mg * exp(-(0.07 * state_v_d + 0.7)));
	J_NMDA_Na_d 	= ( (state_m5 * state_h4 * gNMDA_GHk * Farad * state_v_d * (state_Na_d - (exp(-state_v_d / ph) * state_Na_e))) / (ph * (1 - exp(-state_v_d / ph))) ) / (1 + 0.33 * Mg * exp(-(0.07 * state_v_d + 0.7)));

	// K+ flux through KDR in dendrite
	m6alpha 		= 0.016 * ((state_v_d + 34.9) / (1 - exp(-((0.2 * state_v_d) + 6.98))));
	m6beta 			= 0.25 * exp(-((0.025 * state_v_d) + 1.25));
	J_KDR_d 		= (pow(state_m6,2) * gKDR_GHk * Farad * state_v_d * (state_K_d - (exp(-state_v_d / ph) * state_K_e))) / (ph * (1 - exp(-state_v_d / ph)));

	// K+ flux through KA in dendrite
	m7alpha 		= 0.02 * ((state_v_d + 56.9) / (1 - exp(-((0.1 * state_v_d) + 5.69))));
	m7beta 			= 0.0175 * ((state_v_d + 29.9) / (exp(((0.1 * state_v_d) + 2.99)) - 1));
	h5alpha 		= 0.016 * exp(-((0.056 * state_v_d) + 4.61));
	h5beta 			= 0.5 / (1 + exp(-((0.2 * state_v_d) + 11.98)));
	J_KA_d 			= (pow(state_m7,2) * state_h5 * gKA_GHk * Farad * state_v_d * (state_K_d - (exp(-state_v_d / ph) * state_K_e))) / (ph * (1 - exp(-state_v_d / ph)));

	// ATPase pump
	J_pump1_sa 		= pow((1 + (K_init_e / state_K_e)),-2) * pow((1 + (Na_init_sa / state_Na_sa)),-3);
	J_pump1init_sa 	= 0.0312;				// (1 + (K_init_e / K_init_e))^(-2) * (1 + (Na_init_sa / Na_init_sa)) ^ (-3);
	J_pump1_d 		= pow((1 + (K_init_e / state_K_e)),-2) * pow((1 + (Na_init_d / state_Na_d)),-3);
	J_pump1init_d 	= 0.0312;				// (1 + (K_init_e / K_init_e))^(-2) * (1 + (Na_init_d / Na_init_d)) ^ (-3);

	// Jpump2 uses 02_0 if oxygen is plentiful or O2 if limited. Default is limited, change using O2switch in constants.h
	O2_p        	= O2_0 * (1 - O2switch) + state_O2 * O2switch;
	J_pump2 		= 2 * pow((1 + O2_0 / (((1 - alpha_O2) * O2_p) + alpha_O2 * O2_0)),-1);

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
	f_out 			= pow(state_CBV,(1/d_BOLD)) + tau_TAT * (1/(tau_MTT + tau_TAT) * ( CBF/CBF_init  - pow(state_CBV,(1/d_BOLD)) ));
	CMRO2 			= J_O2_background + J_O2_pump;
	CMRO2_init 		= CBF_init * P_02;
	//OEF 			= CMRO2 * E_0 / CBF;

	// Not actually used here - used in bin_to_vtu
	//BOLD 			= 100 * V_0 * ( a_1 * (1 - state_HbR/HbR_0) - a_2 * (1 - state_CBV/CBV_0) ); // divides HbR, CBV and CBF by initial values, will depend on J_PLC
	//CBF_change	= (CBF - CBF_0)/CBF_0;
	//CBF_N			= CBF/CBF_0;
	//HBT_N			= CBF_N * HbR_N / CMRO2_N;
	//HBO_N			= (HBT_N - 1) - (HbR/HbR_0 - 1) + 1;


// Differential Equations:

    /***********Neuron dynamics*********/

	// Neuron - other
    du[i_Buff_e]   	= Mu * state_K_e * (Buff0 - state_Buff_e) / (1 + exp(-((state_K_e - 5.5) / 1.09))) - (Mu * state_Buff_e); // [mM]
    du[i_O2]       	= J_O2_vascular - J_O2_background - J_O2_pump; 						// [mM]
    du[i_CBV]      	= (1/(tau_MTT + tau_TAT)) * ( CBF/CBF_init  - pow(state_CBV,(1/d_BOLD)) );	// [-]
    du[i_HbR]      	= 1/tau_MTT * ( CMRO2/CMRO2_init - state_HbR/state_CBV * f_out );		// [-]

    // Neuron - ions (mV or mM)
    du[i_v_sa]	 	= 1/Cm * ( -J_tot_sa + 1 / (2 * Ra * pow(dhod,2)) * (state_v_d - state_v_sa) + current_input(t,x,y) );
	du[i_v_d]      	= 1/Cm * (-J_tot_d + 1 / (2 * Ra * pow(dhod,2)) * (state_v_sa - state_v_d));
    du[i_K_sa]     	= -As / (Farad * Vs) * J_K_tot_sa + D_K * (Vd + Vs) / (2 * pow(dhod,2) * Vs) * (state_K_d - state_K_sa);
    du[i_Na_sa]    	= -As / (Farad * Vs) * J_Na_tot_sa + D_Na * (Vd + Vs) / (2 * pow(dhod,2) * Vs) * (state_Na_d - state_Na_sa);
    du[i_K_d]      	= -Ad / (Farad * Vd) * J_K_tot_d + D_K * (Vs + Vd) / (2 * pow(dhod,2) * Vd) * (state_K_sa - state_K_d);
    du[i_Na_d]     	= -Ad / (Farad * Vd) * J_Na_tot_d + D_Na * (Vs + Vd) / (2 * pow(dhod,2) * Vd) * (state_Na_sa - state_Na_d);
    du[i_K_e]      	= 1/(Farad * fe) * (((As * J_K_tot_sa) / Vs)  + ((Ad * J_K_tot_d) / Vd)) - du[i_Buff_e];
    du[i_Na_e]     	= 1/(Farad * fe) * (((As * J_Na_tot_sa) / Vs) + ((Ad * J_Na_tot_d) / Vd));

	// Neuron - gating variables ([-])
    du[i_m1]     	= 1000 * ((m1alpha * (1 - state_m1)) - (m1beta * state_m1));
    du[i_m2]     	= 1000 * ((m2alpha * (1 - state_m2)) - (m2beta * state_m2));
    du[i_m3]     	= 1000 * ((m3alpha * (1 - state_m3)) - (m3beta * state_m3));
    du[i_m4]     	= 1000 * ((m4alpha * (1 - state_m4)) - (m4beta * state_m4));
    du[i_m5]     	= 1000 * ((m5alpha * (1 - state_m5)) - (m5beta * state_m5));
    du[i_m6]     	= 1000 * ((m6alpha * (1 - state_m6)) - (m6beta * state_m6));
    du[i_m7]     	= 1000 * ((m7alpha * (1 - state_m7)) - (m7beta * state_m7));
    du[i_m8]     	= 1000 * ((m8alpha * (1 - state_m8)) - (m8beta * state_m8));
    du[i_h1]     	= 1000 * ((h1alpha * (1 - state_h1)) - (h1beta * state_h1));
    du[i_h2]     	= 1000 * ((h2alpha * (1 - state_h2)) - (h2beta * state_h2));
    du[i_h3]     	= 1000 * ((h3alpha * (1 - state_h3)) - (h3beta * state_h3));
    du[i_h4]     	= 1000 * ((h4alpha * (1 - state_h4)) - (h4beta * state_h4));
    du[i_h5]     	= 1000 * ((h5alpha * (1 - state_h5)) - (h5beta * state_h5));
    du[i_h6]     	= 1000 * ((h6alpha * (1 - state_h6)) - (h6beta * state_h6));

	/***********General dynamics**********/

    du[i_radius] 	= 1 / ETA * (state_r * trans_p / flu_h_r - E * (state_R_dim - R_0)/R_0); // Radius - nondimensional (state_R_dim: dimensional)

    //AC:
    du[i_v_k] 		= gam * ( -flu_J_BK_k - flu_J_K_k - flu_J_Cl_k - flu_J_NBC_k - flu_J_Na_k - flu_J_NaK_k + 2*flu_J_TRPV_k);
    du[i_Na_k] 		= -flu_J_Na_k - 3 * flu_J_NaK_k + flu_J_NKCC1_k + flu_J_NBC_k;
    du[i_K_k] 		= -flu_J_K_k + 2 * flu_J_NaK_k + flu_J_NKCC1_k + flu_J_KCC1_k - flu_J_BK_k;
    du[i_HCO3_k] 	= 2 * flu_J_NBC_k;
    du[i_Cl_k] 		= du[ i_Na_k] + du[ i_K_k] - du[ i_HCO3_k];
    du[i_w_k]		= flu_phi_w * (flu_w_inf - state_w_k);

    //SC:
    du[i_Na_s] 		= 1/VR_sa * ( - du[ i_Na_k] ) - SC_coup * du[i_K_e] * 1e3;
    du[i_K_s] 		= 1/VR_sa * ( flu_J_K_k - 2 * flu_J_NaK_k - flu_J_NKCC1_k - flu_J_KCC1_k ) + SC_coup * du[i_K_e] * 1e3;
    du[i_HCO3_s] 	= 1/VR_sa * ( - du[ i_HCO3_k] );

    //PVS:
    du[i_K_p] 		= flu_J_BK_k / VR_pa + flu_J_KIR_i / VR_ps - R_decay * (state_K_p - K_p_min);

    //SMC:
    du[i_ca_i] 		= flu_c_cpl_i + flu_rho_i * ( flu_ip3_i - flu_SRuptake_i + flu_CICR_i - flu_extrusion_i + flu_leak_i - flu_VOCC_i + flu_NaCa_i - 0.1* flu_J_stretch_i);
    du[i_ca_sr_i] 	= flu_SRuptake_i - flu_CICR_i - flu_leak_i ;
    du[i_v_i] 		= flu_v_cpl_i + gam * ( - flu_NaK_i - flu_Cl_i - 2 * flu_VOCC_i - flu_NaCa_i - flu_K_i - flu_J_stretch_i - flu_J_KIR_i );
    du[i_w_i] 		= lam * (flu_Kactivation_i - state_w_i ) ;
    du[i_ip3_i] 	= flu_I_cpl_i - flu_degrad_i ;
    du[i_K_i]		= - flu_J_KIR_i - flu_K_i + flu_NaK_i;

    //EC:
    du[i_ca_j] 		= flu_c_cpl_j + flu_rho_j * ( flu_ip3_j - flu_ERuptake_j + flu_CICR_j - flu_extrusion_j + flu_leak_j + flu_cation_j + flu_O_j - flu_J_stretch_j ) ;
    du[i_ca_er_j] 	= flu_ERuptake_j - flu_CICR_j - flu_leak_j ;
    du[i_v_j] 		= flu_v_cpl_j - 1/C_m * ( flu_K_j + flu_R_j ) ;
    du[i_ip3_j] 	= flu_I_cpl_j + J_PLC - flu_degrad_j ;  // **

    // Mech:
    du[i_Mp] 		= wallMech * ( K4_c * state_AMp + flu_K1_c * flu_M - (flu_K2_c + K3_c) * state_Mp );
    du[i_AMp] 		= wallMech * ( K3_c * state_Mp + flu_K6_c * state_AM - (K4_c + flu_K5_c) * state_AMp );
    du[i_AM] 		= wallMech * ( flu_K5_c * state_AMp - ( K7_c + flu_K6_c ) * state_AM );

    /***********NO pathway***********/

    // NE:
    du[i_ca_n]      = (flu_I_Ca / (2*F * V_spine) - (k_ex * (state_ca_n - Ca_rest))) / (1 + lambda);
    du[i_nNOS]      = V_maxNOS * flu_CaM / (K_actNOS + flu_CaM) - mu2 * state_nNOS ;
    du[i_NOn]       = flu_p_NO_n - flu_c_NO_n + flu_d_NO_n;

    // AC:
    du[i_NOk]       = (state_NOn - state_NOk) / tau_nk + (state_NOi - state_NOk) / tau_ki - k_O2 * pow(state_NOk,2) * Ok;

    // SMC:
    du[i_NOi]       = (state_NOk - state_NOi) / tau_ki + (state_NOj - state_NOi) / tau_ij - k_dno * state_NOi ;
    du[i_E_b]       = -k1_i * state_E_b * state_NOi + k_1_i * state_E_6c + flu_k4 * E_5c;
    du[i_E_6c]      = k1_i * state_E_b * state_NOi - k_1_i * state_E_6c - k2_i * state_E_6c - k3_i * state_E_6c * state_NOi ;
    du[i_cGMP]      = V_max_sGC * E_5c - V_max_pde * state_cGMP / (K_m_pde + state_cGMP);

    // EC:
    du[i_eNOS]      = gam_eNOS * (K_dis * state_ca_j / (K_eNOS + state_ca_j)) + (1 - gam_eNOS) * (g_max * flu_F_tau_w) - mu2 * state_eNOS;
    du[i_NOj]       = flu_p_NO_j - flu_c_NO_j + flu_d_NO_j;

    /**********Astrocytic Calcium*******/

    // AC:
    du[i_ca_k]		= B_cyt * (flu_ip3_k - flu_pump_k + flu_er_leak + flu_J_TRPV_k/r_buff);
    du[i_s_k]		= -(B_cyt * (flu_ip3_k - flu_pump_k + flu_er_leak)) / (VR_ER_cyt);
    du[i_h_k]		= k_on * (K_inh - (state_ca_k + K_inh) * state_h_k);
    du[i_ip3_k]		= r_h * G - k_deg * state_ip3_k;
    du[i_eet_k]		= V_eet * fmax((state_ca_k - Ca_k_min), 0) - k_eet * state_eet_k;
    du[i_m_k]		= trpv_switch * (minf_k - state_m_k) / (t_Ca_k * state_ca_p);

	//PVS:
    du[i_ca_p]		= -flu_J_TRPV_k / VR_pa + flu_VOCC_i / VR_ps - Ca_decay_k * (state_ca_p - Capmin_k);

    // Curvature variables - don't change with time but added as ODEs for convenience
    du[i_curvature] = 0;
    du[i_coup] = 0;

}

void read_csv(int row, int col, char *filename, double **data){
	FILE *file;
	file = fopen(filename, "r");

	int i = 0;
    char line[4098];
	while (fgets(line, 4098, file) && (i < row))
    {
    	// double row[ssParams->nreal + 1];
        char* tmp = strdup(line);

	    int j = 0;
	    const char* tok;
	    for (tok = strtok(line, "\t"); tok && *tok; j++, tok = strtok(NULL, "\t\n"))
	    {
	        data[i][j] = atof(tok);
	        printf("%f\t", data[i][j]);
	    }
	    printf("\n");

        free(tmp);
        i++;
    }
}

// Theta (curvature) spatially varied ***change
double theta_function(double x, double y, nvu_workspace *nvu_w)
{
	double** data_theta = nvu_w->data_theta; // data from csv mapping
	
	// Number of NVUs along the side of the tissue slice
	int row     = pow(pow(2,NTREE-1),0.5);
	int col     = pow(pow(2,NTREE-1),0.5);
	
	int i,j; // Coordinates for the tissue slice
	
	// Convert x,y into i,j, add small number 1e-9 to somehow fix issue in converting double to int
	i = (x + 0.0002*(row-1))/0.0004 + 1e-9;
	j = (y + 0.0002*(col-1))/0.0004 + 1e-9;

    double theta = data_theta[i][j];
//    printf("x:%f y%f i:%d j:%d theta:%f\n", x, y, i, j, theta);
    return theta;
}


// Time-varying pressure at the root of the tree. 1 is nominal value. If
// you want to work in unscaled units, make sure you *multiply* by P0
// afterwards
double nvu_p0(double t)
{
    //double p0 = 1. * 8000 / P0; // 8000 Pa   original: 1.5 * 8000 / P0;
    //double p0 = (0.5 * sin(t) + 1) * 8000 / P0; //

    double p0 = P_TOP / P0;
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
    double I_strength 	= I_STRENGTH;
    double t_up   		= T_STIM_0;
    double t_down 		= T_STIM_END;

    double ampl = 2;
    double ramp = 0.00075; //0.0035 for N=13;
    double x_centre = 0;
    double y_centre = 0;
    double current_space;

    if (SPATIAL_CHOICE)
	{
    current_space = fmin(1.0, ampl * (exp(-((pow((x - x_centre), 2) + pow((y - y_centre), 2)) / (2 * pow(ramp, 2))))));
	}
    else
    {
    // in centre square with 'radius'=6 blocks
        if ( fmax(fabs(x),fabs(y)) <= (6*0.0004 - 0.0002))
	//    if (x <= 0 && y <= 0)
		{
			current_space = 1;
		}
		else
		{
			current_space = 0;
		}
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
		// Current ICs optimised for J_PLC = 0.11!!
		u0[i_radius]  = 1.1485;

	    u0[i_v_k]     = -88.9;
	    u0[i_Na_k]    = 18268;
	    u0[i_K_k]     = 92708;
	    u0[i_HCO3_k]  = 9131;
	    u0[i_Cl_k]    = 7733;
	    u0[i_w_k]     = 1.703e-4;

	    u0[i_Na_s]    = 150255;
	    u0[i_K_s]     = 2837;
	    u0[i_HCO3_s]  = 16881;

	    u0[i_K_p]       = 3045.1;

	    u0[i_ca_i]      = 0.2637;
	    u0[i_ca_sr_i]   = 1.1686;
	    u0[i_v_i]       = -34.7;
	    u0[i_w_i]       = 0.2206;
	    u0[i_ip3_i]     = 0.275;
	    u0[i_K_i]       = 99994.8;

	    u0[i_ca_j]      = 0.8331;
	    u0[i_ca_er_j]   = 0.6266;
	    u0[i_v_j]       = -68.27;
	    u0[i_ip3_j]     = 0.825;

	    u0[i_Mp]        = 0.0842;
	    u0[i_AMp]       = 0.0622;
	    u0[i_AM]        = 0.2746;

	    // NO pathway
		u0[i_NOn]       = 0.1671;
		u0[i_NOk]       = 0.1106;
		u0[i_NOi]       = 0.0541;
		u0[i_NOj]       = 0.0528;
		u0[i_cGMP]      = 8.2826;
		u0[i_eNOS]      = 0.4479;
		u0[i_nNOS]      = 0.318;
		u0[i_ca_n]      = 0.1;
		u0[i_E_b]       = 0.4077;
		u0[i_E_6c]      = 0.4396;

		// AC Ca2+ pathway
		u0[i_ca_k]     = 0.1612;
		u0[i_s_k]      = 480.8;
		u0[i_h_k]      = 0.3828;
		u0[i_ip3_k]    = 0.048299;
		u0[i_eet_k]    = 0.6123;
		u0[i_m_k]      = 0.5710;
		u0[i_ca_p]     = 1746.4;

		// Neuron - ions
		u0[i_v_sa]     = -70.0337;
		u0[i_v_d]      = -70.0195;
		u0[i_K_sa]     = 134.1858;
		u0[i_Na_sa]    = 9.2691;
		u0[i_K_d]      = 134.4198;
		u0[i_Na_d]     = 9.3203;
		u0[i_K_e]      = 3.493;
		u0[i_Na_e]     = 150;

		// Neuron - other
		u0[i_Buff_e]   = 165.9812;
		u0[i_O2]       = 0.0281;
		u0[i_CBV]      = 1.31557;
		u0[i_HbR]      = 0.667547;

		// Neuron - gating variables
		u0[i_m1]     	= 0.01281;
		u0[i_m2]     	= 0.001209;
		u0[i_m3]     	= 0.1190;
		u0[i_m4]     	= 0.01284;
		u0[i_m5]     	= 0.000869;
		u0[i_m6]     	= 0.001213;
		u0[i_m7]     	= 0.1191;
		u0[i_m8]     	= 0.004962;
		u0[i_h1]     	= 0.9718;
		u0[i_h2]     	= 0.1214;
		u0[i_h3]     	= 0.9718;
		u0[i_h4]     	= 0.9899;
		u0[i_h5]     	= 0.1210;
		u0[i_h6]     	= 0.9961;

		u0[i_curvature]	= cos(theta_function(x,y,nvu_w))/(pow(r_th,2) * (n_th + cos(theta_function(x,y,nvu_w)))); // Gaussian curvature over x,y coordinates
		u0[i_coup]		= 10/(pow(a_th,2)) * pow(( cosh(eta_th) - n_th + pow(a_th,2)* (cos(theta_function(x,y,nvu_w))/(pow(r_th,2) * (n_th + cos(theta_function(x,y,nvu_w))))) / cos(theta_function(x,y,nvu_w)) ) , 2); // Diffusion scaling rate

//		printf("x:%f y:%f theta:%f curvature:%f C:%f\n",x,y,theta_function(x,y,nvu_w),u0[i_curvature],u0[i_coup]);

}
