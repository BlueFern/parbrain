#ifndef SRC_CONSTANTS_H_
#define SRC_CONSTANTS_H_

// Optional command line arguments for parBrainSim: N, T_FINAL, DT_PSEC (in that order). If none specified then the following are used.

/*** Run parameters ***/
    static const double T_FINAL        	= 5;        // Final run time
	static const double T_STIM_0       	= 2;        // Start time for stimulation
	static const double T_STIM_END     	= 6;        // End time for stimulation
    static const int    DT_PSEC       	= 10;       // Time step for writing to file (and screen)
    static const int 	NTREE          	= 13;         // Number of levels in the H-tree (where the tissue slice has 2^(N-1) tissue blocks)
    static const int 	NSUB           	= 1;         // Subtree size (easiest to just keep as 1)
	static const double P_TOP			= 4175;	     // Pressure at the top of the tree, chosen so that the drop over the terminating arterioles is around 18.2 Pa to match with the single NVU model.
									  	  	  	  	 // For NTREE=3, P_TOP=4100 Pa. For NTREE=7, P_TOP=4160 Pa. For NTREE=13, P_TOP=4175
	static const int 	SPATIAL_CHOICE	= 0;	     // 1: current input is a Gaussian plateau into the centre (fixed size), 0: square input

/*** Switches for various pathways ***/
	static const double DIFFUSION_SWITCH 	= 2;		// 2: ECS electrodiffusion, 1: extracellular diffusion between blocks, 0: none
	static const double GJ_SWITCH 			= 0;		// 2: multiple ion astrocytic gap junctions, 1: just K+ astrocytic gap junctions, 0: none
	static const double GluSwitch			= 1;		// 1: glutamate is released with current stimulation, 0: no glutamate
	static const double NOswitch			= 1;		// 1: NO is produced in the NVU, 0: no NO production at all
    static const double trpv_switch	    	= 1;		// 1: TRPV4 channel is active, 0: completely closed (no flux)
    static const double O2switch			= 1;		// 1: Oxygen is limited, 0: oxygen is plentiful (default 1)

/*** Curvature constants ***/
    static const double r_th = 3.1831;	// Minor radius of torus
    static const double n_th = 4;		// Major:Minor radius ratio - Decreasing n (but must have n>1) increases how curvy the surface will be (higher max and min), also modifies a_th and eta_th!!!
    static const double a_th = 12.3281;	// a = r*sqrt(n^2-1)
	static const double eta_th = 2.0634;// eta = atanh(sqrt(n^2-1)/n)
    
/*** Commonly changed model parameters ***/

    /* Normal conditions (neurovascular coupling - stimulation then vasodilation) */
//	static const double I_STRENGTH		= 0.022;	  	// [A] strength of current input (default 0.022)
//    static const double SC_coup	        = 11.5;         // scaling factor for the change in SC K+ concentration based on extracellular K+ concentration (default 11.5)
//    static const double Imax		    = 0.013*6;         // rate of the ATP pump (default 0.013*6)
//    static const double gNaleak_sa	    = 6.2378e-5;       // channel conductances, change depending on Imax, see OO-NVU for other values
//    static const double gKleak_sa	    = 2.1989e-4;
//    static const double gleak_sa	    = 10*6.2378e-5;
//    static const double gNaleak_d	    = 6.2961e-5;
//    static const double gKleak_d	    = 2.1987e-4;
//    static const double gleak_d	        = 10*6.2961e-5;

    /* Cortical spreading depression conditions (stimulation then vasoconstriction), also should change T_STIM_END so that the stimulus is only 1 sec long */
	static const double I_STRENGTH		= 0.006;	  	// [A] strength of current input (default 0.022)
    static const double SC_coup	        = 1;         // scaling factor for the change in SC K+ concentration based on extracellular K+ concentration (default 11.5)
    static const double Imax		    = 0.013;         // rate of the ATP pump (default 0.013*6)
    static const double gNaleak_sa	    = 9.5999e-6;       // channel conductances, change depending on Imax, see OO-NVU for other values
    static const double gKleak_sa	    = 3.4564e-5;
    static const double gleak_sa	    = 10*9.5999e-6;
    static const double gNaleak_d	    = 1.0187e-5;
    static const double gKleak_d	    = 3.4564e-5;
    static const double gleak_d	        = 10*1.0187e-5;

    static const double D_Ke      = 3.8e-9;  // [m^2/s] The diffusion rate of ECS K+
    static const double D_Nae     = 2.5e-9;  // [m^2/s] The diffusion rate of ECS Na+
	static const double D_gap	  = 3.1e-9; // [m^2/s] effective diffusion rate for K+ via gap junctions, D_Kgap = A_ef * P_K / R_k = 3.7e-9 * 5e-8 / 6e-8 = 3.1e-9 
    
    static const double wallMech	    = 1.7;          // rate of wall mechanics (default 1.7)
    static const double J_PLC 		    = 0.11; 	    // 0.11 for steady state or 0.3 for oscillations
	
    // Steady state values used for normalisation of BOLD signal and change in CBF
    //if J_PLC = 0.11
        static const double HbR_0		  = 0.6662;
        static const double CBV_0		  = 1.317;
        static const double CBF_0		  = 0.0637;
    //if J_PLC = 0.3
        //static const double HbR_0		  = 1.0753
        //static const double CBV_0		  = 0.9703
        //static const double CBF_0		  = 0.0295

/*** Model parameters ***/
        
    static const double R_decay         = 0.15;  	    // [s^-1] rate of decay of K+ in the PVS (default 0.15)

// General constants:
    static const double F             	= 96500;          // [C mol-1] Faradays constant
    static const double Farad 		  	= 96.485;         // Faradays constant in different units
    static const double R_gas         	= 8.315;          // [J mol-1K-1]
    static const double Temp          	= 300;            // [K]
    static const double gam		      	= 1970;  		  // mV/uM The change in membrane potential by a scaling factor
    static const double ph			    = 26.6995;		  // RT/Farad where Farad is in [C/mmol]
    static const double z_Na          	= 1;              // [-]
    static const double z_K           	= 1;              // [-]
    static const double z_Cl         	= -1;             // [-]
    static const double z_NBC         	= -1;             // [-]
    static const double z_Ca			= 2;			  // [-]
    static const double z_HCO3			= -1;

// NE & AC constants:
    static const double VR_sa 			= 0.465;			// [-] volume ratio between SC and AC = R_s/R_k = 2.79e-8 / 8e-8

    static const double K_Na_k        	= 10e3;           // [uM]
    static const double K_K_s         	= 1.5e3;          // [uM]

    static const double G_BK_k        	= 10.25;  		  // [uM mV^-1 s^-1]
    static const double G_K_k         	= 6907.77;
    static const double G_KCC1_k      	= 1.728;
    static const double G_NBC_k       	= 130.74;
    static const double G_Cl_k        	= 151.93;
    static const double G_NKCC1_k     	= 9.568;
    static const double G_Na_k        	= 226.94;
    static const double J_NaK_max     	= 2.3667e4;        // [uM s-1]

// Perivascular Space constants:
    static const double K_p_min 	  	= 3e3;  	// uM

// BK channel constants:
    static const double A_ef_k        	= 3.7e-9;  				// m2       Area of an endfoot of an astrocyte, equal to Area astrocyte at synaptic cleft
    static const double v_4           	= 8; 					// mV        A measure of the spread of the distribution
    static const double v_5			    = 15;					// mV
    static const double v_6			    = -55;					// mV
    static const double psi_w         	= 2.664;  				// s-1      A characteristic time
    static const double VR_pa         	= 0.001; 				// [-]       The estimated volume ratio of perivascular space to astrocyte: Model estimation
    static const double VR_ps         	= 0.001;  				// [-]       The estimated volume ratio of perivascular space to SMC: Model Estimation

// SMC constants:
    static const double F_il 		  	= 7.5e2; 		//[-] scaling factor to fit the experimental data of Filosa
    static const double z_1 		  	= 4.5e-3; 			//[-] parameter fitted on experimental data of Filosa
    static const double z_2 		  	= 112; 		//[-] parameter fitted on experimental data of Filosa
    static const double z_3 		  	= 4.2e-4; 		//[-] parameter fitted on experimental data of Filosa
    static const double z_4 		  	= 12.6; 		//[-] parameter fitted on experimental data of Filosa
    static const double z_5 		  	= -7.4e-2;  		//[-] parameter fitted on experimental data of Filosa
    static const double Fmax_i		  	= 0.23; 			// (microM/s)
    static const double Kr_i 		  	= 1;  			// (microM) Half saturation constant for agonist-dependent Ca entry
    static const double G_Ca		  	= 0.00129; 		// (microM/mV/s)
    static const double v_Ca1		  	= 100; 			// (mV)
    static const double v_Ca2		  	= -24; 			// (mV)
    static const double R_Ca		  	= 8.5; 			// (mV)
    static const double G_NaCa		  	= 0.00316; 		// (microM/mV/s)
    static const double c_NaCa		  	= 0.5; 			// (microM)
    static const double v_NaCa		  	= -30;
    static const double B_i		      	= 2.025;
    static const double cb_i		  	= 1;
    static const double C_i		      	= 55;
    static const double sc_i		  	= 2;
    static const double cc_i		  	= 0.9;
    static const double D_i		      	= 0.24;
    static const double vd_i		  	= -100;
    static const double Rd_i		  	= 250;
    static const double L_i		      	= 0.025;
    static const double F_NaK		  	= 0.0432;
    static const double G_Cl		  	= 0.00134;
    static const double v_Cl		  	= -25;
    static const double G_K		      	= 0.00446;
    static const double vK_i		  	= -94;
    static const double lam 		  	= 45;
    static const double v_Ca3		  	= -27;  			// correct
    static const double R_K		      	= 12;
    static const double k_i		      	= 0.1;

// Stretch-activated channels
    static const double G_stretch     	= 0.0061;        // uM mV-1 s-1   (stretch activated channels)
    static const double Esac          	= -18;           // mV
    static const double alpha1        	= 0.0074;
    static const double sig0          	= 500;

// EC constants:
    static const double Fmax_j		  	= 0.23; 			// [microM/s]
    static const double Kr_j		  	= 1;
    static const double B_j 		  	= 0.5;
    static const double cb_j		  	= 1;
    static const double C_j		      	= 5;
    static const double sc_j		  	= 2;
    static const double cc_j		  	= 0.9;
    static const double D_j		      	= 0.24;
    static const double L_j		      	= 0.025;
    static const double G_cat 		  	= 0.66e-3;
    static const double E_Ca		  	= 50;
    static const double m3cat		  	= -0.18;  		//-6.18 changed value!
    static const double m4cat 		  	= 0.37;
    static const double JO_j 		  	= 0.029;  		//constant Ca influx (EC)
    static const double C_m 		  	= 25.8;
    static const double G_tot		  	= 6927;
    static const double vK_j 		  	= -80;
    static const double a1_j		  	= 53.3;
    static const double a2_j		  	= 53.3;
    static const double b_j			  	= -80.8;
    static const double c_j 		  	= -0.4;  		//-6.4 changed value!
    static const double m3b		      	= 1.32e-3;
    static const double m4b		      	= 0.3;
    static const double m3s		      	= -0.28;
    static const double m4s		      	= 0.389;
    static const double G_R		      	= 955;
    static const double v_rest		  	= -31.1;
    static const double k_j		      	= 0.1;
    static const double g_hat         	= 0.5;
    static const double p_hat         	= 0.05;
    static const double p_hatIP3      	= 0.05;
    static const double C_Hillmann    	= 1;
    static const double K3_c          	= 0.4;
    static const double K4_c          	= 0.1;
    static const double K7_c          	= 0.1;
    static const double gam_cross     	= 17;
    static const double LArg_j		  	= 100;

// NO pathway
    static const double LArg          = 100;
    static const double V_spine       = 8e-8;
    static const double k_ex          = 1600;
    static const double Ca_rest       = 0.1;
    static const double lambda        = 20;
    static const double V_maxNOS      = 25e-3;
    static const double V_max_NO_n    = 4.22;
    static const double K_mO2_n 	  = 243;
    static const double K_mArg_n	  = 1.5;
    static const double K_actNOS      = 9.27e-2;
    static const double D_NO 	      = 3300;
    static const double k_O2          = 9.6e-6;
    static const double On            = 200;
    static const double v_n           = -0.04;
    static const double Ok            = 200;
    static const double G_M           = 46000;
    static const double dist_nk       = 25;
    static const double dist_ki       = 25;
    static const double dist_ij       = 3.75;
    static const double tau_nk        = 0.09469697; //pow(25,2)/(2*3300);
    static const double tau_ki        = 0.09469697; //pow(25,2)/(2*3300);
    static const double tau_ij        = 0.00213068; //pow(3.75,2)/(2*3300);
    static const double P_Ca_P_M      = 3.6;
    static const double Ca_ex         = 2e3;
    static const double M_mono        = 1.3e5;
    static const double betA          = 650;
    static const double betB          = 2800;
    static const double K_dis         = 9e-2;
    static const double K_eNOS        = 4.5e-1;
    static const double mu2           = 0.0167;
    static const double g_max         = 0.06;
    static const double alp           = 2;
    static const double W_0           = 1.4;
    static const double delt_wss      = 2.86;
    static const double k_dno         = 0.01;
    static const double k1_i          = 2e3;
    static const double k2_i          = 0.1;
    static const double k3_i          = 3;
    static const double k_1_i         = 100;
    static const double V_max_sGC     = 0.8520;   //\muM s{-1}  (for m   2)
    static const double k_pde         = 0.0195; // s{-1} (for m   2)
    static const double C_4           = 0.011;  // [s{-1} microM{-2}] (note: the changing units are correct!) (for m   2)
    static const double K_m_pde       = 2;            		// [microM]
    static const double k_mlcp_b      = 0.0086;          // [s{-1}]
    static const double k_mlcp_c      = 0.0327;           //[s{-1}]
    static const double K_m_mlcp      = 5.5;         		// [microM]
    static const double bet_i         = 0.13;  // translation factor for membrane potential dependence of KCa channel activation sigmoidal [microM2]
    static const double m_4			  = 2;
    static const double gam_eNOS      = 0.1;  // [-]
    static const double K_mO2_j       = 7.7;
    static const double V_NOj_max     = 1.22;
    static const double K_mArg_j      = 1.5;

// AC Ca2+
    static const double r_buff	        = 0.05;
    static const double G_TRPV_k		= 50;
    static const double g_TRPV_k   	    = 50 * 1e-12 / 3.7e-9;
    static const double J_max			= 2880;
    static const double K_act			= 0.17;
    static const double K_I			    = 0.03;
    static const double P_L			    = 0.0804;
    static const double k_pump			= 0.24;
    static const double V_max			= 20;
    static const double C_astr_k		= 40;
    static const double B_ex 			= 11.35;
    static const double BK_end			= 40;
    static const double K_ex			= 0.26;
    static const double delta			= 1.235e-2;
    static const double K_G			    = 8.82;
    static const double Ca_3			= 0.4;
    static const double Ca_4			= 0.35;
    static const double eet_shift		= 2;
    static const double gam_cae_k		= 200;
    static const double gam_cai_k		= 0.01;
    static const double R_0_passive_k	= 20e-6;
    static const double epshalf_k		= 0.1;
    static const double kappa_k		    = 0.1;
    static const double v1_TRPV_k		= 120;
    static const double v2_TRPV_k		= 13;
    static const double t_TRPV_k		= 0.9;
    static const double VR_ER_cyt		= 0.185;
    static const double K_inh			= 0.1;
    static const double k_on			= 2;
    static const double k_deg			= 1.25;
    static const double r_h			    = 4.8;
    static const double Ca_k_min		= 0.1;
    static const double k_eet			= 7.2;
    static const double V_eet			= 72;
    static const double Ca_decay_k		= 0.5;
    static const double Capmin_k		= 2000;
    static const double m_c			    = 4;

// Glutamate constants
    static const double Glu_max		    = 1846;
    static const double Glu_slope 		= 0.1;
    static const double Ke_switch		= 5.5;
    static const double rho_min		    = 0.1;
    static const double rho_max		    = 0.7;

// Neuron constants
    static const double E_Cl_sa	        = -70;
    static const double E_Cl_d		    = -70;
    static const double Ra			    = 1.83e5;
    static const double dhod 		    = 4.5e-2;
    static const double As 		        = 1.586e-5;
    static const double Ad			    = 2.6732e-4;
    static const double Vs			    = 2.16e-9;
    static const double Vd			    = 5.614e-9;
    static const double fe			    = 0.15;
    static const double Cm			    = 7.5e-7;
    static const double Mu			    = 8e-4;
    static const double Buff0	        = 500;
    static const double gNaP_GHk	    = 2e-6;
    static const double gKDR_GHk	    = 10e-5;
    static const double gKA_GHk	        = 1e-5;
    static const double gNMDA_GHk	    = 1e-5;
    static const double gNaT_GHk        = 10e-5;
    static const double O2_0		    = 0.02;
    static const double alpha_O2        = 0.05;
    static const double D_Na 		    = 1.33e-5;
    static const double D_K 		    = 1.96e-5;
    static const double K_init_e 	    = 2.9;
    static const double Na_init_sa      = 10;
    static const double Na_init_d 	    = 10;
    static const double R_init	        = 1.9341e-5;
    static const double CBF_init 	    = 0.032;
    static const double O2_b 		    = 0.04;
    static const double gamma_O2        = 0.1;
    static const double Mg			    = 1.2;

// BOLD constants
    static const double tau_MTT	        = 3;
    static const double tau_TAT	        = 20;
    static const double d_BOLD		    = 0.4;
    static const double E_0		        = 0.4;
    static const double a_1		        = 3.4;
    static const double a_2		        = 1;
    static const double V_0		        = 0.03;

// Pressure constants
    static const double HRR         = 0.1;             // Nondimensional (thickness to radius ratio)
    static const double RSCALE      = 0.6;            // Dimensionless
    static const double E0          = 66e3;            // Pa
    static const double EPASSIVE    = 66e3;           // Pa
    static const double EACTIVE     = 233e3;           // Pa
    static const double ETA         = 2.8e2;          // Pa s
    static const double T0          = 1;               // s
    static const double PA2MMHG     = 0.00750061683;   // convert from Pa to mmHg

/*** Diffusion constants ***/
// Number of variables stored in diffusion structs.
    static const int NUM_DIFF_VARS = 8; //2;

// Rates of diffusion (characteristic times) for diffusion between tissue blocks
//    static const double tau_Ke      = 4.3;  // (sec) The diffusion rate - characteristic time scale for K+ to travel one NVU block
//    static const double tau_Nae     = 6.4;  // (sec) The diffusion rate - characteristic time scale for Na+ to travel one NVU block
    

// Gap junction constants
    static const double delta_x 	= 1.24e-4; // [m] length/width of one NVU block

/*** H-tree constants ***/
// Constants for the H-tree, don't change
    static const double RMIN                = 10e-6;   					    	// m, radius of smallest vessel
    static const double BIFURCATION_SCALE   = 1.4142135623730951; 	            //   sqrt(2), amount the radius decreases by when going down a level
    static const double L0                  = 200e-6;   						// m (for nondimensionalising), length characteristic value
    static const double LRR                 = 20;   							// Nondimensional, length to radius ratio
    static const double MU                  = 3.5e-3;   						// Pa s, blood viscosity
    static const double R0                  = 10e-6;                           	// m (for nondimensionalising)
    static const double P0                  = 8000;                            	// Pa (scaling factor for nondim)
    static const double PCAP                = 4000;                            	// Pa (capillary bed pressure)

/*** State variable indexing ***/
//only needs to be changed when more variables are added
	static const int i_radius  = 0; // radius has to be 0, this is assumed elsewhere

	// AC
	static const int i_v_k       = 1;
	static const int i_Na_k    = 2;
	static const int i_K_k     = 3;
	static const int i_HCO3_k  = 4;
	static const int i_Cl_k    = 5;
	static const int i_w_k       = 10;

	// SC
	static const int i_Na_s    = 6;
	static const int i_K_s     = 7;
	static const int i_HCO3_s  = 8;

	// PVS
	static const int i_K_p       = 9;

	// SMC
	static const int i_ca_i      = 11;
	static const int i_ca_sr_i   = 12;
	static const int i_v_i       = 13;
	static const int i_w_i       = 14;
	static const int i_ip3_i     = 15;
	static const int i_K_i       = 16;

	// EC
	static const int i_ca_j      = 17;
	static const int i_ca_er_j   = 18;
	static const int i_v_j       = 19;
	static const int i_ip3_j     = 20;

	// Mech
	static const int i_Mp        = 21;
	static const int i_AMp       = 22;
	static const int i_AM        = 23;

	// NO pathway
	static const int i_NOn        = 24;
	static const int i_NOk        = 25;
	static const int i_NOi        = 26;
	static const int i_NOj        = 27;
	static const int i_cGMP       = 28;
	static const int i_eNOS       = 29;
	static const int i_nNOS       = 30;
	static const int i_ca_n       = 31;
	static const int i_E_b        = 32;
	static const int i_E_6c       = 33;

	// AC Ca2+
	static const int i_ca_k       = 34;
	static const int i_s_k        = 35;
	static const int i_h_k        = 36;
	static const int i_ip3_k      = 37;
	static const int i_eet_k      = 38;
	static const int i_m_k        = 39;
	static const int i_ca_p       = 40;

	// Neuron - ions
	static const int i_v_sa	   = 41;
	static const int i_v_d	   = 42;
	static const int i_K_sa	   = 43;
	static const int i_Na_sa   = 44;
	static const int i_K_d	   = 45;
	static const int i_Na_d	   = 46;
	static const int i_K_e	   = 47;
	static const int i_Na_e	   = 48;

	// Neuron - other
	static const int i_Buff_e	   = 49;
	static const int i_O2		   = 50;
	static const int i_CBV	   = 51;
	static const int i_HbR	   = 52;

	// Neuron Gating Variables
	static const int i_m1	   	   = 53;
	static const int i_m2	   	   = 54;
	static const int i_m3	   	   = 55;
	static const int i_m4	   	   = 56;
	static const int i_m5	   	   = 57;
	static const int i_m6	   	   = 58;
	static const int i_m7	   	   = 59;
	static const int i_m8	   	   = 60;
	static const int i_h1	   	   = 61;
	static const int i_h2	   	   = 62;
	static const int i_h3	   	   = 63;
	static const int i_h4	   	   = 64;
	static const int i_h5	   	   = 65;
	static const int i_h6	   	   = 66;
	
	// Curvature variables (not time dependent)
	static const int i_curvature  = 67;
	static const int i_coup = 68;

// Number of ODEs
    static const int NEQ       = 69;



#endif /* SRC_CONSTANTS_H_ */
