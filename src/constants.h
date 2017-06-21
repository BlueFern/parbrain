#ifndef SRC_CONSTANTS_H_
#define SRC_CONSTANTS_H_

// Optional command line arguments for parBrainSim: N, NSUB, T_FINAL (in that order). If none specified then the following are used.

// TODO: switches (TRPV, NO, K, Ca etc), ECS input size, type of spatial input (centre vs corner etc)


/*** Run parameters ***/
    #define T_FINAL         100       // Final run time
	#define T_STIM_0        20        // Start time for stimulation
	#define T_STIM_END      36        // End time for stimulation
    #define DT_WRITE        0.1       // Time step for writing to file (and screen)
    #define NTREE           3         // Number of levels in the H-tree (where the tissue slice has 2^(N-1) tissue blocks)
    #define NSUB            1         // Subtree size (usually 3 but must be smaller when N is small - adjust as needed)
	#define P_TOP			4100	  // Pressure at the top of the tree, chosen so that the drop over the terminating arterioles is around 18.2 Pa to match with the single NVU model.
									  // For NTREE=3, P_TOP=4100 Pa. For NTREE=7, P_TOP=4160 Pa. For NTREE=13, P_TOP=?

/*** Commonly changed model parameters ***/

	#define I_STRENGTH		0.025	  	// [A] strength of current input
    #define wallMech	    1.1         // rate of wall mechanics, 1 for normal (default 1.1)
    #define SC_coup	        11.5        // scaling factor for the change in SC K+ concentration based on extracellular K+ concentration (default 11.5)
    #define J_PLC 		    0.11 	    // 0.11 for steady state or 0.3 for oscillations
    #define R_decay         0.15  	    // [s^-1] rate of decay of K+ in the PVS (default 0.15)

    #define Imax		    0.013*6         // rate of the ATP pump (default 0.013*6)
    #define gNaleak_sa	    6.2378e-5       // channel conductances, change depending on Imax, see OO-NVU for other values
    #define gKleak_sa	    2.1989e-4
    #define gleak_sa	    10*6.2378e-5
    #define gNaleak_d	    6.2961e-5
    #define gKleak_d	    2.1987e-4
    #define gleak_d	        10*6.2961e-5

    // Steady state values used for normalisation of BOLD signal and change in CBF
    //if J_PLC = 0.11
        #define DHG_0		  0.6662
        #define CBV_0		  1.317
        #define CBF_0		  0.0637
    //if J_PLC = 0.3
        //#define DHG_0		  1.0753
        //#define CBV_0		  0.9703
        //#define CBF_0		  0.0295

/*** Model parameters ***/

// general constants:
    #define F             96500          // [C mol-1] Faradays constant
    #define Farad 		  96.485         // Faradays constant in different unit
    #define R_gas         8.315          // [J mol-1K-1]
    #define Temp          300            // [K]
    #define unitcon       1e3            // [-] Factor to convert equations to another unit.

// NE & AC constants:
    #define L_p           2.1e-9         // [m uM-1s-1]
    #define R_tot         8.79e-8        // [m]   total volume surface area ratio AC+SC  **see nvu.h
    #define X_k           12.41e-3       // [uMm]
    #define z_Na          1              // [-]
    #define z_K           1              // [-]
    #define z_Cl          -1             // [-]
    #define z_NBC         -1             // [-]
    #define g_K_k         40             // [ohm-1m-2]
    #define g_KCC1_k      1e-2           // [ohm-1m-2]
    #define g_NBC_k       7.57e-1        // [ohm-1m-2]
    #define g_Cl_k        8.797e-1       // [ohm-1m-2]
    #define g_NKCC1_k     5.54e-2        // [ohm-1m-2]
    #define g_Na_k        1.314          // [ohm-1m-2]
    #define J_NaK_max     1.42e-3        // [uMm s-1]
    #define K_Na_k        10e3           // [uM]
    #define K_K_s         1.5e3          // [uM]

// Perivascular Space constants:
    #define K_p_min 	  3e3  	// uM

// BK channel constants:
    #define A_ef_k        3.7e-9  				// m2       Area of an endfoot of an astrocyte, equal to Area astrocyte at synaptic cleft
    #define v_4           8e-3 				// V        A measure of the spread of the distribution
    #define psi_w         2.664  				// s-1      A characteristic time
    #define G_BK_k        225  				// !!!
    #define g_BK_k        G_BK_k * 1e-12 / A_ef_k  	// ohm-1m-2  Specific capacitance of the BK-Channel in units of Ostby
    #define VR_pa         0.001 				// [-]       The estimated volume ratio of perivascular space to astrocyte: Model estimation
    #define VR_ps         0.001  				// [-]       The estimated volume ratio of perivascular space to SMC: Model Estimation

// SMC constants:
    #define F_il 		  7.5e2 		//[-] scaling factor to fit the experimental data of Filosa
    #define z_1 		  4.5 			//[-] parameter fitted on experimental data of Filosa
    #define z_2 		  -1.12e2 		//[-] parameter fitted on experimental data of Filosa
    #define z_3 		  4.2e-1 		//[-] parameter fitted on experimental data of Filosa
    #define z_4 		  -1.26e1 		//[-] parameter fitted on experimental data of Filosa
    #define z_5 		  -7.4e-2  		//[-] parameter fitted on experimental data of Filosa
    #define Fmax_i		  0.23 			// (microM/s)
    #define Kr_i 		  1  			// (microM) Half saturation constant for agonist-dependent Ca entry
    #define G_Ca		  0.00129 		// (microM/mV/s)
    #define v_Ca1		  100 			// (mV)
    #define v_Ca2		  -24 			// (mV)
    #define R_Ca		  8.5 			// (mV)
    #define G_NaCa		  0.00316 		// (microM/mV/s)
    #define c_NaCa		  0.5 			// (microM)
    #define v_NaCa		  -30
    #define B_i		      2.025
    #define cb_i		  1
    #define C_i		      55
    #define sc_i		  2
    #define cc_i		  0.9
    #define D_i		      0.24
    #define vd_i		  -100
    #define Rd_i		  250
    #define L_i		      0.025
    #define gam		      1970  		// mVmicroM-1 The change in membrane potential by a scaling factor
    #define F_NaK		  0.0432
    #define G_Cl		  0.00134
    #define v_Cl		  -25
    #define G_K		      0.00446
    #define vK_i		  -94
    #define lam 		  45
    #define v_Ca3		  -27  			// correct
    #define R_K		      12
    #define k_i		      0.1

// Stretch-activated channels
    #define G_stretch     0.0061        // uM mV-1 s-1   (stretch activated channels)
    #define Esac          -18           // mV
    #define alpha1        0.0074
    #define sig0          500

// EC constants:
    #define Fmax_j		  0.23 			// [microM/s]
    #define Kr_j		  1
    #define B_j 		  0.5
    #define cb_j		  1
    #define C_j		      5
    #define sc_j		  2
    #define cc_j		  0.9
    #define D_j		      0.24
    #define L_j		      0.025
    #define G_cat 		  0.66e-3
    #define E_Ca		  50
    #define m3cat		  -0.18  		//-6.18 changed value!
    #define m4cat 		  0.37
    #define JO_j 		  0.029  		//constant Ca influx (EC)
    #define C_m 		  25.8
    #define G_tot		  6927
    #define vK_j 		  -80
    #define a1_j		  53.3
    #define a2_j		  53.3
    #define b_j			  -80.8
    #define c_j 		  -0.4  		//-6.4 changed value!
    #define m3b		      1.32e-3
    #define m4b		      0.3
    #define m3s		      -0.28
    #define m4s		      0.389
    #define G_R		      955
    #define v_rest		  -31.1
    #define k_j		      0.1
    #define g_hat         0.5
    #define p_hat         0.05
    #define p_hatIP3      0.05
    #define C_Hillmann    1
    #define K3_c          0.4 * C_Hillmann
    #define K4_c          0.1 * C_Hillmann
    #define K7_c          0.1 * C_Hillmann
    #define gam_cross     17 * C_Hillmann
    #define LArg_j		  100

// NO pathway
    #define LArg          100
    #define V_spine       8e-8
    #define k_ex          1600
    #define Ca_rest       0.1
    #define lambda        20
    #define V_maxNOS      25e-3
    #define V_max_NO_n    4.22
    #define K_mO2_n 	  243
    #define K_mArg_n	  1.5
    #define K_actNOS      9.27e-2
    #define D_NO 	      3300
    #define k_O2          9.6e-6
    #define On            200
    #define v_n           -0.04
    #define Ok            200
    #define G_M           46000
    #define dist_nk       25
    #define dist_ki       25
    #define dist_ij       3.75
    #define tau_nk        pow(dist_nk,2)/(2*D_NO)
    #define tau_ki        pow(dist_ki,2)/(2*D_NO)
    #define tau_ij        pow(dist_ij,2)/(2*D_NO)
    #define P_Ca_P_M      3.6
    #define Ca_ex         2e3
    #define M_mono        1.3e5
    #define betA          650
    #define betB          2800
    #define Oj            200
    #define K_dis         9e-2
    #define K_eNOS        4.5e-1
    #define mu2           0.0167
    #define g_max         0.06
    #define alp           2
    #define W_0           1.4
    #define delt_wss      2.86
    #define k_dno         0.01
    #define k1_i          2e3
    #define k2_i          0.1
    #define k3_i          3
    #define k_1_i         100
    #define V_max_sGC     0.8520   //\muM s{-1}  (for m   2)
    #define k_pde         0.0195 // s{-1} (for m   2)
    #define C_4           0.011  // [s{-1} microM{-2}] (note: the changing units are correct!) (for m   2)
    #define K_m_pde       2            		// [microM]
    #define k_mlcp_b      0.0086          // [s{-1}]
    #define k_mlcp_c      0.0327           //[s{-1}]
    #define K_m_mlcp      5.5         		// [microM]
    #define bet_i         0.13  // translation factor for membrane potential dependence of KCa channel activation sigmoidal [microM2]
    #define m_4			  2
    #define gam_eNOS      0.1  // [-]
    #define K_mO2_j       7.7
    #define V_NOj_max     1.22
    #define K_mArg_j      1.5

// AC Ca2+
    #define r_buff	        0.05
    #define G_TRPV_k		50
    #define g_TRPV_k   	    G_TRPV_k * 1e-12 / A_ef_k
    #define J_max			2880
    #define K_act			0.17
    #define K_I			    0.03
    #define P_L			    0.0804
    #define k_pump			0.24
    #define V_max			20
    #define C_astr_k		40
    #define gamma_k		    834.3
    #define B_ex 			11.35
    #define BK_end			40
    #define K_ex			0.26
    #define delta			1.235e-2
    #define K_G			    8.82
    #define Ca_3			0.4
    #define Ca_4			0.35
    #define v_5			    15e-3
    #define v_7			    -55e-3
    #define eet_shift		2e-3
    #define gam_cae_k		200
    #define gam_cai_k		0.01
    #define R_0_passive_k	20e-6
    #define epshalf_k		0.1
    #define kappa_k		    0.1
    #define v1_TRPV_k		0.12
    #define v2_TRPV_k		0.013
    #define t_TRPV_k		0.9
    #define VR_ER_cyt		0.185
    #define K_inh			0.1
    #define k_on			2
    #define k_deg			1.25
    #define r_h			    4.8
    #define Ca_k_min		0.1
    #define k_eet			7.2
    #define V_eet			72
    #define Ca_decay_k		0.5
    #define Capmin_k		2000
    #define reverseBK		0
    #define switchBK		1
    #define trpv_switch	    1
    #define z_Ca			2
    #define m_c			    4

// Glutamate constants
    #define Glu_max		    1846
    #define Glu_slope 		0.1
    #define Ke_switch		5.5
    #define rho_min		    0.1
    #define rho_max		    0.7

// Neuron constants
    #define E_Cl_sa	        -70
    #define E_Cl_d		    -70
    #define Ra			    1.83e5
    #define dhod 		    4.5e-2
    #define As 		        1.586e-5
    #define Ad			    2.6732e-4
    #define Vs			    2.16e-9
    #define Vd			    5.614e-9
    #define fe			    0.15
    #define Cm			    7.5e-7
    #define ph			    26.6995
    #define Mu			    8e-4
    #define Buff0	        500
    #define gNaP_GHk	    2e-6
    #define gKDR_GHk	    10e-5
    #define gKA_GHk	        1e-5
    #define gNMDA_GHk	    1e-5
    #define gNaT_GHk        10e-5
    #define O2_0		    0.02
    #define alpha_O2        0.05
    #define D_Na 		    1.33e-5
    #define D_K 		    1.96e-5
    #define K_init_e 	    2.9
    #define Na_init_sa      10
    #define Na_init_d 	    10
    #define R_init	        1.9341e-5
    #define CBF_init 	    0.032
    #define O2_b 		    0.04
    #define gamma_O2        0.1
    #define Mg			    1.2

// BOLD constants
    #define tau_MTT	        3
    #define tau_TAT	        20
    #define d_BOLD		    0.4
    #define E_0		        0.4
    #define a_1		        3.4
    #define a_2		        1
    #define V_0		        0.03

// Pressure constants
    #define HRR         0.1             // Nondimensional (thickness to radius ratio)
    #define RSCALE      0.6             // Dimensionless
    #define E0          66e3            // Pa
    #define EPASSIVE    66e3            // Pa
    #define EACTIVE     233e3           // Pa
    #define ETA         2.8e2           // Pa s
    #define T0          1               // s
    #define PA2MMHG     0.00750061683   // convert from Pa to mmHg

// Rates of diffusion (characteristic times) for diffusion between tissue blocks
    #define tau_Ke      4.3  // (sec) The diffusion rate - characteristic time scale for K+ to travel one NVU block
    #define tau_Nae     6.4  // (sec) The diffusion rate - characteristic time scale for Na+ to travel one NVU block

/*** H-tree constants ***/
// Constants for the H-tree, don't change
    #define RMIN                10e-6   					    // m, radius of smallest vessel
    #define BIFURCATION_SCALE   1.4142135623730951 	            //   sqrt(2), amount the radius decreases by when going down a level
    #define L0                  200e-6   						// m (for nondimensionalising), length characteristic value
    #define LRR                 20   							// Nondimensional, length to radius ratio
    #define MU                  3.5e-3   						// Pa s, blood viscosity
    #define R0                  10e-6                           // m (for nondimensionalising)
    #define P0                  8000                            // Pa (scaling factor for nondim)
    #define PCAP                4000                            // Pa (capillary bed pressure)

/*** State variable indexing ***/
//only needs to be changed when more variables are added
    #define i_radius    0  // radius has to be 0, this is assumed elsewhere

// AC
    #define R_k         1
    #define N_Na_k      2
    #define N_K_k       3
    #define N_HCO3_k    4
    #define N_Cl_k      5
    #define w_k         10

// SC
    #define N_Na_s      6
    #define N_K_s       7
    #define N_HCO3_s    8

// PVS
    #define K_p         9

// SMC
    #define ca_i        11
    #define ca_sr_i     12
    #define v_i         13
    #define w_i         14
    #define ip3_i       15
    #define K_i         16

// EC
    #define ca_j        17
    #define ca_er_j     18
    #define v_j         19
    #define ip3_j       20

// Mech
    #define Mp          21
    #define AMp         22
    #define AM          23

// NO pathway
    #define NOn         24
    #define NOk         25
    #define NOi         26
    #define NOj         27
    #define cGMP        28
    #define eNOS        29
    #define nNOS        30
    #define ca_n        31
    #define E_b         32
    #define E_6c        33

// AC Ca2+
    #define ca_k        34
    #define s_k         35
    #define h_k         36
    #define ip3_k       37
    #define eet_k       38
    #define m_k         39
    #define ca_p        40

// Neuron - ions
    #define v_sa	    41
    #define v_d	        42
    #define K_sa	    43
    #define Na_sa	    44
    #define K_d	        45
    #define Na_d	    46
    #define K_e	        47
    #define Na_e	    48

// Neuron - other
    #define Buff_e	    49
    #define O2		    50
    #define CBV	        51
    #define DHG	        52

// Neuron Gating Variables
    #define m1	   	    53
    #define m2	   	    54
    #define m3	   	    55
    #define m4	   	    56
    #define m5	   	    57
    #define m6	   	    58
    #define m7	   	    59
    #define m8	   	    60
    #define h1	   	    61
    #define h2	   	    62
    #define h3	   	    63
    #define h4	   	    64
    #define h5	   	    65
    #define h6	   	    66

// Number of ODEs
    #define NEQ         67



#endif /* SRC_CONSTANTS_H_ */
