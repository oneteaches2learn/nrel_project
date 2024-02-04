#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include "mechanism.H"
#include <GPU_misc.H>
#include <PelePhysics.H>

#include <CoolProp.h>

int
main(int argc, char* argv[])
{
	amrex::Initialize(argc, argv);

	// STEP 1: SET UP GRID
	// declare variables
    int i,j, T_steps, P_steps;
    double T_START, T_END, P_START, P_END, T_delta, P_delta;
    std::string substance, discard_param_name;
	std::string P_unit, R_unit, mu_unit, lam_unit;

	// Read Inputs from Config File
    std::ifstream inputs("parameters.txt");
    inputs >> discard_param_name >> substance
		   >> discard_param_name >> T_delta
           >> discard_param_name >> T_steps
           >> discard_param_name >> T_START
		   >> discard_param_name >> T_END
           >> discard_param_name >> P_delta
           >> discard_param_name >> P_steps
           >> discard_param_name >> P_START
		   >> discard_param_name >> P_END
		   >> discard_param_name >> P_unit
		   >> discard_param_name >> R_unit
		   >> discard_param_name >> mu_unit
		   >> discard_param_name >> lam_unit;
 
    // Complete missing data
    if (T_steps == 0) {T_steps = (T_END - T_START) / T_delta + 1;}
    if (T_END == 0)   {T_END   = T_START + T_steps * T_delta;}
    if (T_delta == 0) {T_delta = (T_END - T_START) / T_steps;}
    if (P_steps == 0) {P_steps = (P_END - P_START) / P_delta + 1;}
    if (P_END == 0)   {P_END   = P_START + P_steps * P_delta;} 
    if (P_delta == 0) {P_delta = (P_END - P_START) / P_steps;}

	// Print state to command line
    std::cout << "\n Running PelePhysics with the following inputs..."
			  << "\n   substance: \t" << substance
              << "\n   T_delta:   \t" << T_delta
              << "\n   T_steps:   \t" << T_steps
              << "\n   T_start:   \t" << T_START
			  << "\n   T_end:     \t" << T_END
			  << "\n   T_unit:    \t" << "[K]"
              << "\n   P_delta:   \t" << P_delta
              << "\n   P_steps:   \t" << P_steps
              << "\n   P_start:   \t" << P_START
			  << "\n   P_end:     \t" << P_END
			  << "\n   P_unit:    \t" << "[" << P_unit   << "]"
			  << "\n   R_unit:    \t" << "[" << R_unit   << "]"
			  << "\n   mu_unit:   \t" << "[" << mu_unit  << "]"
			  << "\n   lam_unit:  \t" << "[" << lam_unit << "]" << std::endl
              << std::endl << std::endl;
 

	// STEP 2: SET UP PROGRAM
    // Declare Variables
	amrex::Real P_unit2cgs, R_unit2cgs, mu_unit2cgs, lam_unit2cgs;
	amrex::Real P_unit2SI, R_unit2SI, mu_unit2SI, lam_unit2SI;

	amrex::Real P, T, am, bm, Wbar, Cp, xi;
	amrex::Real R_SRK[P_steps][T_steps]; 
	amrex::Real R_SW[P_steps][T_steps]; 
	amrex::Real Z[P_steps][T_steps];
	amrex::Real lam_SRK_Chung[P_steps][T_steps]; 
	amrex::Real	lam_SRK_Huber[P_steps][T_steps]; 
	amrex::Real	lam_SW_Chung[P_steps][T_steps]; 
	amrex::Real	lam_SW_Huber[P_steps][T_steps];
	amrex::Real mu_SRK_Chung[P_steps][T_steps]; 
	amrex::Real	mu_SRK_LM[P_steps][T_steps]; 
	amrex::Real	mu_SW_Chung[P_steps][T_steps];
	amrex::Real	mu_SW_LM[P_steps][T_steps];
    amrex::Real dmudrho_SRK[P_steps][T_steps]; 
	amrex::Real	dlamdrho_SRK[P_steps][T_steps]; 
	amrex::Real	dmudrho_SW[P_steps][T_steps]; 
	amrex::Real	dlamdrho_SW[P_steps][T_steps];
	amrex::Real	R_Delta[P_steps][T_steps];
	amrex::Real	mu_Delta[P_steps][T_steps];
	amrex::Real	mu_Delta_EoS[P_steps][T_steps];
	amrex::Real	mu_Delta_cor[P_steps][T_steps];
	amrex::Real	lam_Delta[P_steps][T_steps];
	amrex::Real	lam_Delta_EoS[P_steps][T_steps];
	amrex::Real	lam_Delta_cor[P_steps][T_steps];
	amrex::Real	mu_Delta_approx1[P_steps][T_steps];
	amrex::Real	mu_product[P_steps][T_steps];
	amrex::Real	mu_Delta_approx2[P_steps][T_steps];
	amrex::Real	lam_Delta_approx1[P_steps][T_steps];
	amrex::Real	lam_product[P_steps][T_steps];
	amrex::Real	lam_Delta_approx2[P_steps][T_steps];
	
	amrex::Real H2Ofrac;
    amrex::Real Y[NUM_SPECIES];
    amrex::Real Ddiag[NUM_SPECIES]   = {0.0};
    amrex::Real chi_mix[NUM_SPECIES] = {0.0};
    const bool wtr_get_xi    = true;
    const bool wtr_get_mu    = true;
    const bool wtr_get_lam   = true;
    const bool wtr_get_Ddiag = false;
    const bool wtr_get_chi   = false;
            
    // Instantiate eos and transport objects
    pele::physics::eos::SRK srk;
    pele::physics::transport::SimpleTransport trans;

    // Get transport parameters
    amrex::ParmParse pp;
    pele::physics::transport::TransportParams<
        pele::physics::PhysicsType::transport_type>
        trans_parms;
    trans_parms.allocate();
    auto const* ltransparm = trans_parms.device_trans_parm();

    // Set mass fraction
    for (int i = 0; i < NUM_SPECIES; i++) {
        Y[i] = 0;
    }
	if (substance == "CO2") {Y[CO2_ID] = 1;}
	if (substance == "H2O") {Y[H2O_ID] = 1;}
	if (substance == "O2")  {Y[O2_ID]  = 1;}
	if (substance == "CH4") {Y[CH4_ID] = 1;}
//	if (substance == "H2")  {Y[H2_ID]  = 1;}
//	if (substance == "N2")  {Y[N2_ID]  = 1;}
//	if (substance == "HF")  {Y[HF_ID]  = 1;}

    std::cout << "\n   Y[CO2]:   \t" << Y[CO2_ID]
              << "\n   Y[H2O]:   \t" << Y[H2O_ID]
              << "\n   Y[O2]:    \t" << Y[O2_ID]
              << "\n   Y[CH4]:   \t" << Y[CH4_ID]
              << "\n   H2O_ID:   \t" << H2O_ID
//            << "\n   Y[H2]:    \t" << Y[H2_ID]
//            << "\n   Y[N2]:    \t" << Y[N2_ID]
//			  << "\n   Y[HF]:    \t" << Y[HF_ID]
			  << std::endl << std::endl;

    // Set Wbar
    srk.Y2WBAR(Y,Wbar);
	
	// SETUP UNIT CONVERSIONS
	// pressure
	if (P_unit == "barye")    {P_unit2cgs = 1; P_unit2SI = 1e-1;}
	else if (P_unit == "Pa")  {P_unit2cgs = 10; P_unit2SI = 1;}
	else if (P_unit == "MPa") {P_unit2cgs = 1e7; P_unit2SI = 1e6;}
	else if (P_unit == "atm") {P_unit2cgs = 1013250; P_unit2SI = 101325;}
	else {P_unit2cgs = 1013250; P_unit2SI = 101325;}

	// density
    if (R_unit == "g/cm3")      {R_unit2cgs = 1; R_unit2SI = 1e3;}
    else if (R_unit == "kg/m3") {R_unit2cgs = 1e-3; R_unit2SI = 1;}
    else if (R_unit == "g/mL")  {R_unit2cgs = 1; R_unit2SI = 1e3;}
    else if (R_unit == "g/L")   {R_unit2cgs = 1e-3; R_unit2SI = 1;}
    else if (R_unit == "kg/L")  {R_unit2cgs = 1; R_unit2SI = 1e3;}
    else {R_unit2cgs = 1e-3; R_unit2SI = 1;} // presumed SI units if R_unit not specified

	// viscosity
	if (mu_unit == "poise")      {mu_unit2cgs = 1; mu_unit2SI = 1e-1;}
	else if (mu_unit == "Pa_s")  {mu_unit2cgs = 10; mu_unit2SI = 1;}
	else if (mu_unit == "mPa_s") {mu_unit2cgs = 0.01; mu_unit2SI = 1e-3;}
	else if (mu_unit == "uPa_s") {mu_unit2cgs = 1e-5; mu_unit2SI = 1e-6;}
	else {mu_unit2cgs = 10; mu_unit2SI = 1;}

	// thermal conductivity
	if (lam_unit == "W/m/K")       {lam_unit2cgs = 1e6; lam_unit2SI = 1;}
	else if (lam_unit == "W/cm/K") {lam_unit2cgs = 1e5; lam_unit2SI = 1e-1;}
	else {lam_unit2cgs = 1e6; lam_unit2SI = 1;}


    // STEP 3: WRITE INDEPENDENT VARIABLES
    // write temperature
    std::ofstream myfile;
    myfile.open("IV_temperature.txt");
    for (i = 0; i < T_steps; i++){
        T = T_START + i*T_delta;
        myfile << T; 
        if (i != T_steps-1) {myfile << ",";}
    }
    myfile.close();

    // write pressure
    myfile.open ("IV_pressure.txt");
    for (i = 0; i < P_steps; i++){
        P = P_START + i*P_delta; 
        myfile << P;
        if (i != P_steps-1) {myfile << ",";}
    }
    myfile.close();


	// STEP 4: STORE DEPENDENT VARIABLES
    // Loop over pressure and temperature. write results.
    for (int j = 0; j < P_steps; j++) {
        P = P_START + j*P_delta;
        for (int i = 0; i < T_steps; i++) {
            T = T_START + i*T_delta;

			// step 1: compute pelephysics values from SRK densities
            srk.MixingRuleAmBm(T,Y,am,bm);
            srk.Calc_CompressFactor_Z(Z[j][i],am,bm,P*P_unit2cgs,T,Wbar);
            srk.PYT2R(P*P_unit2cgs,Y,T,R_SRK[j][i]); // R in cgs, stays in cgs until "unit conversions" below          
            srk.RTY2Cp(R_SRK[j][i],T,Y,Cp);

            trans.transport(
                wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, wtr_get_chi, T, R_SRK[j][i],
                Y, Ddiag, chi_mix, mu_SRK_Chung[j][i], xi, lam_SRK_Chung[j][i], 
				dmudrho_SRK[j][i], dlamdrho_SRK[j][i], ltransparm);

			// unit conversion from CGS to desired unit
			R_SRK[j][i] 		= R_SRK[j][i] / R_unit2cgs;
			mu_SRK_Chung[j][i]  = mu_SRK_Chung[j][i]   / mu_unit2cgs;
			lam_SRK_Chung[j][i] = lam_SRK_Chung[j][i]  / lam_unit2cgs;
			dmudrho_SRK[j][i]   = dmudrho_SRK[j][i]    * R_unit2cgs / mu_unit2cgs;  
			dlamdrho_SRK[j][i]  =	dlamdrho_SRK[j][i] * R_unit2cgs / lam_unit2cgs;	
			
			// step 2: compute coolprop values 
	    	R_SW[j][i]          = CoolProp::PropsSI("Dmass","T",T,"P",P * P_unit2SI,substance);
	    	lam_SW_Huber[j][i]  = CoolProp::PropsSI("L","T",T,"P",P * P_unit2SI,substance);
	    	lam_SRK_Huber[j][i] = CoolProp::PropsSI("L","T",T,"D",R_SRK[j][i] * R_unit2SI,substance);
	    	mu_SW_LM[j][i]      = CoolProp::PropsSI("V","T",T,"P",P * P_unit2SI,substance);
	    	mu_SRK_LM[j][i]     = CoolProp::PropsSI("V","T",T,"D",R_SRK[j][i] * R_unit2SI,substance);
           
			// unit conversion from CoolProp units to desired unit
			R_SW[j][i] 			= R_SW[j][i] 		  / R_unit2SI;
			lam_SW_Huber[j][i] 	= lam_SW_Huber[j][i]  / lam_unit2SI;
			lam_SRK_Huber[j][i] = lam_SRK_Huber[j][i] / lam_unit2SI;
			mu_SW_LM[j][i] 		= mu_SW_LM[j][i] 	  / mu_unit2SI;
			mu_SRK_LM[j][i] 	= mu_SRK_LM[j][i] 	  / mu_unit2SI;
 
			// step 3: compute pelephysics values from SW densities
            trans.transport(
                wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, wtr_get_chi, T, 
				R_SW[j][i] * R_unit2cgs, Y, Ddiag, chi_mix, mu_SW_Chung[j][i], xi, 
				lam_SW_Chung[j][i], dmudrho_SW[j][i], dlamdrho_SW[j][i], ltransparm);
			
			// unit conversion from CGS to desired unit
			mu_SW_Chung[j][i]  = mu_SW_Chung[j][i]  / mu_unit2cgs;
			lam_SW_Chung[j][i] = lam_SW_Chung[j][i] / lam_unit2cgs;
			dmudrho_SW[j][i]   = dmudrho_SW[j][i]   * R_unit2cgs / mu_unit2cgs;  
			dlamdrho_SW[j][i]  = dlamdrho_SW[j][i]  * R_unit2cgs / lam_unit2cgs;	

			// step 4: compute errors
			R_Delta[j][i]       = R_SRK[j][i] 		- R_SW[j][i];
			mu_Delta[j][i]      = mu_SRK_Chung[j][i]  - mu_SW_LM[j][i];
			mu_Delta_EoS[j][i]  = mu_SRK_LM[j][i]     - mu_SW_LM[j][i];
			mu_Delta_cor[j][i]  = mu_SW_Chung[j][i]   - mu_SW_LM[j][i];
			lam_Delta[j][i]     = lam_SRK_Chung[j][i] - lam_SW_Huber[j][i];
			lam_Delta_EoS[j][i] = lam_SRK_Huber[j][i] - lam_SW_Huber[j][i];
			lam_Delta_cor[j][i] = lam_SW_Chung[j][i]  - lam_SW_Huber[j][i];
			
			// step 5: compute additional quatities of interest
			mu_Delta_approx1[j][i] = mu_Delta_EoS[j][i] + mu_Delta_cor[j][i];
			mu_product[j][i]       = R_Delta[j][i] * dmudrho_SRK[j][i];
			mu_Delta_approx2[j][i] = mu_product[j][i] + mu_Delta_cor[j][i];	
			lam_Delta_approx1[j][i] = lam_Delta_EoS[j][i] + lam_Delta_cor[j][i];
			lam_product[j][i]       = R_Delta[j][i] * dlamdrho_SRK[j][i];
			lam_Delta_approx2[j][i] = lam_product[j][i] + lam_Delta_cor[j][i];	
		} 
	} 

    // STEP 5: WRITE DEPENDENT VARIABLES
    // Initialize output streams
    std::ofstream result_compressibility;
    std::ofstream result_density_SRK;
    std::ofstream result_density_SW;
    std::ofstream result_specificHeat;
    std::ofstream result_lambda_SRK_Chung;
    std::ofstream result_lambda_SRK_Huber;
    std::ofstream result_lambda_SW_Chung;
    std::ofstream result_lambda_SW_Huber;
    std::ofstream result_mu_SRK_Chung;
    std::ofstream result_mu_SRK_LM;
    std::ofstream result_mu_SW_Chung;
    std::ofstream result_mu_SW_LM;
    std::ofstream result_dmudrho_SRK;
    std::ofstream result_dlamdrho_SRK;
    std::ofstream result_dmudrho_SW;
    std::ofstream result_dlamdrho_SW;
	std::ofstream result_R_Delta;
	std::ofstream result_mu_Delta;
	std::ofstream result_mu_Delta_EoS;
	std::ofstream result_mu_Delta_cor;
	std::ofstream result_lam_Delta;
	std::ofstream result_lam_Delta_EoS;
	std::ofstream result_lam_Delta_cor;
	std::ofstream result_mu_Delta_approx1;
	std::ofstream result_mu_product;
	std::ofstream result_mu_Delta_approx2;
	std::ofstream result_lam_Delta_approx1;
	std::ofstream result_lam_product;
	std::ofstream result_lam_Delta_approx2;

    // Create and open output text files
    result_compressibility.open   ("DV_compressibility.txt");
    result_density_SRK.open       ("DV_density_SRK.txt");
    result_density_SW.open        ("DV_density_SW.txt");
    result_lambda_SRK_Chung.open  ("DV_lambda_SRK_Chung.txt");
    result_lambda_SRK_Huber.open  ("DV_lambda_SRK_Huber.txt");
    result_lambda_SW_Chung.open   ("DV_lambda_SW_Chung.txt");
    result_lambda_SW_Huber.open   ("DV_lambda_SW_Huber.txt");
    result_mu_SRK_Chung.open      ("DV_mu_SRK_Chung.txt");
    result_mu_SRK_LM.open         ("DV_mu_SRK_LM.txt");
    result_mu_SW_Chung.open       ("DV_mu_SW_Chung.txt");
    result_mu_SW_LM.open          ("DV_mu_SW_LM.txt");
    result_dmudrho_SRK.open       ("DV_dmudrho_SRK.txt");
    result_dlamdrho_SRK.open      ("DV_dlamdrho_SRK.txt");
    result_dmudrho_SW.open        ("DV_dmudrho_SW.txt");
    result_dlamdrho_SW.open       ("DV_dlamdrho_SW.txt");
	result_R_Delta.open			  ("DV_R_Delta.txt");
	result_mu_Delta.open		  ("DV_mu_Delta.txt");
	result_mu_Delta_EoS.open	  ("DV_mu_Delta_EoS.txt");
	result_mu_Delta_cor.open	  ("DV_mu_Delta_cor.txt");
	result_lam_Delta.open		  ("DV_lam_Delta.txt");
	result_lam_Delta_EoS.open	  ("DV_lam_Delta_EoS.txt");
	result_lam_Delta_cor.open	  ("DV_lam_Delta_cor.txt");
	result_mu_Delta_approx1.open  ("DV_mu_Delta_approx1.txt");
	result_mu_product.open		  ("DV_mu_product.txt");
	result_mu_Delta_approx2.open  ("DV_mu_Delta_approx2.txt");
	result_lam_Delta_approx1.open ("DV_lam_Delta_approx1.txt");
	result_lam_product.open		  ("DV_lam_product.txt");
	result_lam_Delta_approx2.open ("DV_lam_Delta_approx2.txt");

    // Write title to output files
    result_compressibility   << substance << " Compressibility [dimensionless]\n";
    result_density_SRK       << substance << " Density_SRK [" << R_unit << "]\n";
    result_density_SW        << substance << " Density_SW [" << R_unit << "]\n";
    result_lambda_SRK_Chung  << substance << " Conductivity_SRK_Chung [" << lam_unit << "]\n";
    result_lambda_SRK_Huber  << substance << " Conductivity_SRK_Huber [" << lam_unit << "]\n";
    result_lambda_SW_Chung   << substance << " Conductivity_SW_Chung [" << lam_unit << "]\n";
    result_lambda_SW_Huber   << substance << " Conductivity_SW_Huber [" << lam_unit << "]\n";
    result_mu_SRK_Chung		 << substance << " Viscosity_SRK_Chung [" << mu_unit << "]\n";
    result_mu_SRK_LM		 << substance << " Viscosity_SRK_LM [" << mu_unit << "]\n";
    result_mu_SW_Chung		 << substance << " Viscosity_SW_Chung [" << mu_unit << "]\n";
    result_mu_SW_LM 		 << substance << " Viscosity_SW_LM [" << mu_unit << "]\n";
    result_dmudrho_SRK       << substance << " dmu/drho_SRK [" << mu_unit << "/" << R_unit << "]\n";
    result_dmudrho_SW        << substance << " dmu/drho_SW [" << mu_unit << "/" << R_unit << "]\n";
    result_dlamdrho_SRK      << substance << " dlam/drho_SRK [" << lam_unit << "/" << R_unit << "]\n";
    result_dlamdrho_SW       << substance << " dlam/drho_SW [" << lam_unit << "/" << R_unit << "]\n";
	result_R_Delta			 << substance << " Density_Delta [" << R_unit << "]\n"; ;
	result_mu_Delta			 << substance << " Viscosity_Delta [" << mu_unit << "]\n";
	result_mu_Delta_EoS		 << substance << " Viscosity_Delta_EoS [" << mu_unit << "]\n";
	result_mu_Delta_cor		 << substance << " Viscosity_Delta_corr [" << mu_unit << "]\n"; 
	result_lam_Delta		 << substance << " Conductivity_Delta [" << lam_unit << "]\n";
	result_lam_Delta_EoS	 << substance << " Conductivity_Delta_EoS [" << lam_unit << "]\n";
	result_lam_Delta_cor	 << substance << " Conductivity_Delta_corr [" << lam_unit << "]\n";
	result_mu_Delta_approx1	 << substance 
		<< " Viscosity_Delta_EoS + Viscosity_Delta_corr [" << mu_unit << "]\n";
	result_mu_product		 << substance 
		<< " Density_Delta * dmu/drho_SRK" << mu_unit << "]\n";
	result_mu_Delta_approx2	 << substance 
		<< " Density_Delta * dmu/drho_SRK + Viscosity_Delta_corr"  << mu_unit << "]\n";
	result_lam_Delta_approx1 << substance 
		<< " Conductivity_Delta_EoS + Conducitivity_Delta_corr [" << lam_unit << "]\n";
	result_lam_product		 << substance 
		<< " Density_Delta * dmu/drho_SRK" << lam_unit << "]\n";
	result_lam_Delta_approx2 << substance 
		<< " Density_Delta * dmu/drho_SRK + Viscosity_Delta_corr" << lam_unit << "]\n";

    // Loop over pressure and temperature. write results.
    for (int j = 0; j < P_steps; j++) {
        P = P_START + j*P_delta;
        for (int i = 0; i < T_steps; i++) {
            T = T_START + i*T_delta;

			// write values
		    result_compressibility   << Z[j][i];
            result_density_SRK       << R_SRK[j][i];   
            result_density_SW        << R_SW[j][i];   
            result_lambda_SRK_Chung  << lam_SRK_Chung[j][i];
            result_lambda_SRK_Huber  << lam_SRK_Huber[j][i];
            result_lambda_SW_Chung   << lam_SW_Chung[j][i];
            result_lambda_SW_Huber   << lam_SW_Huber[j][i];
            result_mu_SRK_Chung      << mu_SRK_Chung[j][i];
            result_mu_SRK_LM         << mu_SRK_LM[j][i];
            result_mu_SW_Chung       << mu_SW_Chung[j][i];
            result_mu_SW_LM          << mu_SW_LM[j][i];
			result_dmudrho_SRK		 << dmudrho_SRK[j][i];
			result_dlamdrho_SRK		 << dlamdrho_SRK[j][i];
			result_dmudrho_SW		 << dmudrho_SW[j][i];
			result_dlamdrho_SW		 << dlamdrho_SW[j][i];
			result_R_Delta 			 << R_Delta[j][i];
			result_mu_Delta 		 << mu_Delta[j][i];
			result_mu_Delta_EoS 	 << mu_Delta_EoS[j][i];
			result_mu_Delta_cor 	 << mu_Delta_cor[j][i];
			result_lam_Delta 		 << lam_Delta[j][i];
			result_lam_Delta_EoS 	 << lam_Delta_EoS[j][i];
			result_lam_Delta_cor 	 << lam_Delta_cor[j][i];
			result_mu_Delta_approx1  << mu_Delta_approx1[j][i];
			result_mu_product 		 << mu_product[j][i];
			result_mu_Delta_approx2  << mu_Delta_approx2[j][i];
			result_lam_Delta_approx1 << lam_Delta_approx1[j][i];
			result_lam_product 		 << lam_product[j][i];
			result_lam_Delta_approx2 << lam_Delta_approx2[j][i];
            if (i != T_steps-1) {
                result_compressibility   << ",";
                result_density_SRK       << ",";
                result_density_SW        << ",";
                result_lambda_SRK_Chung  << ",";
                result_lambda_SRK_Huber  << ",";
                result_lambda_SW_Chung   << ",";
                result_lambda_SW_Huber   << ",";
                result_mu_SRK_Chung      << ",";
                result_mu_SRK_LM         << ",";
                result_mu_SW_Chung       << ",";
                result_mu_SW_LM          << ",";
                result_dmudrho_SRK       << ",";
                result_dlamdrho_SRK      << ",";
                result_dmudrho_SW        << ",";
                result_dlamdrho_SW       << ",";
				result_R_Delta 			 << ",";
				result_mu_Delta 		 << ",";
				result_mu_Delta_EoS 	 << ",";
				result_mu_Delta_cor 	 << ",";
				result_lam_Delta 		 << ",";
				result_lam_Delta_EoS 	 << ",";
				result_lam_Delta_cor 	 << ",";
				result_mu_Delta_approx1  << ",";
				result_mu_product 		 << ",";
				result_mu_Delta_approx2  << ",";
				result_lam_Delta_approx1 << ",";
				result_lam_product 		 << ",";
				result_lam_Delta_approx2 << ",";
            }
        }	    
        result_compressibility   << std::endl;
        result_density_SRK       << std::endl;
        result_density_SW        << std::endl;
        result_lambda_SRK_Chung  << std::endl;
        result_lambda_SRK_Huber  << std::endl;
        result_lambda_SW_Chung   << std::endl;
        result_lambda_SW_Huber   << std::endl;
        result_mu_SRK_Chung      << std::endl;
        result_mu_SRK_LM         << std::endl;
        result_mu_SW_Chung       << std::endl;
        result_mu_SW_LM          << std::endl;
        result_dmudrho_SRK       << std::endl;
        result_dlamdrho_SRK      << std::endl;
        result_dmudrho_SW        << std::endl;
        result_dlamdrho_SW       << std::endl;
		result_R_Delta 			 << std::endl;
		result_mu_Delta 		 << std::endl;
		result_mu_Delta_EoS 	 << std::endl;
		result_mu_Delta_cor 	 << std::endl;
		result_lam_Delta 		 << std::endl;
		result_lam_Delta_EoS 	 << std::endl;
		result_lam_Delta_cor 	 << std::endl;
		result_mu_Delta_approx1  << std::endl;
		result_mu_product 		 << std::endl;
		result_mu_Delta_approx2  << std::endl;
		result_lam_Delta_approx1 << std::endl;
		result_lam_product 		 << std::endl;
		result_lam_Delta_approx2 << std::endl;
    }    
    result_compressibility.close();  
    result_density_SRK.close();
    result_density_SW.close();
    result_lambda_SRK_Chung.close();
    result_lambda_SRK_Huber.close();
    result_lambda_SW_Chung.close();
    result_lambda_SW_Huber.close();
    result_mu_SRK_Chung.close();
    result_mu_SRK_LM.close();
    result_mu_SW_Chung.close();
    result_mu_SW_LM.close();
    result_dmudrho_SRK.close();
    result_dlamdrho_SRK.close();
    result_dmudrho_SW.close();
    result_dlamdrho_SW.close();
	result_R_Delta.close();
	result_mu_Delta.close();
	result_mu_Delta_EoS.close();
	result_mu_Delta_cor.close();
	result_lam_Delta.close();
	result_lam_Delta_EoS.close();
	result_lam_Delta_cor.close();
	result_mu_Delta_approx1.close();
	result_mu_product.close();
	result_mu_Delta_approx2.close();
	result_lam_Delta_approx1.close();
	result_lam_product.close();
	result_lam_Delta_approx2.close();

/*
	// RECORD ERROR IN DENSITY
	// Store SRK Densities
	std::ifstream file1("DV_density_SRK.txt");
	double rho_SRK[P_steps][T_steps];
	std::istringstream input1[P_steps + 1];
	i = 0;
	for (std::string line; std::getline(file1,line,'\n');) {
		if (i == 0) {input1[i].str(line);}
		else {
			j = 0;
			input1[i].str(line);
			for (std::string value; std::getline(input1[i],value,',');) {
				rho_SRK[i-1][j] = std::stod(value);
				j++;
			}
		}
		i++;
	}

    // Store SW Densities
	std::ifstream file2("DV_density_SW.txt");
	double rho_SW[P_steps][T_steps];
	std::istringstream input2[P_steps + 1];
	i = 0;
	for (std::string line; std::getline(file2,line,'\n');) {
		if (i == 0) {input2[i].str(line);}
		else {
			j = 0;
			input2[i].str(line);
			for (std::string value; std::getline(input2[i],value,',');) {
				rho_SW[i-1][j] = std::stod(value);
				j++;
			}
		}
		i++;
	}

	// Store Difference in Densities
	double rho_Delta[P_steps][T_steps];
    std::ofstream result_density_Delta;
    result_density_Delta.open  ("DV_density_Delta.txt");
    for (int i = 0; i < P_steps; i++) {
        for (int j = 0; j < T_steps; j++) {
			rho_Delta[i][j]	= rho_SRK[i][j] - rho_SW[i][j];
		    result_density_Delta << rho_Delta[i][j];
            if (j != T_steps-1) {
                result_density_Delta << ",";
            }
        }	    
        result_density_Delta << std::endl;
    }    
    result_density_Delta.close();  


	// RECORD PRODUCT OF DENSITY ERROR AND PARTIAL DERIVATIVE
	// Store dlamdrho 
	std::ifstream file3("DV_dlamdrho_SRK.txt");
	double dlamdrho_SRK_array[P_steps][T_steps];
	std::istringstream input3[P_steps + 1];
	i = 0;
	for (std::string line; std::getline(file3,line,'\n');) {
		if (i == 0) {input3[i].str(line);}
		else {
			j = 0;
			input3[i].str(line);
			for (std::string value; std::getline(input3[i],value,',');) {
				dlamdrho_SRK_array[i-1][j] = std::stod(value);
				j++;
			}
		}
		i++;
	}

	// Store Product of dlamdrho and density error
    std::ofstream result_product;
    result_product.open  ("DV_product_dlamdrho_SRK_density_Delta.txt");
    for (int i = 0; i < P_steps; i++) {
        for (int j = 0; j < T_steps; j++) {
			
		    result_product << dlamdrho_SRK_array[i][j] * rho_Delta[i][j];
            if (j != T_steps-1) {
                result_product << ",";
            }
        }	    
        result_product << std::endl;
    }    
    result_product.close();  
*/


	amrex::Finalize();

    return 0;
}
