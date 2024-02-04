#include <iostream>
#include <vector>
#include <fstream>
#include <string>

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

    // Conversion Factors
    const amrex::Real atm2barye       = 1013250;
    const amrex::Real atm2pa          = 101325;
    const amrex::Real gcm32kgm3       = 1000; 
    const amrex::Real SIheat2cgsheat  = 10000;
    const amrex::Real SIcond2cgscond  = 10000;
    const amrex::Real SIvisc2cgsvisc  = 10;
	amrex::Real P_unit2cgs, mu_unit2cgs;

    // Declare Variables
    int i,j, T_steps, T_START, T_END, P_steps, P_START, P_END;
    double T_delta, P_delta, H2Ofrac;
    std::string substance, discard_param_name;
	std::string P_unit,mu_unit;
    amrex::Real P, T, Z, R, am, bm, Wbar, Cp, mu, xi, lam;
    amrex::Real output;
    amrex::Real Y[NUM_SPECIES];
    amrex::Real Ddiag[NUM_SPECIES]   = {0.0};
    amrex::Real chi_mix[NUM_SPECIES] = {0.0};
    const bool wtr_get_xi    = true;
    const bool wtr_get_mu    = true;
    const bool wtr_get_lam   = true;
    const bool wtr_get_Ddiag = false;
    const bool wtr_get_chi   = false;
   
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
		   >> discard_param_name >> mu_unit;
 
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
			  << "\n   P_unit:    \t" << "[" << P_unit  << "]"
			  << "\n   mu_unit:   \t" << "[" << mu_unit << "]"
              << std::endl << std::endl;
             
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

	// SETUP UNIT CONVERSIONS
	// pressure
	if (P_unit == "barye") {P_unit2cgs = 1;}
	else if (P_unit == "Pa") {P_unit2cgs = 10;}
	else if (P_unit == "MPa") {P_unit2cgs = 1e7;}
	else if (P_unit == "atm") {P_unit2cgs = 1013250;}
	else {P_unit2cgs = 1013250;}

	// viscosity
	if (mu_unit == "poise") {mu_unit2cgs = 1;}
	else if (mu_unit == "Pa_s") {mu_unit2cgs = 10;}
	else if (mu_unit == "mPa_s") {mu_unit2cgs = 0.01;}
	else if (mu_unit == "uPa_s") {mu_unit2cgs = 1e-5;}
	else {mu_unit2cgs = 10;}


	// SETUP FOR LOOP
	// Initialize output streams
	std::ofstream file;
    file.open  ("PelePhys_DV_result.txt");
    file << "H2O-CO2_mix Property [dimensionless] PelePhys\n";

	// Set pressure and temperature
    P = (P_START) * P_unit2cgs;
    T = T_START;
	double num_steps = 100.0;

	// LOOP OVER MASS FRACTION
	for (int i = 0; i < num_steps + 1; i++) {

		// Set H2O fraction
		H2Ofrac =  i / num_steps;

	    // Set mass fraction
	    for (int i = 0; i < NUM_SPECIES; i++) {
	   	     Y[i] = 0;
   		}
		Y[H2O_ID] = H2Ofrac;
		Y[CO2_ID] = 1 - H2Ofrac; 

	    // Set Wbar
	    srk.Y2WBAR(Y,Wbar);

	
    	// WRITE DEPENDENT VARIABLES
        // compute values
        srk.MixingRuleAmBm(T,Y,am,bm);
        srk.Calc_CompressFactor_Z(Z,am,bm,P,T,Wbar);
        srk.PYT2R(P,Y,T,R);            
        srk.RTY2Cp(R,T,Y,Cp);

        trans.transport(
            wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, wtr_get_chi, T,
            R, Y, Ddiag, chi_mix, mu, xi, lam, output, ltransparm);

        // write values
//        file << output;
		file << (i*1.0)/100;
        if (i != 100) {file << ",";} 
    }    
    file.close();  

    amrex::Finalize();

    return 0;
}
