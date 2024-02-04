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

    // Declare Variables
    int i,j, T_steps, T_START, T_END, R_steps, R_START, R_END;
    double T_delta, R_delta;
    std::string substance, discard_param_name, P_unit, R_unit, mu_unit;
    amrex::Real P, T, Z, R, am, bm, Wbar, Cp, mu, xi, lam;
    amrex::Real dmudrho, P_unit2cgs, R_unit2cgs, mu_unit2cgs;
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
           >> discard_param_name >> R_delta
           >> discard_param_name >> R_steps
           >> discard_param_name >> R_START
		   >> discard_param_name >> R_END
		   >> discard_param_name >> P_unit
           >> discard_param_name >> R_unit
           >> discard_param_name >> mu_unit;


	// Complete missing data
	if (T_steps == 0) {T_steps = (T_END - T_START) / T_delta + 1;}
	if (T_END == 0)   {T_END   = T_START + T_steps * T_delta;}
	if (T_delta == 0) {T_delta = (T_END - T_START) / T_steps;}
	if (R_steps == 0) {R_steps = (R_END - R_START) / R_delta + 1;}
	if (R_END == 0)   {R_END   = R_START + R_steps * R_delta;} 
	if (R_delta == 0) {R_delta = (R_END - R_START) / R_steps;}

    // Print state to command line 
    std::cout << "\n Running PelePhysics with the following inputs..."
			  << "\n   substance: \t" << substance
              << "\n   T_delta:   \t" << T_delta
              << "\n   T_steps:   \t" << T_steps
              << "\n   T_start:   \t" << T_START
			  << "\n   T_end:     \t" << T_END
			  << "\n   T_unit:    \t" << "[K]"
              << "\n   R_delta:   \t" << R_delta
              << "\n   R_steps:   \t" << R_steps
              << "\n   R_start:   \t" << R_START 
			  << "\n   R_end:     \t" << R_END
			  << "\n   P_unit:    \t" << "[" << P_unit << "]"
			  << "\n   R_unit:    \t" << "[" << R_unit << "]"
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

    // Set mass fraction
    for (int i = 0; i < 16; i++) {
        Y[i] = 0;
    }
	if (substance == "CO2") {Y[CO2_ID] = 1;}
	if (substance == "CH4") {Y[CH4_ID] = 1;}
	if (substance == "O2")  {Y[O2_ID]  = 1;}
	if (substance == "H2O") {Y[H2O_ID] = 1;}
//	if (substance == "CH3O2H") {Y[CH3O2H_ID] = 1;}
//	if (substance == "OH")  {Y[OH_ID] = 1;}

    // Set Wbar
    srk.Y2WBAR(Y,Wbar);

    // SETUP UNIT CONVERSIONS
    // pressure
    if (P_unit == "barye") {P_unit2cgs = 1;}
    else if (P_unit == "Pa") {P_unit2cgs = 10;}
    else if (P_unit == "MPa") {P_unit2cgs = 1e7;}
    else if (P_unit == "atm") {P_unit2cgs = 1013250;}
    else {P_unit2cgs = 1013250;}

    // density
    if (R_unit == "g/cm3")      {R_unit2cgs = 1;}
    else if (R_unit == "kg/m3") {R_unit2cgs = 1e-3;}
    else if (R_unit == "g/mL")  {R_unit2cgs = 1;}
    else if (R_unit == "g/L")   {R_unit2cgs = 1e-3;}
    else if (R_unit == "kg/L")  {R_unit2cgs = 1;}
    else {R_unit2cgs = 1e-3;} // presumed SI units if R_unit not specified

    // viscosity
    if (mu_unit == "poise") {mu_unit2cgs = 1;}
    else if (mu_unit == "Pa_s") {mu_unit2cgs = 10;}
    else if (mu_unit == "mPa_s") {mu_unit2cgs = 0.01;}
    else if (mu_unit == "uPa_s") {mu_unit2cgs = 1e-5;}
    else {mu_unit2cgs = 10;}


    // WRITE INDEPENDENT VARIABLES
    // write temperatures
    std::ofstream myfile;
    myfile.open("PelePhys_IV_temperature.txt");
    for (i = 0; i < T_steps; i++){
        T = T_START + i*T_delta;
        myfile << T; 
        if (i != T_steps-1) {myfile << ",";}
    }
    myfile.close();

    // write densities
    myfile.open ("PelePhys_IV_density.txt");
    for (i = 0; i < R_steps; i++){
        R = R_START + i*R_delta; 
        myfile << R;
        if (i != R_steps-1) {myfile << ",";}
    }
    myfile.close();
 

    // WRITE DEPENDENT VARIABLES
    // Initialize output streams
    std::ofstream result_compressibility;
    std::ofstream result_pressure;
    std::ofstream result_specificHeat;
    std::ofstream result_conductivity;
    std::ofstream result_viscosity;
    std::ofstream result_dmudrho;

    // Create and open output text files
    result_compressibility.open  ("PelePhys_DV_compressibility.txt");
    result_pressure.open         ("PelePhys_DV_pressure.txt");
    result_specificHeat.open     ("PelePhys_DV_specificHeat.txt");
    result_conductivity.open     ("PelePhys_DV_conductivity.txt");
    result_viscosity.open 		 ("PelePhys_DV_viscosity.txt");
    result_dmudrho.open          ("PelePhys_DV_dmudrho(analytic).txt");


    // Write title to output files
    result_compressibility  << substance << " Compressibility [dimensionless] PelePhys\n";
    result_pressure         << substance << " Pressure [" << P_unit << "] PelePhys\n";
    result_specificHeat     << substance << " Specific Heat [J kg^-1 K^-1] PelePhys\n";
    result_conductivity     << substance << " Thermal Conductivity [W m^-1 K^-1] PelePhys\n";
    result_viscosity		<< substance << " Viscosity [" << mu_unit << "] PelePhys\n";
    result_dmudrho          << substance << 
             " dmudrho [(" <<  mu_unit << ") / (" << R_unit << ")] PelePhys\n";


    // Loop over pressure and temperature. write results.
    for (int i = 0; i < R_steps; i++) {
        R = R_START + i*R_delta;
        for (int j = 0; j < T_steps; j++) {
            T = T_START + j*T_delta;

            // compute values
            srk.MixingRuleAmBm(T,Y,am,bm);
			srk.RTY2P(R * R_unit2cgs,T,Y,P);
            srk.Calc_CompressFactor_Z(Z,am,bm,P,T,Wbar);
            srk.RTY2Cp(R * R_unit2cgs,T,Y,Cp);

            trans.transport(
                wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, wtr_get_chi, T,
                R * R_unit2cgs, Y, Ddiag, chi_mix, mu, xi, lam, dmudrho, ltransparm);

            // write values
		    result_compressibility  << Z;
            result_pressure         << P   / P_unit2cgs;
            result_specificHeat     << Cp  / SIheat2cgsheat;
            result_conductivity     << lam / SIcond2cgscond;
            result_viscosity        << mu  / mu_unit2cgs;
            result_dmudrho          << dmudrho / mu_unit2cgs * R_unit2cgs;
            if (j != T_steps-1) {
                result_compressibility  << ",";
                result_pressure         << ",";
                result_specificHeat     << ",";
                result_conductivity     << ",";
                result_viscosity        << ",";
                result_dmudrho          << ",";
            }
        }	    
        result_compressibility  << std::endl;
        result_pressure         << std::endl;
        result_specificHeat     << std::endl;
        result_conductivity     << std::endl;
        result_viscosity        << std::endl;
        result_dmudrho          << std::endl;
    }    
    result_compressibility.close();  
    result_pressure.close();
    result_specificHeat.close();
    result_conductivity.close();
    result_viscosity.close();
    result_dmudrho.close();

    amrex::Finalize();

    return 0;
}
