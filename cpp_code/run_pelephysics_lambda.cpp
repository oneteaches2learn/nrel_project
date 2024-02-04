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
	amrex::Real P_unit2cgs, R_unit2cgs, mu_unit2cgs;

    // Declare Variables
    int i,j, T_steps, P_steps;
    double T_START, T_END, P_START, P_END, T_delta, P_delta;
    std::string substance, discard_param_name;
	std::string P_unit, R_unit, mu_unit;
    amrex::Real P, T, Z, R, am, bm, Wbar, Cp, mu, xi, lam;
    amrex::Real dmudrho;
	amrex::Real H2Ofrac;
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
		   >> discard_param_name >> R_unit
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
			  << "\n   R_unit:    \t" << "[" << R_unit  << "]"
			  << "\n   mu_unit:   \t" << "[" << mu_unit << "]" << std::endl
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
    // write temperature
    std::ofstream myfile;
    myfile.open("PelePhys_IV_temperature.txt");
    for (i = 0; i < T_steps; i++){
        T = T_START + i*T_delta;
        myfile << T; 
        if (i != T_steps-1) {myfile << ",";}
    }
    myfile.close();

    // write pressure
    myfile.open ("PelePhys_IV_pressure.txt");
    for (i = 0; i < P_steps; i++){
        P = P_START + i*P_delta; 
        myfile << P;
        if (i != P_steps-1) {myfile << ",";}
    }
    myfile.close();
 

    // WRITE DEPENDENT VARIABLES
    // Initialize output streams
    std::ofstream result_compressibility;
    std::ofstream result_density;
    std::ofstream result_specificHeat;
    std::ofstream result_conductivity;
    std::ofstream result_viscosity;
    std::ofstream result_dmudrho;

    // Create and open output text files
    result_compressibility.open  ("PelePhys_DV_compressibility.txt");
    result_density.open          ("PelePhys_DV_density.txt");
    result_specificHeat.open     ("PelePhys_DV_specificHeat.txt");
    result_conductivity.open     ("PelePhys_DV_conductivity.txt");
    result_viscosity.open        ("PelePhys_DV_viscosity.txt");
    result_dmudrho.open          ("PelePhys_DV_dmudrho.txt");

    // Write title to output files
    result_compressibility  << substance << " Compressibility [dimensionless] PelePhys\n";
    result_density          << substance << " Density [kg m^-3] PelePhys\n";
    result_specificHeat     << substance << " Specific Heat [J kg^-1 K^-1] PelePhys\n";
    result_conductivity     << substance << " Thermal Conductivity [W m^-1 K^-1] PelePhys\n";
    result_viscosity  		<< substance << " Viscosity [" << mu_unit << "] PelePhys\n";
    result_dmudrho          << substance << " dmu_drho [" << mu_unit << "/" << R_unit << "] PelePhys\n";

    // Loop over pressure and temperature. write results.
    for (int i = 0; i < P_steps; i++) {
        P = (P_START + i*P_delta) * P_unit2cgs;
        for (int j = 0; j < T_steps; j++) {
            T = T_START + j*T_delta;

            // compute values
            srk.MixingRuleAmBm(T,Y,am,bm);
            srk.Calc_CompressFactor_Z(Z,am,bm,P,T,Wbar);
            srk.PYT2R(P,Y,T,R);            
            srk.RTY2Cp(R,T,Y,Cp);

            trans.transport(
                wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, wtr_get_chi, T,
                R, Y, Ddiag, chi_mix, mu, xi, lam, dmudrho, ltransparm);

            // write values
		    result_compressibility  << Z;
            result_density          << R   * gcm32kgm3;
            result_specificHeat     << Cp  / SIheat2cgsheat;
            result_conductivity     << lam / SIcond2cgscond;
            result_viscosity 	    << mu  / mu_unit2cgs;
			result_dmudrho			<< dmudrho / mu_unit2cgs * R_unit2cgs;
            if (j != T_steps-1) {
                result_compressibility  << ",";
                result_density          << ",";
                result_specificHeat     << ",";
                result_conductivity     << ",";
                result_viscosity        << ",";
                result_dmudrho          << ",";
            }
        }	    
        result_compressibility  << std::endl;
        result_density          << std::endl;
        result_specificHeat     << std::endl;
        result_conductivity     << std::endl;
        result_viscosity        << std::endl;
        result_dmudrho          << std::endl;
    }    
    result_compressibility.close();  
    result_density.close();
    result_specificHeat.close();
    result_conductivity.close();
    result_viscosity.close();
    result_dmudrho.close();

    amrex::Finalize();

    return 0;
}
