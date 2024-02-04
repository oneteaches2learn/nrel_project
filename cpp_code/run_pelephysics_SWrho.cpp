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
    int i, j, T_steps, P_steps;
    double T_delta, P_delta, T_START, T_END, P_START, P_END;
    std::string substance, discard_param_name;
    std::string P_unit, R_unit, mu_unit;
    amrex::Real P, Z, am, bm, Wbar, Cp, mu, xi, lam;
    amrex::Real dmudrho, dlamdrho;
    amrex::Real Y[NUM_SPECIES];
    amrex::Real Ddiag[NUM_SPECIES]   = {0.0};
    amrex::Real chi_mix[NUM_SPECIES] = {0.0};
    const bool wtr_get_xi    = true;
    const bool wtr_get_mu    = true;
    const bool wtr_get_lam   = true;
    const bool wtr_get_Ddiag = false;
    const bool wtr_get_chi   = false;
  
	double H2Ofrac;

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
	if (substance == "CH4") {Y[CH4_ID] = 1;}
	if (substance == "O2")  {Y[O2_ID]  = 1;}
	if (substance == "H2O") {Y[H2O_ID] = 1;}
//	if (substance == "HF")  {Y[HF_ID]  = 1;}
//	if (substance == "H2")  {Y[H2_ID]  = 1;}
//	if (substance == "N2")	{Y[N2_ID]  = 1;}

//	Y[CO2_ID] = 1 - H2Ofrac;
//	Y[H2O_ID] = H2Ofrac;	

    // Set Wbar
    srk.Y2WBAR(Y,Wbar);

    // WRITE INDEPENDENT VARIABLES
    // store temperature and write temperature to file
	double T[T_steps];
    std::ofstream myfile;
    myfile.open("PelePhys_IV_temperature.txt");
    for (i = 0; i < T_steps; i++){
        T[i] = T_START + i*T_delta;
        myfile << T[i]; 
        if (i != T_steps-1) {myfile << ",";}
    }
    myfile.close();

    // write pressures
    myfile.open ("PelePhys_IV_pressure.txt");
    for (i = 0; i < P_steps; i++){
        P = P_START + i*P_delta; 
        myfile << P;
        if (i != P_steps-1) {myfile << ",";}
    }
    myfile.close();

	// Store Densities
	std::ifstream file("CoolProp_DV_density.txt");
	double R[P_steps][T_steps];
	std::istringstream input[P_steps];
	i = 0;
	for (std::string line; std::getline(file,line,'\n');) {
		j = 0;
		input[i].str(line);
		for (std::string value; std::getline(input[i],value,',');) {
			R[i][j] = std::stod(value);
			j++;
		}
		i++;
	}

    // Print state to command line 
    std::cout << "\n Computing dmudrho from CoolProp densities for the following grid..."
			  << "\n   substance: \t" << substance
              << "\n   T_delta:   \t" << T_delta
              << "\n   T_steps:   \t" << T_steps
              << "\n   T_start:   \t" << T_START
			  << "\n   T_end:     \t" << T_END
              << "\n   P_delta:   \t" << P_delta
              << "\n   P_steps:   \t" << P_steps
              << "\n   P_start:   \t" << P_START 
			  << "\n   P_end:     \t" << P_END
              << "\n   P_unit:    \t" << "[" << P_unit  << "]"
              << "\n   R_unit:    \t" << "[" << R_unit  << "]"
              << "\n   mu_unit:   \t" << "[" << mu_unit << "]" << std::endl
			  << "\n   Y[CH4]:    \t" << Y[CH4_ID]
			  << "\n   Y[CO2]:    \t" << Y[CO2_ID]
			  << "\n   Y[H2O]:    \t" << Y[H2O_ID]
			  << "\n   Y[O2]:     \t" << Y[O2_ID]
//			  << "\n   Y[N2]:     \t" << Y[N2_ID]
//			  << "\n   Y[H2]:     \t" << Y[H2_ID]
              << std::endl << std::endl;
 

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


    // WRITE DEPENDENT VARIABLES
    // Initialize output streams
    std::ofstream result_dmudrho;
	std::ofstream result_viscosity;

    // Create and open output text files
    result_dmudrho.open          ("PelePhys_DV_dmudrho_rho(SW).txt");
    result_viscosity.open        ("PelePhys_DV_viscosity_rho(SW).txt");


    // Write title to output files
    result_dmudrho          << substance << 
            " dmudrho [(" << mu_unit << ") / (" << R_unit << ")] PelePhys\n";
    result_viscosity        << substance << " Viscosity [" << mu_unit << "] PelePhys\n";


    // Loop over pressure and temperature. write results.
    for (int i = 0; i < P_steps; i++) {
        for (int j = 0; j < T_steps; j++) {

            // compute values
//			srk.MixingRuleAmBm(T,Y,am,bm);
//			srk.RTY2P(R[i][j] * R_unit2cgs,T,Y,P * P_unit2cgs);
//			srk.Calc_CompressFactor_Z(Z,am,bm,P * P_unit2cgs,T,Wbar);
//			srk.RTY2Cp(R[i][j] * R_unit2cgs,T,Y,Cp);

            trans.transport(
                wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, wtr_get_chi, T[j],
                R[i][j] * R_unit2cgs, Y, Ddiag, chi_mix, mu, xi, lam, dmudrho, dlamdrho, ltransparm);

            // write values
			result_dmudrho   << dmudrho / mu_unit2cgs * R_unit2cgs;
			result_viscosity << mu / mu_unit2cgs;
            
			if (j != T_steps-1) {
                result_dmudrho          << ",";
				result_viscosity		<< ",";
            }
        }	    
        result_dmudrho          << std::endl;
		result_viscosity		<< std::endl;
    }    
    result_dmudrho.close();
	result_viscosity.close();

    amrex::Finalize();

    return 0;
}
