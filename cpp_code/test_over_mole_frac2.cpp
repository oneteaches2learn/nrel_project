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
	// STEP 1: READ IN SPECIES PARAMETERS
	// Determine number of species
	int NUM_SPEC = 2; // Note: potentially update this to be automated later.
	
	// Read in species properties
	std::string substance[NUM_SPEC]; // substance name
	amrex::Real v_c[NUM_SPEC]; 		// critical volumes [cm3 / mol]
	amrex::Real T_c[NUM_SPEC];		// critical temperatures [K]
	amrex::Real sigma[NUM_SPEC];		// potential distance parameters [dimensionless]
	amrex::Real epsilon[NUM_SPEC];	// potential energy parameters [dimensionless]
	amrex::Real d[NUM_SPEC];			// dipoles [Debye]
	
	std::string discard_param_name;
    std::ifstream inputs("parameters.txt");
	for (int i = 0; i < NUM_SPEC; i++) {
		inputs >> discard_param_name >> substance[i]
			   >> discard_param_name >> v_c[i]
		   	   >> discard_param_name >> T_c[i]
		   	   >> discard_param_name >> sigma[i]
		   	   >> discard_param_name >> epsilon[i]
		   	   >> discard_param_name >> d[i];
	} 


	// STEP 2: SETUP INTERACTION MATRICES
	// Instantiate interaction matrices
	int MATRIX_SIZE = NUM_SPEC * NUM_SPEC;
	amrex::Real v_c_ij[MATRIX_SIZE];
	amrex::Real T_c_ij[MATRIX_SIZE];
	amrex::Real sigma_ij[MATRIX_SIZE];
	amrex::Real epsilon_ij[MATRIX_SIZE];
	amrex::Real d_ij[MATRIX_SIZE];

	// Compute interaction matrix entries
	for (int i = 0; i < NUM_SPEC; i++) {
		for (int j = 0; j < NUM_SPEC; j++) {
			v_c_ij[i + j]     = std::pow(std::pow(v_c[i],1.0/3.0) + std::pow(v_c[j],1.0/3.0),3.0) / 8.0;
			T_c_ij[i + j]     = std::sqrt(T_c[i] * T_c[j]);
			sigma_ij[i + j]   = std::sqrt(sigma[i] * sigma[j]);
			epsilon_ij[i + j] = std::sqrt(epsilon[i] * epsilon[j]);
			d_ij[i + j]       = d[i] * d[i] * d[j] * d[j]; 
		}
	}	


	// STEP 3: LOOP OVER MOLE FRACTION
	// Print state to command line
    std::cout << "\n Calculating effective reduced dipole for the mixture...";
	for (int i = 0; i < NUM_SPEC; i++) {
		std::cout << "\n   substance: \t" << substance[i]
				  << "\n   v_c:       \t" << v_c[i]
				  << "\n   T_c:       \t" << T_c[i]
				  << "\n   sigma:     \t" << sigma[i]
				  << "\n   epsilon:   \t" << epsilon[i]
				  << "\n   dipole:    \t" << d[i];
	}
	std::cout << std::endl
			  << "\n Results:"
			  << "\n   X[CO2]\tX[H2O]\t\td_rm";
//			  << "\t\tv_cm\t\tT_cm\t\tsigma_m\t\tepsilon_m\td_m";

	// Initialize output streams
	std::ofstream reduced_dipole;
    reduced_dipole.open  ("reduced_dipole.txt");

	// Begin loop
	int NUM_DIVISIONS = 100;
	for (int k = 0; k < NUM_DIVISIONS + 1; k++) {
		
		// Setup mole fractions	
		amrex::Real X[NUM_SPEC]; 		// mole fractions [dimensionless]
		X[1] = k * 1.0 / (NUM_DIVISIONS * 1.0);
		X[0] = 1 - X[1];


		// STEP 4: COMPUTE EFFECTIVE MIXTURE PARAMETERS
		// Instantiate mixture parameters
		amrex::Real v_cm = 0;
		amrex::Real T_cm = 0;
		amrex::Real sigma_m = 0;
		amrex::Real epsilon_m = 0;
		amrex::Real d_m = 0;

		// Compute effective mixture parameters
		for (int i = 0; i < NUM_SPEC; i++) {
			for (int j = 0; j < NUM_SPEC; j++) {
				v_cm += X[i] * X[j] * v_c_ij[i + j];
				T_cm += X[i] * X[j] * T_c_ij[i + j] * v_c_ij[i + j];
				sigma_m += X[i] * X[j] * std::pow(sigma_ij[i + j],3.0);
				epsilon_m += X[i] * X[j] * epsilon_ij[i+j] * std::pow(sigma_ij[i+j],3.0);
				d_m += X[i] * X[j] * d[i] * d[i] * d[j] * d[j] 
							/ (epsilon_ij[i + j] * std::pow(sigma_ij[i + j],3.0));
			}
		}
		T_cm = T_cm / v_cm;
		sigma_m = std::pow(sigma_m,1.0/3.0);
		epsilon_m = epsilon_m / std::pow(sigma_m,3.0);
		d_m = d_m * std::pow(sigma_m,3.0) * epsilon_m;			
		d_m = std::pow(d_m, 1.0 / 4.0);


		// STEP 5: COMPUTE EFFECTIVE REDUCED DIPOLE
		amrex::Real d_rm;
		d_rm = 131.3 * d_m / std::sqrt(v_cm * T_cm);


		// STEP 6: WRITE RESULTS
		reduced_dipole << d_rm;
        if (k != NUM_DIVISIONS) {reduced_dipole << ",";} 
		std::cout << "\n   " << X[0] 
                  << "\t\t"    << X[1] 
                  << "\t\t"    << d_rm;
//                << "\t\t"    << v_cm
//                << "\t\t"    << T_cm
//                << "\t\t"    << sigma_m
//                << "\t\t"    << epsilon_m
//                << "\t\t"    << d_m;
    }    
    reduced_dipole.close();  
	std::cout << std::endl << std::endl;
}				
