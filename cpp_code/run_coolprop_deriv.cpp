#include <CoolProp.h>
#include <iostream>
#include <fstream>

int main()
{
    // Conversion Factors
    int atm2pa      = 101325;    

    // Define Variables
    float T, T_delta, P, P_delta;
    int T_steps, T_START, P_steps, P_START, i, j, dump_counter, dump_reset;
    std::string substance, discard_param_name;
    std::string derivative = "d(viscosity)/d(Dmass)|T";
    std::ofstream myfile;

    // Read Inputs from Config File
    std::ifstream inputs("parameters.txt");
    inputs >> discard_param_name >> substance
           >> discard_param_name >> T_delta
           >> discard_param_name >> T_steps
           >> discard_param_name >> T_START
           >> discard_param_name >> P_delta
           >> discard_param_name >> P_steps
           >> discard_param_name >> P_START
           >> discard_param_name >> derivative;
 
    std::cout << "\n Running CoolProp with the following inputs..."
              << "\n   substance: \t" << substance
              << "\n   T_delta:   \t" << T_delta
              << "\n   T_steps:   \t" << T_steps
              << "\n   T_start:   \t" << T_START
              << "\n   P_delta:   \t" << P_delta
              << "\n   P_steps:   \t" << P_steps
              << "\n   P_start:   \t" << P_START 
              << "\n   derivative:\t" << derivative
              << std::endl;
   
    // Setup progress bars
    dump_reset  = P_steps / 50;
    

    // Write Independent Variables
    myfile.open ("CoolProp_temperatures.txt");
    for (i = 0; i < T_steps; i++) {
	T = T_START + i*T_delta;
	myfile << T;
	if (i != T_steps-1) {
  	    myfile << ",";
	}
    }       
    myfile.close();    

    myfile.open("CoolProp_pressures.txt");
    for (i = 0; i < P_steps; i++) {
	P = P_START + i*P_delta;
	myfile << P;
	if (i != P_steps-1) {
  	    myfile << ",";
	}
    }
    myfile.close();
 

    // Write Derivative
    dump_counter = 0;
    std::cout << "\n Writing Derivative " << derivative << std::endl;
    std::cout << "   v                                                v" << std::endl;
    std::cout << "   ";
    myfile.open ("CoolProp_derivative.txt");
    myfile << substance << " " << derivative << " [units?] CoolProp\n";
    for (int i = 0; i < P_steps; i++) {
        P = (P_START + i*P_delta)*atm2pa;
        for (int j = 0; j < T_steps; j++) {
            T = T_START + j*T_delta;
  	        myfile <<  CoolProp::PropsSI(derivative,"T",T,"P",P,substance) * atm2pa;
	        if (j != T_steps-1) {
                myfile << ",";
			}
        }	    
        myfile << std::endl;
        dump_counter++;
        if (dump_counter == dump_reset) {
            std::cout << "|";
            std::cout.flush();
            dump_counter = 0;
        }        
    }    
    myfile.close();
    std::cout << std::endl << std::endl;

/*
    // Write Viscosity
    dump_counter = 0;
    std::cout << "\n Writing Viscosity" << std::endl;
    std::cout << "   v                                                v" << std::endl;
    std::cout << "   ";
    myfile.open ("CoolProp_viscosityDynamic.txt");
    myfile << "CO2 Dynamic Viscosity [Pa s] CoolProp\n";
    for (int i = 0; i < P_steps; i++) {
        P = (P_START + i*P_delta)*atm2pa;
        for (int j = 0; j < T_steps; j++) {
            T = T_START + j*T_delta;
  	        myfile << CoolProp::PropsSI("viscosity","T",T,"P",P,substance);
	        if (j != T_steps-1) {
                myfile << ",";
			}
        }	    
        myfile << std::endl;
        dump_counter++;
        if (dump_counter == dump_reset) {
            std::cout << "|";
            std::cout.flush();
            dump_counter = 0;
        }        
    }    
    myfile.close();
    std::cout << std::endl;


    // Write Density
    dump_counter = 0;
    std::cout << "\n Writing Density" << std::endl;
    std::cout << "   v                                                v" << std::endl;
    std::cout << "   ";
    // note: mass density used
    myfile.open ("CoolProp_density.txt");
    myfile << "CO2 Density [kg m^-3] CoolProp\n";
    for (int i = 0; i < P_steps; i++) {
        P = (P_START + i*P_delta)*atm2pa;
        for (int j = 0; j < T_steps; j++) {
            T = T_START + j*T_delta;
	    	myfile << CoolProp::PropsSI("Dmass","T",T,"P",P,substance);
	    	if (j != T_steps-1) {
            	    myfile << ",";
		    }
        }	    
        myfile << std::endl;
        dump_counter++;
        if (dump_counter == dump_reset) {
            std::cout << "|";
            std::cout.flush();
            dump_counter = 0;
        }        
    }    
    myfile.close(); 
    std::cout << std::endl;


    // Write Thermal Conductivity
    dump_counter = 0;
    std::cout << "\n Writing Conductivity" << std::endl;
    std::cout << "   v                                                v" << std::endl;
    std::cout << "   ";
    myfile.open ("CoolProp_conductivity.txt");
    myfile << "CO2 Thermal Conductivity [W m^-1 K^-1] CoolProp\n";
    for (int i = 0; i < P_steps; i++) {
        P = (P_START + i*P_delta)*atm2pa;
        for (int j = 0; j < T_steps; j++) {
            T = T_START + j*T_delta;
		    myfile << CoolProp::PropsSI("conductivity","T",T,"P",P,substance);
	    	if (j != T_steps-1) {
            	    myfile << ",";
		    }
        }	    
        myfile << std::endl;
        dump_counter++;
        if (dump_counter == dump_reset) {
            std::cout << "|";
            std::cout.flush();
            dump_counter = 0;
        }        
    }    
    myfile.close();
    std::cout << std::endl;


    // Write Specific Heat
    // note: mass specific constant pressure specific heat used
    dump_counter = 0;
    std::cout << "\n Writing Specific Heat" << std::endl;
    std::cout << "   v                                                v" << std::endl;
    std::cout << "   ";
    myfile.open ("CoolProp_specificHeat.txt");
    myfile << "CO2 Specific Heat [J kg^-1 K^-1] CoolProp\n";
    for (int i = 0; i < P_steps; i++) {
        P = (P_START + i*P_delta)*atm2pa;
        for (int j = 0; j < T_steps; j++) {
            T = T_START + j*T_delta;
		    myfile << CoolProp::PropsSI("Cpmass","T",T,"P",P,substance);
	    	if (j != T_steps-1) {
            	    myfile << ",";
		    }
        }	    
        myfile << std::endl;
        dump_counter++;
        if (dump_counter == dump_reset) {
            std::cout << "|";
            std::cout.flush();
            dump_counter = 0;
        }        
    }    
    myfile.close();
    std::cout << std::endl;


    // Write Compressibility
    dump_counter = 0;
    std::cout << "\n Writing Compressibility" << std::endl;
    std::cout << "   v                                                v" << std::endl;
    std::cout << "   ";
    myfile.open ("CoolProp_compressibility.txt");
    myfile << "CO2 Compressibility Factor [dimensionless] CoolProp\n";
    for (int i = 0; i < P_steps; i++) {
        P = (P_START + i*P_delta)*atm2pa;
        for (int j = 0; j < T_steps; j++) {
            T = T_START + j*T_delta;
		    myfile << CoolProp::PropsSI("Z","T",T,"P",P,substance);
	    	if (j != T_steps-1) {
	                myfile << ",";
		    }
        }	    
        myfile << std::endl;
        dump_counter++;
        if (dump_counter == dump_reset) {
            std::cout << "|";
            std::cout.flush();
            dump_counter = 0;
        }        
    }    
    myfile.close();
    std::cout << std::endl;
    std::cout << std::endl;
*/

    return EXIT_SUCCESS;
}
