#include <CoolProp.h>
#include <iostream>
#include <fstream>

int main()
{
    // Conversion Factors
    int atm2pa      = 101325;    

    // Define Variables
    double T, T_START, T_END, T_delta, R, R_START, R_END, R_delta;
	double P_unit2SI, R_unit2SI, mu_unit2SI;
    int T_steps, R_steps, i, j, dump_counter, dump_reset;
    std::string substance, discard_param_name;
	std::string P_unit, R_unit, mu_unit;

    std::ofstream myfile;

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


    std::cout << "\n Running CoolProp with the following inputs..."
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
              << "\n   R_unit:    \t" << "[" << R_unit << "]"
			  << "\n   P_unit:    \t" << "[" << P_unit << "]"
			  << "\n   mu_unit:   \t" << "[" << mu_unit << "]"
              << std::endl << std::endl;
   
    // SETUP UNIT CONVERSIONS
    // pressure
    if (P_unit == "barye")    {P_unit2SI = 1e-1;}
    else if (P_unit == "Pa")  {P_unit2SI = 1;}
    else if (P_unit == "MPa") {P_unit2SI = 1e6;}
    else if (P_unit == "atm") {P_unit2SI = 101325;}
    else {P_unit2SI = 101325;} // presumed units of atm if P_unit not specified

    // density
    if (R_unit == "g/cm3")      {R_unit2SI = 1e3;}
    else if (R_unit == "kg/m3") {R_unit2SI = 1;}
    else if (R_unit == "g/mL")  {R_unit2SI = 1e3;}
    else if (R_unit == "g/L")   {R_unit2SI = 1;}
    else if (R_unit == "kg/L")  {R_unit2SI = 1e3;}
    else {R_unit2SI = 1;} // presumed SI units if R_unit not specified

    // viscosity
    if (mu_unit == "poise")      {mu_unit2SI = 1e-1;}
    else if (mu_unit == "Pa_s")  {mu_unit2SI = 1;}
    else if (mu_unit == "mPa_s") {mu_unit2SI = 1e-3;}
    else if (mu_unit == "uPa_s") {mu_unit2SI = 1e-6;}
    else {mu_unit2SI = 1;} // presumed SI units if mu_unit not specified


    // Setup progress bars
    dump_reset  = R_steps / 50;
    
    // Write Independent Variables
    myfile.open ("CoolProp_IV_temperature.txt");
    for (i = 0; i < T_steps; i++) {
	T = T_START + i * T_delta;
	myfile << T;
	if (i != T_steps-1) {
  	    myfile << ",";
	}
    }       
    myfile.close();    

    myfile.open("CoolProp_IV_density.txt");
    for (i = 0; i < R_steps; i++) {
	R = R_START + i * R_delta;
	myfile << R;
	if (i != R_steps-1) {
  	    myfile << ",";
	}
    }
    myfile.close();
 

    // Write Viscosity
    dump_counter = 0;
    std::cout << "\n Writing Viscosity" << std::endl;
    std::cout << "   v                                                v" << std::endl;
    std::cout << "   ";
    myfile.open ("CoolProp_DV_viscosity.txt");
    myfile << substance << " Viscosity [" << mu_unit << "] CoolProp\n";
    for (int i = 0; i < R_steps; i++) {
        R = (R_START + i * R_delta);
        for (int j = 0; j < T_steps; j++) {
            T = T_START + j * T_delta;
  	        myfile << CoolProp::PropsSI("viscosity","T",T,"Dmass",R * R_unit2SI,substance) 
						/ mu_unit2SI;
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


    // Write Pressure
    dump_counter = 0;
    std::cout << "\n Writing Pressure" << std::endl;
    std::cout << "   v                                                v" << std::endl;
    std::cout << "   ";
    myfile.open ("CoolProp_DV_pressure.txt");
    myfile << substance << " Pressure [" << P_unit << "] CoolProp\n";
    for (int i = 0; i < R_steps; i++) {
        R = (R_START + i * R_delta);
        for (int j = 0; j < T_steps; j++) {
            T = T_START + j*T_delta;
	    	myfile << CoolProp::PropsSI("P","T",T,"Dmass",R * R_unit2SI,substance)
						/ P_unit2SI;
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
    myfile.open ("CoolProp_DV_conductivity.txt");
    myfile << substance << " Thermal Conductivity [W m^-1 K^-1] CoolProp\n";
    for (int i = 0; i < R_steps; i++) {
        R = (R_START + i * R_delta);
        for (int j = 0; j < T_steps; j++) {
            T = T_START + j*T_delta;
		    myfile << CoolProp::PropsSI("conductivity","T",T,"Dmass",R * R_unit2SI,substance);
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
    myfile.open ("CoolProp_DV_specificHeat.txt");
    myfile << substance << " Specific Heat [J kg^-1 K^-1] CoolProp\n";
    for (int i = 0; i < R_steps; i++) {
        R = (R_START + i * R_delta);
        for (int j = 0; j < T_steps; j++) {
            T = T_START + j*T_delta;
		    myfile << CoolProp::PropsSI("Cpmass","T",T,"Dmass",R * R_unit2SI,substance);
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
    myfile.open ("CoolProp_DV_compressibility.txt");
    myfile << substance << " Compressibility [dimensionless] CoolProp\n";
    for (int i = 0; i < R_steps; i++) {
        R = (R_START + i * R_delta);
        for (int j = 0; j < T_steps; j++) {
            T = T_START + j*T_delta;
		    myfile << CoolProp::PropsSI("Z","T",T,"Dmass",R * R_unit2SI,substance);
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


    return EXIT_SUCCESS;
}
