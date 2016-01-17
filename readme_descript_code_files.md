# Description of code files

The files of the code are. 

- adiab_cool_ph_and_el.h
- electron_properties.h
- electron_reheat_properties.h
- filename_bin_save_logfile.h
- func_proto_tuple_decl.h
- fund_const.h
- IC_Compton_scatt.h
- input_parameters.h
- main_protospheric_code.cpp
- photon_properties.h
- print_sim_prop.h


A description of what each file does and the code design is given below.

### input_parameters.h
File contains the input parameters of the code. See
https://github.com/Santana9937/Photospheric_Code/blob/master/readme_code_input_param.md
for a description of this code file and the code parameters.

### fund_const.h
File contains the fundamental in cgs units used in the code.

### electron_properties.h
Variable declarations for electrons and void function which initializes
the direction and gamma_e of the electrons. There is another function in
this header file, which updates the gamma_e of the electrons we are tracking

### photon_properties.h
Variable declarations for photons and void function which initializes
the direction and energy of the photons. There is another function which
calculates the new location of the photon and the elapsed time of the photon
(both in the observer frame) after a s_disp drawn in the jet-comoving frame.

### electron_reheat_properties.h
This header file contains the function which draws new directions and gamma_e
for the electrons after a re-heating episode. The variable declarations for
electron re-heating are also in this header file

### adiab_cool_ph_and_el.h
// Contains a function which updates photon energy and electron
// gam_e, beta_e due to adiabatic cooling

### IC_Compton_scatt.h
This header file contains 3 function related to the IC/Compton scattering.
The first function calculates the photon energy and direction after scattering
and calculates an acceptance rejection quantity. The second function checks to
see if the energy and direction of the photon after scattering are accepted. The
last function updates the energy and direction of the electrons after scattering

### filename_bin_save_logfile.h
Function which converts simulation parameters to strings for file-saving.
A function which bins the photon and electron energies and finds the f_nu
spectrum is also included. Thus function also writes the spectrum to file.
Lastly, this header file also contains a function which creates a log-file
which saves the simulation properties.

### func_proto_tuple_decl.h
This header-file includes the function prototypes of all the functions in the program
Tuples to extract values from function are also declared in this file

### print_sim_prop.h
This header file includes a function which prints the simulation parameters
and properties to the command prompt

### main_protospheric_code.h
xxx
