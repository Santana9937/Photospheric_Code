# Description of code files and Simulation

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
This file contains the input parameters of the code. See
https://github.com/Santana9937/Photospheric_Code/blob/master/readme_code_input_param.md
for a description of this code file and the code parameters.

### fund_const.h
This file contains the fundamental in cgs units used in the code.

### electron_properties.h
This file initializes the directions of the electrons (drawn randomly) and
and the gamma_e and beta_e for the electrons, according to the electron
distribution specified in input_parameters.h. In addition, if the user has
decided to track 3 electrons (elec_tracking_knob=1), this file contains the
functions that track the gamma_e of the 3 electrons tracked after each 
scattering event. 

### photon_properties.h
This file initializes the directions of the photons (drawn randomly) and
the energy of the photons according to the distribution specified in 
input_parameters.h. This file also initializes the photons in a hemisphere
in the direction of the observer at a distance r_initial (corresponding
to tau_initial). In addition, the function that propagates the photons
after a s^prime is drawn for a photon is in this file (see Section 2.2
of the paper)

### electron_reheat_properties.h
This file contains the functions that reheat the elecrons if 
the user decides to include electron reheating in input_parameters.h.
For each reheating event, the directions of the electrons (drawn randomly) and
and the gamma_e and beta_e for the electrons, according to the electron
distribution specified in input_parameters

### adiab_cool_ph_and_el.h
This file contains the function that updates the photon energy and
beta_e, gamma_e of the electrons due to adiabatic cooling.

### IC_Compton_scatt.h
This file contains the functions that update the energy and direction
of the photons after a scattering event. This file also contains the
functions that update the beta_e, gamma_e and direction of the electrons
after a scattering event.

### filename_bin_save_logfile.h
This file converts the parameters specified in input_parameters.h into
strings. The strings are then concatenated to make the filenames for
the files outputed at the end of the simulation. This file also contains
the function that produces the log-file for the code. 

### func_proto_tuple_decl.h
This file contains the function protoypes of all the function in the code
and the declarations of the all the tuples used to store the 
values returned from the functions used in the code.

### print_sim_prop.h
This file contains the function that prints out all the simulation parameters
specified in input_parameters.h to the terminal at the start of the simulation.

### main_protospheric_code.h
This file contains the main function for the code. This function first checks
if the simulation was initialized with proper input parameter in input_parameters.h.
If not, the improper parameter is printed to the terminal and the code is
stopped. If the code has proper input parameter, the input parameter are printed
to the terminal and the simulation continues. The function the initializes the electrons
and photons for the simulation. Once all the photons are initialized, a message is
printed to the terminal. The code then does the first propagation for each
photon and builds the priority queue (see Section 2.4 of paper). Once the priority
queue is built, a message is printed and another message is printed, indicating that
the simulation has now begun. The simulation then proceeds until N_photon_collect
photons have escaped the photosphere. Once the simulation ends, the time the simulation
took is and the average number of scatterings for the photons are printed to the terminal.
The simulation results are then stored to binary files and the program ends.