# Code input parameter file

The file input_parameters.h contains all the input parameters for the the code. 
To run the code with different parameters, simply change the parameters
in this file, compile the code, and run the new executable file. Below is a
description of each of the parameters in this file.

### Gamma-ray Bursts (GRB) jet parameters

- **L**: Isotropic equivalent luminosity of the jet in the observer frame in ergs/sec.
- **bulk_gamma**: Bulk Lorentz factor of the GRB jet.

### Number of Photons, Electrons, and Collected Photons

- **Nphotons**: Number of photons in a simulation. An important note is that a 
simulation with 10^8 photons requires 9GB of RAM. For these simulations, a
computer with large memory or a computer cluster is needed.
- **Nelectrons**: Number of electrons in a simulation. An important note is that
at least 10^2 - 10^3 electrons are needed so that the output photon spectrum wont
be too noisy due to fluctuations in the Monte Carlo process.
- **N_photon_collect**: Number of photons collected for the output spectrum. See Section 2.1
of paper for more discussion on this parameter.

### Optical Depths

- **tau_initial**: Optical depth where all the photons are initilized.
- **tau_photosphere**: Optical depth of photosphere, where the photons escape.

### Initial Photon Distribution

The code can handle 2 photon distributions: a Blackbody (BB) photon distrubution
and a Power-Law (PL) photon distribution.

- **photon_distr_type**: Choose the seed photon spectrum you want for the photons. Enter
1 for BB distribution and enter 2 for PL distribution.
- **T_phot_comv**: Temperature of the seed photons in the jet-comoving frame.
In the numerator, enter the photon temperature in eV. The photon temperature
will then be converted to Kelvin by dividing by Boltzmann's constant 8.6173e-5 eV/K.
- **E_1_phot_PL**: Energy where the PL distribution of photons begins. Enter energy
in eV. Energy then converted to ergs with the conversion factor (1.6022e-12 ergs)/eV.
- **E_2_phot_PL**: Energy where the PL distribution of photons ends. Enter energy
in eV. Energy then converted to ergs with the conversion factor (1.6022e-12 ergs)/eV.
- **p_phot**: p_phot determines the slope of the seed f_nu photon spectrum. p_phot
is choosen so that f_nu propto nu^(1-p_phot). For example, if f_nu propto nu^-0.5 is
desired, enter p_phot = 1.5 .

### Initial Electron Distribution

The code can handle 3 electron distributions: Mono-energetic electrons (all
the electrons initialized to the same electron Lorentz factor), 
Maxwell-Boltzman (MB) Distribution of Electrons, and a Power-Law (PL) distribution
of electrons. 

- **electron_distr_type**: Choose the seed electron distribution you want 
for the electrons. Enter 1 for mono-energetic electrons, enter 2 for MB electrons,
and enter 3 for PL electrons.
- **gamma_e_init_mono**: Electron Lorentz factor all the electrons are initialized
to for mono-energetic electrons.
- **gamma_e_init_MB**: gamma_e corresponding to electron temperature for a MB distribution.
- **ge_1_PL**: Electron Lorentz factor where the power-law distribution of electrons begins.
- **ge_2_PL**: Electron Lorentz factor where the power-law distribution of electrons ends.
- **p_elec**: Power-law index for electron power-law distribution, 
i.e. dN_e/d_gamma_e propto gamma_e^-p

### Adiabatic Cooling

The user has the option deciding whether or not to include the adiabatic cooling
of photons and electrons as the jet expands outwards.

- **adiab_cool_knob**: Enter 1 to include adiabatic cooling. 
Enter 0 to NOT include adiabatic cooling

### Scattering ON/OFF

The user has the option deciding whether or not to include the exchange
of energy between photons and electrons due to IC/Compton scattering.

- **scatt_cool_knob**: Enter 1 to include the exchange of energy due to
IC/Compton scattering. Enter 0 to NOT include the exchange of energy 
due to IC/Compton scattering. 

### Electron Tracking

The user has the option deciding whether or not to track the gamma_e
of 3 electrons after every scattering. See Section 5.2.4 and Figure 6
of paper. 

- **elec_tracking_knob**: Enter 0 to NOT track the gamma_e of
3 electrons throughout the simulation. Enter 1 to track the gamma_e of
3 electrons throughout the simulation. Note: if 1 is entered for
elec_tracking_knob, 3 additional binary files are outputed by the code.

### Bining Photon and Electron Energies

The bining of the output spectrum of the photons and electrons 
is done in the code. 

- **min_log10_bin_ener_eV**: minimum value of log10(energy in eV) 
considered for bining
- **max_log10_bin_ener_eV**: maximum value of log10(energy in eV) 
considered for bining
- **bin_size_log10_eV**: bin size in log10(energy in eV) considered
for the bining

### Writing all final Photon energies and electron gamma_e to file

The user has the option deciding whether or not to save all the energies
of the N_photon_collect photons that escaped the photosphere at the end
of the simulation and the gamma_e of the electrons at the end of the 
simulation.

- **ssave_all_Eph_gam_e_knob**: Enter 0 to NOT save all energies
of photons that escaped the photosphere and all electron gamma_e at
the end of the simulation. Enter 1 to save all energies
of photons that escaped the photosphere and all electron gamma_e at
the end of the simulation. Note: if 1 is entered for
ssave_all_Eph_gam_e_knob, 2 additional binary files are outputed by the code.


