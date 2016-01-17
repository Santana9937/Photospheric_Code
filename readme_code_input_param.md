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

### Photon Distribution

The code can handle to photon distributions: a Blackbody (BB) photon distrubution
and a Power-Law (PL) photon distribution.

- **photon_distr_type**: Choose the seed spectrum you want for the photons. Enter
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

### Electron Distribution

