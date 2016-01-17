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





