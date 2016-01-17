# Photospheric Code 

This is a Monte Carlo code that simulates the Photospheric process, one of the main
mechanisms studied to explain the prompt gamma-ray emission from Gamma-ray Bursts. The
code is described in this paper: 
http://adsabs.harvard.edu/abs/2016MNRAS.456.1049S

### Requirements for running code and plotting code results

The code is written in C++11. A compiler capable of interpreting C++11 is 
needed to run the code. The code results are outputed as binary files. 
Any programming language can be used to plot the code results. In this repository,
in the file python_plotting_script, we include a simple Python 2.7 script to plot 
the code results.

### Downloading the code

This repository can be cloned with the command:

git clone https://github.com/Santana9937/Photospheric_Code.git

Or, you can click Download ZIP in the top right to download this repository.

### Commands for running the code with GCC compiler

Open a terminal and cd into the directory where you downloaded the code files.
To compile the code, type

g++ -Wall -std=c++11 -O3 main_protospheric_code.cpp -o main_protospheric_code

and hit enter. Once the code is successfully compiled, an executable file named
main_protospheric_code is created. Compiling the code should not take 
longer than 5-10 seconds. We also note that we use O3 optimization flag to 
compile the code. Then type

./main_protospheric_code

and hit enter to run the code. Running the code with the parameters in the repository
should take under a minute.

### Output files of code

Once the code finishes running, 3 binary files and 1 log file are outputed into 
your working directory:

File beginning with sim_fnu_ph: Binary file that contains the f_nu photon spectrum
(peak normalized) in the observer frame.

File beginning with sim_fnu_ph: Binary file that contains the f_nu photon spectrum
(peak normalized) at the end of the simulation in the observer frame.

File beginning with sim_fnu_el: Binary file that contains the f_nu electron spectrum
(peak normalized) at the end of the simulation in the observer frame.

File beginning with sim_en_eV: Binary file that contains the energy values in eV
corresponding to the f_nu photon spectrum and f_nu electron spectrum.

File beginning with sim_log: Log file containing the parameters used for the simulation.

### Plotting code results

In the terminal, cd in to the python_plotting_script directory.

