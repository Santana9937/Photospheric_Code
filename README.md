# Photospheric Code 

This is a Monte Carlo code that simulates the Photospheric process, one of the main
mechanisms studied to explain the prompt gamma-ray emission from Gamma-ray Bursts. The
code is described in this paper: 
http://adsabs.harvard.edu/abs/2016MNRAS.456.1049S

Feel free to e-mail me at Santana9937@gmail.com if you have any questions.

### Requirements for running code and plotting code results

The code is written in C++11. A compiler capable of interpreting C++11 is 
needed to run the code. The code results are outputted as binary files. 
Any programming language can be used to plot the code results. In this repository, in
https://github.com/Santana9937/Photospheric_Code/tree/master/python_plotting_script
we include a simple Python 2.7 script to plot the code results. The Python libraries
Numpy and Matplotlib are needed to run this Python plotting script.

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

and hit enter to run the code (executable file). Running the code with 
the parameters in the repository should take under a minute.

### Output files of code

Once the code finishes running, 3 binary files and 1 log file are outputed into 
your working directory:

File beginning with sim_fnu_ph: Binary file that contains the f_nu photon spectrum
(peak normalized) at the end of the simulation in the observer frame.

File beginning with sim_fnu_el: Binary file that contains the f_nu electron spectrum
(peak normalized) at the end of the simulation in the observer frame.

File beginning with sim_en_eV: Binary file that contains the energy values in eV
corresponding to the f_nu photon spectrum and f_nu electron spectrum.

File beginning with sim_log: Log file containing the parameters used for the simulation.

### Plotting code results

In the terminal, cd in to the python_plotting_script directory. This directory 
contains the Python file py_phot_plot_script.py , used for plotting, and the 4 files
outputed at the end of the simulation. To run the plotting script, in the terminal type

python py_phot_plot_script.py

and hit enter. This python script saves the simulation results figure as a png and
eps file with the filenames fig_phot_code.png and fig_phot_code.eps, respectively.
The script then displays a figure with the simulation results. Once you are done 
looking at results, close the figure window to escape the Matplotlib figure environment. 


### Running the code with different parameters

Click on the markdown file 
https://github.com/Santana9937/Photospheric_Code/blob/master/readme_code_input_param.md
in this repository to see details on the code parameters and how to change 
the code parameters. 

### Description of code files and Simulation

Click on the markdown file
https://github.com/Santana9937/Photospheric_Code/blob/master/readme_descript_code_files.md    
to see a description on what each file of the code does and the design of the code
in the main function of the program.

