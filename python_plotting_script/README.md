# Python script for plotting simulation results

In this directory, we include a simple Python 2.7 script to plot 
the code results. The Python libraries Numpy and Matplotlib are
needed to run the plotting script.

This directory contains the Python file py_phot_plot_script.py used for 
plotting and the 4 files outputed at the end of the simulation (3 binary
files and 1 log-file). 

In the file file py_phot_plot_script.py , we first load the f_nu photon
spectrum, f_nu electron spectrum, and energy values in eV to arrays. 
The energy values are then converted from eV to keV. The f_nu photon
spectrum and the f_nu electron spectrum are then plotted in the same
figure, the figure is saved to a png file and eps file in your working
directory, and then the figure is shown in the computer. 

To run the plotting script, in the terminal type

python py_phot_plot_script.py

and hit enter. The simulation results figure is saved as a png and
eps file with the filenames fig_phot_code.png and fig_phot_code.eps, respectively.
The script then displays a window with the simulation results. Once done looking at
results, close the figure window to escape the Matplotlib figure environment. 

