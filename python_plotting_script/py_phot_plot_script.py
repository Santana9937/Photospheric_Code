import matplotlib.pyplot as plt
import numpy as np

# ---------- Loading code results from C++ into Python ------------

# Peak normalized photon spectrum in eV in observer frame
f_nu_photon = np.fromfile("sim_fnu_ph_Tph_1.0e+02_bG_300_MB_Te_2.6e+05_tau_2_Nph_1.0e+06_Ne_1.0e+03_frac_0.33_RH_no_AD_yes_IC_yes.bin")

# Peak normalized electron spectrum in eV in observer frame
f_nu_electron = np.fromfile("sim_fnu_el_Tph_1.0e+02_bG_300_MB_Te_2.6e+05_tau_2_Nph_1.0e+06_Ne_1.0e+03_frac_0.33_RH_no_AD_yes_IC_yes.bin")

# Midpoint of energy bins in eV in observer frame
energy_eV_bin_mid_pt = np.fromfile("sim_en_eV_Tph_1.0e+02_bG_300_MB_Te_2.6e+05_tau_2_Nph_1.0e+06_Ne_1.0e+03_frac_0.33_RH_no_AD_yes_IC_yes.bin")

# Converting energies from eV to keV
energy_keV_bin_mid_pt = energy_eV_bin_mid_pt/1000.0

# ---------- Plotting code results in a loglog f_nu spectrum -------------

plt.figure(1)
#
plt.loglog(energy_keV_bin_mid_pt, f_nu_photon, '-', linewidth=3, color="green", label = r"$f_{\nu}$" + " photon")
#
plt.hold(True)
#
plt.loglog(energy_keV_bin_mid_pt, f_nu_electron, '--', linewidth=3, color="blue", label = r"$f_{\nu}$" + " electron")
#
plt.hold(False)
#
plt.xlabel(r"$E $ (keV)", fontsize=22)
plt.ylabel(r"$f_{\nu}$" + " (peak normalized)", fontsize=22)
plt.minorticks_on()
plt.axis([1.0e0, 1.0e6, 1.0e-5, 1.0e1])
plt.legend(loc="upper right", fontsize=14)

# Save figure to a png and eps file and then show figure
plt.savefig("fig_phot_code.png", bbox_inches="tight")
plt.savefig("fig_phot_code.eps", bbox_inches="tight")
plt.show()

