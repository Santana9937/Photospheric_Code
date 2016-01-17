// **************************************************
// --------------- GRB Jet Parameters ---------------
// **************************************************

// Jet Isotropic equivalent luminosity in observer frame in ergs/sec
const double L = 1.0e51 ;

// Bulk Lorentz factor of jet (bulk_gamma)
const double bulk_gamma = 300.0 ;
// beta_jet (speed_jet/c) corresponding to bulk_gamma
const double beta_jet = sqrt( 1.0 - 1.0/(bulk_gamma*bulk_gamma) ) ;

// ***************************************************************************
// ------ Number of Photons, Electrons, Collected Photons for Spectrum -------
// -------------------- tau_initial and tau_photosphere ----------------------
// ***************************************************************************

// Number of photons in simulation
const int Nphotons = 1e6 ;

// Number of electrons in simulation
const int Nelectrons = 1e3 ;

// Number of photons collected for the spectrum
const int N_photon_collect = Nphotons/3 ;

// tau_initial defines the radius where all the photons are injected
const double tau_initial = 2.0 ;

// tau_photosphere defines the radius where photons escape
const double tau_photosphere = 1.0 ;

// **************************************************
// --------- Photon Distribution Properties ---------
// **************************************************

// Specify Photon Distribution Type:
// 1 = Blackbody Photons
// 2 = Power-Law Distribution of Photons
// If neither 1 or 2 is specified, default choice of BB photons executed
const int photon_distr_type = 1 ;

// For Blackbody photons, the input parameter is the photon temperature
// in the jet-comoving frame. Enter Tphoton in eV in the numerator
// and then temperature converted to Kelvin by dividing by kB = 8.6173e-5 eV/K
const double T_phot_comv = 1.0e2/8.6173e-5 ;

// To have peak energy at 300 keV with adiabatic cooling, keep this option
//const double T_phot_comv = ( 1.0e5*bulk_gamma*(1.0-beta_jet)*(pow( tau_initial/tau_photosphere , 2.0/3.0)) )/8.6173e-5 ;

// To have peak energy at 300 keV without adiabatic cooling, keep this option
//const double T_phot_comv = ( 1.7e5*bulk_gamma*(1.0-beta_jet) )/8.6173e-5 ;


// For power-law distribution photons,

// For photons from a Power-Law distribution, Power-Law begins at the
// energy E_1_phot_PL (jet-comv frame) and ends at the energy E_2_phot_PL (jet-comv frame)
// Give photon energies in eV and then energies converted to ergs with (1.6022e-12 ergs)/eV
// The photon index for PL is p_phot. f_nu propto nu^(1-p_phot)
const double E_1_phot_PL = 1.00e3*1.6022e-12 ;
const double E_2_phot_PL = 1.00e4*1.6022e-12 ;
const double p_phot = 1.5 ;

// **************************************************
// --------- Electron Distribution Properties -------
// **************************************************

// Specify Electron Distribution Type:
// 1 = mono-energetic electrons,
// 2 = Maxwell-Boltzman Distribution of Electrons
// 3 = Power-Law Distribution of Electrons
// If neither 1,2, or 3 specified, default choice of mono electrons executed
const int electron_distr_type = 2 ;

//  Mono-chromatic electrons, i.e. all electrons have same gamma_e
// Specify the initial electron Lorentz factor all the electrons have
const double gamma_e_init_mono = 2.0 ;

// Maxwell-Boltzmann distribution of electrons. Give gamma_e corresponding
// to electron temperature
const double gamma_e_init_MB = 2.0 ;

// Power-law distribution of electrons. Power-law begins
// at ge_1_PL and ends at ge_2_PL. The electron index is p_elec.
const double ge_1_PL =  2.0 ;
const double ge_2_PL = 30.0 ;
const double p_elec = 2.4 ;
const double E_1_PL = m_e*c*c*ge_1_PL, E_2_PL = m_e*c*c*ge_2_PL ;

// **************************************************
// --------------- Electron Re-heating --------------
// **************************************************

// Turning Electron Re-heating On/Off:
// 0 = OFF
// 1 = ON
// In neither 0 or 1 specified, default choice is Electron Re-heating OFF
const int reheating_knob = 0 ;

// Specify the number of electron re-heating events
const int N_tot_reheat_events = 10 ;

// NOTE: IF ELECTRON RE-HEATING TURNED ON, ELECTRONS
// RE-HEATED TO THE SAME ELECTRON DISTRIBUTION
// SPECIFIED IN ELECTRON DISTRIBUTION PROPERTIES

// Mono electron re-heating, electrons
// re-heated to gamma_e_reheat_mono
const double gamma_e_reheat_mono = gamma_e_init_mono ;

// MB electron re-heating, Give gamma_e corresponding
// to electron temperature
const double gamma_e_reheat_MB = gamma_e_init_MB ;

// Power-law distribution of electrons. Power-law begins
// at E_1_PL and ends at E_2_PL. The electron index is p.
// Give E_1_PL and E_2_PL in eV and will be converted to ergs
const double ge_1_PL_reheat =  ge_1_PL ;
const double ge_2_PL_reheat = ge_2_PL ;
const double E_1_PL_reheat = m_e*c*c*ge_1_PL_reheat, E_2_PL_reheat = m_e*c*c*ge_2_PL_reheat ;
const double p_elec_reheat = p_elec ;

// **************************************************
// ---------------- Adiabatic Cooling ---------------
// **************************************************

// Turning Adiabatic Cooling On/Off:
// 0 = OFF
// 1 = ON
// In neither 0 or 1 specified, default choice is Adiabatic cooling ON
const int adiab_cool_knob = 1 ;

// **************************************************
// ---------------- Scattering ON/OFF ---------------
// **************************************************

// Turning Scattering On/Off:
// 0 = OFF, Photon Energy and Electron gamma_e NOT updated due to IC/Comp Scattering
// 1 = ON, Photon Energy and Electron gamma_e updated due to IC/Comp Scattering
// In neither 0 or 1 specified, default choice is Scattering ON
const int scatt_cool_knob = 1 ;

// **************************************************
// ---------------- Electron Tracking ----------------
// **************************************************

// Electron Tracking On/Off for tracking 3 electrons:
// 0 = OFF, the gamma_e of 3 electrons throughout the simulation NOT tracked
// 1 = ON, the gamma_e of 3 electrons throughout the simulation tracked
// In neither 0 or 1 specified, default choice is Tracking OFF
const int elec_tracking_knob = 0 ;

// NOTE: If 3 electrons are tracked, the gamma_e of the electrons
// is updated after every electron re-heating episode, every adiabatic
// cooling episode, and every IC/Compton scattering episode

// Indicies of the electrons we are tracking
const int elec_track_1 = (Nelectrons-1) ;
const int elec_track_2 = (Nelectrons-1)/2 ;
const int elec_track_3 = (Nelectrons-1)/3 ;

// **************************************************
// ---- Bining The Photon And Electron Energies -----
// **************************************************

// Declaring the min and log10(energy in eV) considered for bining
// and the bin size in log10
const double min_log10_bin_ener_eV = -5.0 ;
const double max_log10_bin_ener_eV = 15.0 ;
const double bin_size_log10_eV = 0.1 ;

// **************************************************
// ---- Writing all E_ph and gamma_e at end of ------
// ---------- to file as previously done ------------
// **************************************************

// Saving All E_ph that escape photosphere and gamma_e to file:
// 0 = OFF, all photon and electron gamma_e not saved
// 1 = ON, all photon and electron gamma_e not saved
// Note: Recommended to Turn Off, for Nphotons = 1e8,
// file with photon energies has a size of about 0.3 GB
const int save_all_Eph_gam_e_knob = 0 ;





