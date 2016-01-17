

// This function returns a string with the simulation properties
// for the filename
string str_sim_properties_func(void)
{
  // Below, for each parameter, variable declaration for char and string.
  // Then, formating char with sprint and then converting char to string

  char char_T_phot_comv_eV[50] ;
  char char_E_1_phot_PL[50], char_E_2_phot_PL[50], char_p_phot[50] ;
  string str_T_phot_comv_eV ;
  string str_E_1_phot_PL, str_E_2_phot_PL, str_p_phot ;
  switch(photon_distr_type)
  {
    case 1: // Case 1 is BB photons
    {
      sprintf(char_T_phot_comv_eV, "%2.1e", T_phot_comv*8.6173e-5) ;
      str_T_phot_comv_eV = string(char_T_phot_comv_eV) ;
      break ;
    }
    case 2: // Case 2 is PL photons
    {
      sprintf(char_E_1_phot_PL, "%2.1e", E_1_phot_PL/(1.6022e-12)) ;
      sprintf(char_E_2_phot_PL, "%2.1e", E_2_phot_PL/(1.6022e-12)) ;
      sprintf(char_p_phot, "%.2f", p_phot ) ;
      str_E_1_phot_PL = string(char_E_1_phot_PL) ;
      str_E_2_phot_PL = string(char_E_2_phot_PL) ;
      str_p_phot = string(char_p_phot) ;
      break ;
    }
  }

  char char_bulk_gamma[50] ;
  string str_bulk_gamma ;
  sprintf(char_bulk_gamma, "%1i", int(bulk_gamma));
  str_bulk_gamma = string(char_bulk_gamma) ;


  char char_gamma_e_init_mono[50] ;
  char char_T_e_initial[50] ;
  char char_ge_1_pl[50], char_ge_2_pl[50], char_p_elec[50] ;
  string str_gamma_e_init_mono ;
  string str_T_e_initial ;
  string str_ge_1_pl, str_ge_2_pl, str_p_elec ;
  switch(electron_distr_type)
  {
    case 1: // Case 1 is mono electrons
    {
      sprintf(char_gamma_e_init_mono, "%1i", int(gamma_e_init_mono)) ;
      str_gamma_e_init_mono = string(char_gamma_e_init_mono) ;
      break ;
    }

    case 2: // Case 2 is MB electrons
    {
      sprintf(char_T_e_initial, "%2.1e", T_e_comv*8.6173e-5) ;
      str_T_e_initial = string(char_T_e_initial) ;
      break ;
    }

    case 3: // Case 3 is PL electrons
    {
      sprintf(char_ge_1_pl, "%1i", int(ge_1_PL)) ;
      sprintf(char_ge_2_pl, "%1i", int(ge_2_PL)) ;
      sprintf(char_p_elec, "%.2f", p_elec ) ;
      str_ge_1_pl = string(char_ge_1_pl) ;
      str_ge_2_pl = string(char_ge_2_pl) ;
      str_p_elec = string(char_p_elec) ;
      break ;
    }
  }

  char char_tau_initial[50] ;
  string str_tau_initial ;
  sprintf(char_tau_initial, "%1i", int(tau_initial));
  str_tau_initial = string(char_tau_initial) ;

  char char_Nphotons[50] ;
  string str_Nphotons ;
  sprintf(char_Nphotons, "%1.1e", double(Nphotons)) ;
  str_Nphotons = string(char_Nphotons) ;

  char char_Nelectrons[50] ;
  string str_Nelectrons ;
  sprintf(char_Nelectrons, "%1.1e", double(Nelectrons)) ;
  str_Nelectrons = string(char_Nelectrons) ;

  char char_N_photon_fraction[50] ;
  string str_N_photon_fraction ;
  sprintf(char_N_photon_fraction, "%1.2f", double(N_photon_collect)/Nphotons) ;
  str_N_photon_fraction = string(char_N_photon_fraction) ;


  // Putting all the simulation properties in a string. Declaring the string sim_param
  string sim_param ;

  // Photon distribution and photon distribution input parameters
  switch(photon_distr_type)
  {
    case 1: // Case 1 is BB photons
    {
      sim_param = "_Tph_" + str_T_phot_comv_eV ;
      break ;
    }
    case 2: // Case 2 is PL photons
    {
      sim_param = "_E1ph_" + str_E_1_phot_PL ;
      sim_param += "_E2ph_" + str_E_2_phot_PL ;
      sim_param += "_p_ph_" + str_p_phot ;
      break ;
    }
  }

  // Input value for bulk Lorentz factor
  sim_param += "_bG_" + str_bulk_gamma ;

  // Electron distribution and electron distribution input parameters
  switch(electron_distr_type)
  {
    case 1: // Case 1 is mono electrons
    {
      sim_param += "_mono_ge_" + str_gamma_e_init_mono ;
      break ;
    }

    case 2: // Case 2 is MB electrons
    {
      sim_param += "_MB_Te_" + str_T_e_initial ;
      break ;
    }

    case 3: // Case 3 is PL electrons
    {
      sim_param += "_PL_el_ge1_" + str_ge_1_pl ;
      sim_param += "_ge2_" + str_ge_2_pl ;
      sim_param += "_p_el_" + str_p_elec ;
      break ;
    }
  }

  sim_param += "_tau_" +  str_tau_initial ;
  sim_param += "_Nph_" + str_Nphotons ;
  sim_param += "_Ne_" + str_Nelectrons ;
  sim_param += "_frac_" + str_N_photon_fraction ;


  string sim_knobs ;


  char char_Nreheat[50] ;
  string str_Nreheat ;
  sprintf(char_Nreheat, "%1i", N_tot_reheat_events);
  str_Nreheat = string(char_Nreheat) ;

  switch(reheating_knob)
  {
    case 0: // Case 0 is electron re-heating OFF
    {
      sim_knobs = "_RH_no" ;
      break ;
    }
    case 1: // Case 1 is electron re-heating ON
    {
      sim_knobs = "_RH_yes_Nrh_" + str_Nreheat ;
      break ;
    }
  }

  switch(adiab_cool_knob)
  {
    case 0: // Case 0 is adiabatic cooling OFF
    {
      sim_knobs += "_AD_no" ;
      break ;
    }
    case 1: // Case 1 is adiabatic cooling ON
    {
      sim_knobs += "_AD_yes" ;
      break ;
    }
  }

  switch(scatt_cool_knob)
  {
    case 0: // Case 0 is IC/Comp Scattering OFF
    {
      sim_knobs += "_IC_no" ;
      break ;
    }
    case 1: // Case 1 is IC/Comp Scattering ON
    {
      sim_knobs += "_IC_yes" ;
      break ;
    }
  }

  // Declaring string with simulation properties
  string sim_properties = sim_param + sim_knobs ;

  return sim_properties ;

}


void fnu_bin_and_file_save(string sim_properties)
{

  // ***********************************************************************
  // **** Bining the data to get f_nu photon and f_nu electron spectrum ****
  // **** Saving f_nu photon, f_nu electron, and energies in eV to file ****
  // ***********************************************************************

  // Declaring the size of the vector that holds the bin edges and the vector itself
  const int size_log10_bin_edges_eV_vect = int( (max_log10_bin_ener_eV-min_log10_bin_ener_eV)/bin_size_log10_eV ) + 1 ;
  vector<double> vect_log10_bin_edges_eV(size_log10_bin_edges_eV_vect) ;

  // Filling in the vector that hold the edges of the bins
  for(int i=0 ; i<size_log10_bin_edges_eV_vect ; i++)
  {
    vect_log10_bin_edges_eV[i] = min_log10_bin_ener_eV + i*bin_size_log10_eV ;
    //cout << "vect_log10_bin_edges_eV[" << i << "]: " << vect_log10_bin_edges_eV[i] << endl ;
  }

  // Converting electron gamma_e to kinetic energy in the observer frame
  // Declaring vector that will hold electron KE in the observer frame
  vector<double> fin_elec_KE_obs_eV(Nelectrons) ;

  // Filling in the final kinetic energies of the electrons in the observer frame
  for(int i=0 ; i<Nelectrons ; i++)
  {
    fin_elec_KE_obs_eV[i] = 5.11e5*(gamma_e[i]-1.0)*bulk_gamma ;
  }

  // Now, declaring the vector dN, which holds the count for the number of photons/electrons in each bin
  // Subtracting 1 for the vector length since you lose 1 point when bining the data
  vector<int> dN_num_ph(size_log10_bin_edges_eV_vect-1) , dN_num_elec(size_log10_bin_edges_eV_vect-1) ;

  // Initializing all the entries of photon and electron vector to 0
  fill( dN_num_ph.begin() , dN_num_ph.end() , 0 ) ;
  fill( dN_num_elec.begin() , dN_num_elec.end() , 0 ) ;

  // Declaring an iterator to find the index of the bin where the photon/electron
  // belongs with lower_bound. Also, declaring an integer to store the array index
  vector<double>::iterator low ;
  int int_low ;

  // Finding the bin the photon belongs in and adding the photon to this bin
  for(int i=0 ; i<N_photon_collect ; i++)
  {
    low = lower_bound( vect_log10_bin_edges_eV.begin() , vect_log10_bin_edges_eV.end() , log10(fin_ph_ener_obs_eV[i]) ) ;
    int_low = low - vect_log10_bin_edges_eV.begin() -1 ;
    dN_num_ph[int_low] ++ ;
  }

  // Finding the bin the electron belongs in and adding the photon to this bin
  for(int i=0 ; i<Nelectrons ; i++)
  {
    low = lower_bound( vect_log10_bin_edges_eV.begin() , vect_log10_bin_edges_eV.end() , log10(fin_elec_KE_obs_eV[i]) ) ;
    int_low = low - vect_log10_bin_edges_eV.begin() -1 ;
    dN_num_elec[int_low] ++ ;
  }

  // Declaring the vectors that will hold the photon and electron spetrum
  vector<double> f_nu_ph(size_log10_bin_edges_eV_vect-1), f_nu_elec(size_log10_bin_edges_eV_vect-1) ;

  // Defining f_nu spectrum, which is E*(dN/dE) = dN/d ln E =
  // log_10 (exp(1)) * dN/d log_10 (E) . With bin size = d log_10 (E)
  for(unsigned int i=0 ; i< f_nu_ph.size() ; i++)
  {
    f_nu_ph[i] = (log10(exp(1.0)))*(dN_num_ph[i]/bin_size_log10_eV) ;
    f_nu_elec[i] = (log10(exp(1.0)))*(dN_num_elec[i]/bin_size_log10_eV) ;
  }

  // Now, finding the index of the max value in the vector to peak normalize
  vector<double>::iterator res_max_elem;
  int ind_of_max_fnuph, ind_of_max_fnuelec ;
  double max_val_fnuph, max_val_fnuelec ;

  res_max_elem = max_element(f_nu_ph.begin(), f_nu_ph.end()) ;
  ind_of_max_fnuph = res_max_elem - f_nu_ph.begin() ;
  max_val_fnuph = f_nu_ph[ind_of_max_fnuph] ;

  res_max_elem = max_element(f_nu_elec.begin(), f_nu_elec.end()) ;
  ind_of_max_fnuelec = res_max_elem - f_nu_elec.begin() ;
  max_val_fnuelec = f_nu_elec[ind_of_max_fnuelec] ;

  // Peak normalizing f_nu_ph and f_nu_elec
  for(unsigned int i=0 ; i< f_nu_ph.size() ; i++)
  {
    f_nu_ph[i] = f_nu_ph[i]/max_val_fnuph ;
    f_nu_elec[i] = f_nu_elec[i]/max_val_fnuelec ;
  }

  // Creating a vector that holds the energies in eV. This vector
  // will hold the energies at the midpoint of the bin edges
  vector<double> vect_ener_eV(size_log10_bin_edges_eV_vect-1) ;

  for(unsigned int i=0 ; i<vect_ener_eV.size() ; i++)
  {
    vect_ener_eV[i] = pow( 10.0 , min_log10_bin_ener_eV + bin_size_log10_eV/2.0 + i*(bin_size_log10_eV) ) ;
    //cout << "vect_ener_eV[" << i << "]: " << vect_ener_eV[i] << endl ;
  }

  string fnu_ph_filename = "sim_fnu_ph" + sim_properties + ".bin" ;
  string fnu_el_filename = "sim_fnu_el" + sim_properties + ".bin" ;
  string ener_eV_filename = "sim_en_eV" + sim_properties + ".bin" ;

  // Writting values of f_nu_ph, f_nu_elec, and vect_ener_eV
  ofstream fout_fnuph ;
  fout_fnuph.open(fnu_ph_filename, ios::binary | ios::out) ;
  for(unsigned int i=0; i < f_nu_ph.size() ; i++)
  {
    fout_fnuph.write((char*) &f_nu_ph[i], sizeof(f_nu_ph[i])) ;
  }
  fout_fnuph.close() ;

  ofstream fout_fnuelec ;
  fout_fnuelec.open(fnu_el_filename, ios::binary | ios::out) ;
  for(unsigned int i=0; i < f_nu_elec.size() ; i++)
  {
    fout_fnuelec.write((char*) &f_nu_elec[i], sizeof(f_nu_elec[i])) ;
  }
  fout_fnuelec.close() ;

  ofstream fout_eneV ;
  fout_eneV.open(ener_eV_filename, ios::binary | ios::out) ;
  for(unsigned int i=0; i < vect_ener_eV.size() ; i++)
  {
    fout_eneV.write((char*) &vect_ener_eV[i], sizeof(vect_ener_eV[i])) ;
  }
  fout_eneV.close() ;


  // ***********************************************************************
  // ******* Saving all E_ph in eV in observer frame to binary file ********
  // ***************** Saving all electron to binary file ******************
  // ***********************************************************************

  string ph_filename = "sim_ph_en" + sim_properties + ".bin" ;
  string elec_filename = "sim_el_ge" + sim_properties + ".bin" ;
  string log_filename = "sim_ph_en_log" + sim_properties + ".log" ;

  // -- Storing Photon energy in eV in electron gamma_e after simulation over --

  switch(save_all_Eph_gam_e_knob)
  {
    case 0: // Case 0 is Saving all E_ph and gamma_e OFF
    {
      break ;
    }
    case 1: // Case 0 is Saving all E_ph and gamma_e OFN
    {
      ofstream fout ;
      fout.open(ph_filename, ios::binary | ios::out) ;
      for(int i=0; i<N_photon_collect; i++)
      {
        fout.write((char*) &fin_ph_ener_obs_eV[i], sizeof(fin_ph_ener_obs_eV[i])) ;
      }
      fout.close() ;

      ofstream fout1 ;
      fout1.open(elec_filename, ios::binary | ios::out) ;
      for(int i=0; i<Nelectrons; i++)
      {
        fout1.write((char*) &gamma_e[i], sizeof(gamma_e[i])) ;
      }
      fout1.close() ;

      break ;
    }
  }


  // ***********************************************************************
  // ****** Storring gamma_e values of the electrons we are tracking *******
  // ***********************************************************************


  // Declaring strings for filename of electron tracking
  string gamma_1_filename = "sim_gam_e_1_track" + sim_properties + ".bin" ;
  string gamma_2_filename = "sim_gam_e_2_track" + sim_properties + ".bin" ;
  string gamma_3_filename = "sim_gam_e_3_track" + sim_properties + ".bin" ;

  // If electron-tracking is on, update the gamma_e of the 3 electrons we are tracking
  switch(elec_tracking_knob)
  {
    case 0: // Case 0 is electron tracking OFF
    {
      break ;
    }
    case 1: // Case 1 is electron tracking ON
    {
      ofstream fout3 ;
      fout3.open(gamma_1_filename, ios::binary | ios::out) ;
      for(unsigned int i=0; i<vec_gam_e_track_1.size(); i++)
      {
        fout3.write((char*) &vec_gam_e_track_1[i], sizeof(vec_gam_e_track_1[i])) ;
      }
      fout3.close() ;

      ofstream fout4 ;
      fout4.open(gamma_2_filename, ios::binary | ios::out) ;
      for(unsigned int i=0; i<vec_gam_e_track_2.size(); i++)
      {
        fout4.write((char*) &vec_gam_e_track_2[i], sizeof(vec_gam_e_track_2[i])) ;
      }
      fout4.close() ;

      ofstream fout5 ;
      fout5.open(gamma_3_filename, ios::binary | ios::out) ;
      for(unsigned int i=0; i<vec_gam_e_track_3.size(); i++)
      {
        fout5.write((char*) &vec_gam_e_track_3[i], sizeof(vec_gam_e_track_3[i])) ;
      }
      fout5.close() ;

      break ;
    }

  }


}




void logfile_func(string sim_properties, double fin_num_scatt_ph_esc_sum,
                  double N_scatt_each_ph_sum, double timeSec)
{

    string log_filename = "sim_log" + sim_properties + ".log" ;

    // Formating Nphotons and Nelectrons to print to screen in scientific notation
    char char_Nphot[50] ;
    string str_Nphot ;
    sprintf(char_Nphot, "%1.1e", double(Nphotons)) ;
    str_Nphot = string(char_Nphot) ;

    char char_Nelec[50] ;
    string str_Nelec ;
    sprintf(char_Nelec, "%1.1e", double(Nelectrons)) ;
    str_Nelec = string(char_Nelec) ;


    ofstream fout(log_filename) ;
    fout << "*Log file*" << endl ;
    fout << endl ;
    fout << "************** Simulation Parameters and Properties ************** " << endl ;
    fout << "------------------------------------------------------------------ " << endl ;
    fout << "- Isotropic Luminosity: L_iso = " << L << " ergs/sec" << endl ;
    fout << "- Photon number: Nphotons = " << str_Nphot << endl ;
    fout << "- Electron number: Nelectrons = " << str_Nelec << endl ;
    fout << "- Collected photon fraction: N_photon_fraction = " << double(N_photon_collect)/Nphotons << endl ;
    fout << "- Initial optical depth: tau_intial = " << tau_initial << endl ;
    fout << "------------------------------------------------------------------ " << endl ;

    switch(photon_distr_type)
    {
      case 1: // Case 1 is BB photons
      {
        fout << "- BB Photon Distr., jet-comv T_photon = " << T_phot_comv*8.6173e-5 << " eV" << endl ;
        break ;
      }

      case 2: // Case 2 is PL photons
      {
        fout << "- PL Photon Distr., E_1_phot_PL (jet-comv) = " << E_1_phot_PL/(1.6022e-12)
             << " eV," << endl << "  E_2_phot_PL (jet-comv) = " << E_2_phot_PL/(1.6022e-12)
             << " eV, p_phot = " << p_phot << endl ;
        break ;
      }
    }

    fout << "- Bulk gamma: Gamma = " << bulk_gamma << endl ;
    fout << "------------------------------------------------------------------ " << endl ;

    switch(electron_distr_type)
    {
      case 1: // Case 1 is mono electron
      {
        fout << "- Mono Electrons, Initial gamma_e: gamma_e_init_mono = " << gamma_e_init_mono << endl ;
        break ;
      }
      case 2:
      {
        fout << "- MB Electrons, jet-comv T_e_comv = " << T_e_comv*8.6173e-5
             << " eV, corresp. gamma_e_init_MB = " << gamma_e_init_MB << endl ;
        break ;
      }
      case 3:
      {
      fout << "PL Electron Distr., ge_1_PL = " << ge_1_PL
           << ", ge_2_PL = " << ge_2_PL << ", p_elec = " << p_elec << endl ;
      break ;
    }


    }

    fout << "------------------------------------------------------------------ " << endl ;

    switch(reheating_knob)
    {
      case 0: // Case 0 is electron re-heating OFF
      {
        fout << "- Electron Re-heating: OFF" << endl ;
        break ;
      }
      case 1: // Case 1 is electron re-heating ON
      {
       fout << "- Electron Re-heating: ON, N_reheat = " << N_tot_reheat_events << endl ;
       break ;
      }
    }


    switch(adiab_cool_knob)
    {
      case 0: // Case 0 is adiabatic cooling OFF
      {
        fout << "- Adiabatic Cooling: OFF" << endl ;
        break ;
      }
      case 1: // Case 1 is adiabatic cooling ON
      {
        fout << "- Adiabatic Cooling: ON" << endl ;
        break ;
      }
    }


    switch(scatt_cool_knob)
    {
      case 0: // Case 0 is IC/Comp Scattering OFF
      {
        fout << "IC/Comp Scattering: OFF" << endl ;
        break ;
      }
      case 1: // Case 1 is IC/Comp Scattering ON
      {
        fout << "IC/Comp Scattering: ON" << endl ;
        break ;
      }
    }


    switch(elec_tracking_knob)
    {
      case 0: // Case 0 is electron tracking OFF
      {
        fout << "Tracking 3 electrons: OFF" << endl ;
        break ;
      }
      case 1: // Case 1 is electron tracking ON
      {
        fout << "Tracking 3 electrons: ON" << endl ;
        break ;
      }
    }


    switch(save_all_Eph_gam_e_knob)
    {
      case 0: // Case 0 is Saving all E_ph and gamma_e OFF
      {
        fout << "Saving all E_ph that escape and final gamma_e values: OFF" << endl ;
        break ;
      }
      case 1: // Case 0 is Saving all E_ph and gamma_e OFN
      {
        fout << "Saving all E_ph that escape and final gamma_e values: ON" << endl ;
        break ;
      }
    }




    fout << "------------------------------------------------------------------" << endl ;
    fout << "- Ave. Num Scatt for Escape Phot: " << fin_num_scatt_ph_esc_sum/N_photon_collect << endl ;
    fout << "- Ave. Scatt All Photons: " << N_scatt_each_ph_sum/Nphotons << endl ;

    fout << "------------------------------------------------------------------" << endl ;
    fout << "- Code took: " << timeSec << " sec" << endl ;

    fout.close();

}


