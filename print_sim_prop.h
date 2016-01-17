

bool print_sim_param_and_props(void)
{

  // Formating Nphotons and Nelectrons to print to screen in scientific notation
  char char_Nphot[50] ;
  string str_Nphot ;
  sprintf(char_Nphot, "%1.1e", double(Nphotons)) ;
  str_Nphot = string(char_Nphot) ;

  char char_Nelec[50] ;
  string str_Nelec ;
  sprintf(char_Nelec, "%1.1e", double(Nelectrons)) ;
  str_Nelec = string(char_Nelec) ;


  // Printing input parameters with no knobs
  cout << endl ;
  cout << "Simulation Parameters: " << endl ;
  cout << "Nphotons = " << str_Nphot << ", Nelectrons = " << str_Nelec << ", bulk_gamma = " << bulk_gamma << endl ;
  cout << "tau_initial = " << tau_initial <<  ", Fraction of Collected Photons = " << double(N_photon_collect)/double(Nphotons) << endl ;

  // *******************************************************
  // ******* Printing Photon Distribution Properties *******
  // *******************************************************

  switch(photon_distr_type)
  {
    case 1: // Case 1 is BB photons
    {
      cout << "BB Photon Distr., jet-comv T_photon = " << T_phot_comv*8.6173e-5 << " eV" << endl ;
      break ;
    }
    case 2: // Case 2 is PL photons
    {
      cout << "PL Photon Distr., E_1_phot_PL (jet-comv) = " << E_1_phot_PL/(1.6022e-12)
           << " eV, E_2_phot_PL (jet-comv) = " << E_2_phot_PL/(1.6022e-12)
           << " eV, p_phot = " << p_phot << endl ;
      break ;
    }
  }

  if( (photon_distr_type!=1) && (photon_distr_type!=2))
  {
   cout << endl << "UNACCEPTABLE INTEGER ENTERED FOR: photon_distr_type" << endl ;
   cout << "Enter 1 for Blackbody Photons or 2 for Power-Law Photons" << endl << endl ;
   // Enter return to exit out of the main function and terminate code
   return false ;
  }

  // *******************************************************
  // ****** Printing Electron Distribution Properties ******
  // *******************************************************

  switch(electron_distr_type)
  {
    case 1: // Case 1 is mono electrons
    {
      cout << "Mono-Electron Distr., gamma_e_initial = " << gamma_e_init_mono << endl ;
      break ;
    }

    case 2: // Case 2 is MB electrons
    {
      cout << "MB Electron Distr., jet-comv T_electron = " << T_e_comv*8.6173e-5
           << " eV, corresp. gamma_e_init_MB = " << gamma_e_init_MB << endl ;
      break ;
    }

    case 3: // Case 2 is PL electrons
    {
      cout << "PL Electron Distr., ge_1_PL = " << ge_1_PL
           << ", ge_2_PL = " << ge_2_PL << ", p_elec = " << p_elec << endl ;
      break ;
    }
  }

  if( (electron_distr_type!=1) && (electron_distr_type!=2) && (electron_distr_type!=3))
  {
   cout << endl << "UNACCEPTABLE INTEGER ENTERED FOR: electron_distr_type" << endl ;
   cout << "Enter 1 for mono electrons, 2 for MB electrons, or 3 for PL distribution electron" << endl << endl ;
   // Enter return to exit out of the main function and terminate code
   return false ;
  }

  // *******************************************************
  // ******* Printing Electron Re-Heating Properties *******
  // *******************************************************

  switch(reheating_knob)
  {
    case 0: // Case 0 is electron re-heating OFF
    {
     cout << "Electron Re-heating: OFF" << endl ;
     break ;
    }
    case 1: // Case 1 is electron re-heating ON
    {
     cout << "Electron Re-heating: ON, N_reheat = " << N_tot_reheat_events << endl ;
     break ;
    }
  }

  if( (reheating_knob!=0) && (reheating_knob!=1) )
  {
   cout << endl << "UNACCEPTABLE INTEGER ENTERED FOR: reheating_knob" << endl ;
   cout << "Enter 0 for electron re-heating OFF, Enter 1 for electron re-heating ON" << endl << endl ;
   // Enter return to exit out of the main function and terminate code
   return false ;
  }

  // *******************************************************
  // ******** Printing Adiabatic Cooling Properties ********
  // *******************************************************

  switch(adiab_cool_knob)
  {
    case 0: // Case 0 is adiabatic cooling OFF
    {
     cout << "Adiabatic Cooling: OFF" << endl ;
     break ;
    }
    case 1: // Case 1 is adiabatic cooling ON
    {
     cout << "Adiabatic Cooling: ON" << endl ;
     break ;
    }
  }

  if( (adiab_cool_knob!=0) && (adiab_cool_knob!=1) )
  {
   cout << endl << "UNACCEPTABLE INTEGER ENTERED FOR: adiab_cool_knob" << endl ;
   cout << "Enter 0 for adiabatic cooling OFF, Enter 1 for adiabatic cooling ON" << endl << endl ;
   // Enter return to exit out of the main function and terminate code
   return false ;
  }

  // *******************************************************
  // ******** Printing IC/Comp Scattering Properties *******
  // *******************************************************

  switch(scatt_cool_knob)
  {
    case 0: // Case 0 is IC/Comp Scattering OFF
    {
      cout << "IC/Comp Scattering: OFF" << endl ;
      break ;
    }
    case 1: // Case 1 is IC/Comp Scattering ON
    {
      cout << "IC/Comp Scattering: ON" << endl ;
      break ;
    }
  }

  if( (scatt_cool_knob!=0) && (scatt_cool_knob!=1) )
  {
   cout << endl << "UNACCEPTABLE INTEGER ENTERED FOR: scatt_cool_knob" << endl ;
   cout << "Enter 0 for IC/Comp Scattering OFF, Enter 1 for IC/Comp Scattering ON" << endl << endl ;
   // Enter return to exit out of the main function and terminate code
   return false ;
  }

  // *******************************************************
  // ******** Printing Electron Tracking Properties ********
  // *******************************************************

  switch(elec_tracking_knob)
  {
    case 0: // Case 0 is electron tracking OFF
    {
      cout << "Tracking 3 electrons: OFF" << endl ;
      break ;
    }
    case 1: // Case 1 is electron tracking ON
    {
      cout << "Tracking 3 electrons: ON" << endl ;
      break ;
    }
  }


  if( (elec_tracking_knob!=0) && (elec_tracking_knob!=1) )
  {
   cout << endl << "UNACCEPTABLE INTEGER ENTERED FOR: elec_tracking_knob" << endl ;
   cout << "Enter 0 for Electron Tracking OFF, Enter 1 for Electron Tracking ON" << endl << endl ;
   // Enter return to exit out of the main function and terminate code
   return false ;
  }


  // *******************************************************
  // ******** Printing Storring all E_ph and gamma_e *******
  // *******************************************************

  switch(save_all_Eph_gam_e_knob)
  {
    case 0: // Case 0 is Saving all E_ph and gamma_e OFF
    {
      cout << "Saving all E_ph that escape and final gamma_e values: OFF" << endl ;
      break ;
    }
    case 1: // Case 0 is Saving all E_ph and gamma_e OFN
    {
      cout << "Saving all E_ph that escape and final gamma_e values: ON" << endl ;
      break ;
    }
  }


  if( (save_all_Eph_gam_e_knob!=0) && (save_all_Eph_gam_e_knob!=1) )
  {
   cout << endl << "UNACCEPTABLE INTEGER ENTERED FOR: save_all_Eph_gam_e_knob" << endl ;
   cout << "Enter 0 for E_ph and gamma_e saving OFF, Enter 1 for E_ph and gamma_e saving ON" << endl << endl ;
   // Enter return to exit out of the main function and terminate code
   return false ;
  }

  cout << endl ;


  return true ;

}



