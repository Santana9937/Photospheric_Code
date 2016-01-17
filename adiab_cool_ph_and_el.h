

tuple<double, double, double> func_adiab_cool_ph_el(double t_init_this_ph_obs_fr,
                                                    double t_delta_ph, double t_elec_init,
                                                    double E_ph_init, double gam_e_init)
{
  // Declaring variables returned by function
  double E_ph_fin, gam_e_fin, beta_e_fin ;

  // Final photon energy with adiabatic cooling
  E_ph_fin = E_ph_init*(pow( (r_initial + (t_init_this_ph_obs_fr + t_delta_ph)*c*beta_jet)/
                             (r_initial + (t_init_this_ph_obs_fr)*c*beta_jet) , -2.0/3.0 )) ;

  // Declaring exponent for adiabatic cooling for electron
  double elec_adiab_index = 2.0 - 2.0*( (4.0*gam_e_init + 1.0)/(3.0*gam_e_init) ) ;

  // Final gam_e of electron with adiabatic cooling
  gam_e_fin = 1.0 + (gam_e_init-1.0)*(pow( (r_initial + (t_init_this_ph_obs_fr + t_delta_ph)*c*beta_jet)/
                                           (r_initial + t_elec_init*c*beta_jet), elec_adiab_index )) ;

  // Calculating final beta_e from final gam_e
  beta_e_fin = sqrt(1.0 - 1.0/(gam_e_fin*gam_e_fin)) ;

  return make_tuple(E_ph_fin, gam_e_fin, beta_e_fin) ;

}
