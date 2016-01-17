/* ------ Declaration of Tuples That Store Values Returned From Functions ------ */

// Tuple for Photon Propagation Function: phot_propag_func
tuple<double, double, double, double> phot_propag_step ;

// Tuple for returning photon energy and electron gam_e, beta_e after adiabatic cooling
tuple<double, double, double> abiab_cool_ph_el_return ;

// Tuple for funct_scatt_acc_rej declared in IC_Compton_scatt.h

// Tuple for returning photon properties after scattering from funct_IC_Comp_scattering
tuple<double, double, double, double> phot_IC_scatt_return ;

// Tuple for returning electron properties after scattering from update_elec_func
tuple<double, double, double, double, double> elec_IC_scatt_return ;


/* ---------- Function Prototypes before the main function ----------- */

// This void function initializes electrons
void initialize_electrons() ;

// This void function initializes photons
void initialize_photons() ;

// This function propagates photons forward in observer frame and computes
// time elapsed in the observer frame for this particular s_disp_pr drawn
tuple<double, double,
      double, double> phot_propag_func(double beta_jet, double bulk_gamma,
                                       double r_ph_pos_x_init, double r_ph_pos_y_init,
                                       double r_ph_pos_z_init,
                                       double s_disp_pr_fin, double om_one_ph_fin,
                                       double om_two_ph_fin, double om_three_ph_fin) ;

// This void function fills in the electron gamma_e, beta_e, and electron direction
// for each re-heating episode from the distribution specified
void reheat_electrons(void) ;

// Function prototype for adiabatic cooling of photons and electrons
tuple<double, double, double> func_adiab_cool_ph_el(double t_init_this_ph_obs_fr,
                                                    double t_delta_ph, double t_elec_init,
                                                    double E_ph_init, double gam_e_init) ;

// This function in scatt_acc_rej_func.cpp returns energy and direction of photon
// after scattering accepted and also
tuple<double, double, double,
      double, double> funct_scatt_acc_rej(double rnum_mu_pr, double rphi_mu_pr,
                                          double rnum_ener_acc_rej,
                                          double beta_e_init, double gamma_e_init,
                                          double v10elec_init, double v20elec_init,
                                          double v30elec_init, double E_phot_comv_init,
                                          double x_dim_ph_ener, double om_one_ph_init,
                                          double om_two_ph_init, double om_three_ph_init) ;

// This function in IC_Compton_scatt.cpp returns energy and direction of photon after scattering accepted
tuple<double, double,
      double, double> funct_IC_Comp_scattering(double beta_e_init,
                                               double gamma_e_init, double E_phot_comv_init,
                                               double v10elec_init, double v20elec_init,
                                               double v30elec_init, double om_one_ph_init,
                                               double om_two_ph_init, double om_three_ph_init) ;

// This function in upd_elec_aft_scatt.cpp updates electron energy and direction after scattering
tuple<double, double, double,
      double, double> update_elec_func(double gamma_e_init, double beta_e_init,
                                       double v10elec_init, double v20elec_init,
                                       double v30elec_init, double E_phot_comv_init,
                                       double E_phot_comv_fin, double om_one_ph_init,
                                       double om_two_ph_init, double om_three_ph_init,
                                       double om_one_ph_fin, double om_two_ph_fin,
                                       double om_three_ph_fin) ;

// Function for file saving
string str_sim_properties_func(void) ;

void fnu_bin_and_file_save(string sim_properties) ;

// Function prototype for logfile
void logfile_func(string sim_properties, double fin_num_scatt_ph_esc_sum,
                  double N_scatt_each_ph_sum, double timeSec) ;




bool print_sim_param_and_props(void) ;



void elec_track_switch_push(int curr_elec_ind) ;





