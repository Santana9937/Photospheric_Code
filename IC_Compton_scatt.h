
// Function which calculates photon energy and direction after scattering and
// calculates quantity for acceptance rejection
tuple<double, double, double,
      double, double> funct_scatt_acc_rej(double rnum_mu_pr, double rphi_mu_pr,
                                          double rnum_ener_acc_rej,
                                          double beta_e_init, double gamma_e_init,
                                          double v10elec_init, double v20elec_init,
                                          double v30elec_init, double E_phot_comv_init,
                                          double x_dim_ph_ener, double om_one_ph_init,
                                          double om_two_ph_init, double om_three_ph_init)
{

  double mu_pr, phi_pr, rho_elec ;
  double om_one_ph_fin, om_two_ph_fin, om_three_ph_fin ;
  double scat_angle, x_pr_over_x, E_phot_comv_fin, capitalx, phot_ener_acc_rej ;

  // phi_prime and mu_prime of photon after scattering
  mu_pr = ( (beta_e_init + 2.0*rnum_mu_pr -1.0)/
             (1.0 + beta_e_init*(2.0*rnum_mu_pr -1.0)) ) ;
  phi_pr = 2.0*PI*rphi_mu_pr ;

  // Defining rho_elec; used to express photon direction after scattering in cartesian coordinates
  rho_elec = sqrt( v10elec_init*v10elec_init + v20elec_init*v20elec_init ) ;

  // om_one_ph_pr, om_one_ph_pr, and om_one_ph_pr are photon directions
  // after scattering in cartesian coordinates in the jet comoving frame
  om_one_ph_fin = ( mu_pr*v10elec_init + (sqrt(1-mu_pr*mu_pr)/rho_elec)*
                   (v20elec_init*cos(phi_pr) + v10elec_init*v30elec_init*sin(phi_pr)) ) ;

  om_two_ph_fin = ( mu_pr*v20elec_init + (sqrt(1-mu_pr*mu_pr)/rho_elec)*
                   (-v10elec_init*cos(phi_pr) + v20elec_init*v30elec_init*sin(phi_pr)) ) ;

  // Decrease om_one_ph_fin,om_two_ph_fin if round-off causes normalization condition to be imaginary
  if( om_one_ph_fin*om_one_ph_fin + om_two_ph_fin*om_two_ph_fin >= 1.0 )
  {
    om_one_ph_fin = 0.9999*om_one_ph_fin ;
    om_two_ph_fin = 0.9999*om_two_ph_fin ;
  }

  om_three_ph_fin = copysign(sqrt( 1.0 - om_one_ph_fin*om_one_ph_fin - om_two_ph_fin*om_two_ph_fin) ,
                             mu_pr*v30elec_init - sqrt(1-mu_pr*mu_pr)*rho_elec*sin(phi_pr) ) ;

  // Scattering angle; dot product of intial photon direction and photon direction after scattering
  scat_angle = ( om_one_ph_init*om_one_ph_fin +
                 om_two_ph_init*om_two_ph_fin +
                 om_three_ph_init*om_three_ph_fin ) ;

  // x_prime_over_x; ratio of dimensionless photon energy in the electron rest frame after scattering
  // to dimensionless photon energy in the electron rest frame before scattering
  x_pr_over_x = 1.0/( 1.0 +
                     (E_phot_comv_init/(gamma_e_init*m_e*c*c))*
                     ((1.0-scat_angle)/(1.0 - mu_pr*beta_e_init)) ) ;

  // Photon energy after scattering in the jet comoving frame
  E_phot_comv_fin = ( (x_pr_over_x*x_dim_ph_ener*m_e*c*c)/
                      (2.0*gamma_e_init*(1.0 - mu_pr*beta_e_init)) ) ;

  // Quantity given in Pozdnyakov et al 1983 Pg. 325 to determine if
  // direction of photon after scattering is accepted
  capitalx = (  1.0/x_pr_over_x + x_pr_over_x +
                (4.0/x_dim_ph_ener)*(1.0 - 1.0/x_pr_over_x) +
                (4.0/(x_dim_ph_ener*x_dim_ph_ener))*(1.0 - 1.0/x_pr_over_x)*(1.0 - 1.0/x_pr_over_x)  ) ;

  // Condition to determine if outgoing photon direction is accepted or rejected
  phot_ener_acc_rej = 2.0*rnum_ener_acc_rej - x_pr_over_x*x_pr_over_x*capitalx ;

  return make_tuple(phot_ener_acc_rej, E_phot_comv_fin,
                    om_one_ph_fin, om_two_ph_fin, om_three_ph_fin) ;

}


// Function which returns photon energy and direction after energy and direction of photon
// after scattering is accepted. Function in scatt_acc_rej_func.h for acceptance-rejection
// condition called here.
tuple<double, double,
      double, double> funct_IC_Comp_scattering(double beta_e_init,
                                               double gamma_e_init, double E_phot_comv_init,
                                               double v10elec_init, double v20elec_init,
                                               double v30elec_init, double om_one_ph_init,
                                               double om_two_ph_init, double om_three_ph_init)
{

  // Variable declarations for quantities used in this function
  double mu_pr_rnum, phi_pr_rnum, scatt_prob_acc_rej_rnum, ener_acc_rej_rnum ;
  double sigma_hat , scatt_acc_rej_cond ;
  double phot_ener_acc_rej, E_phot_comv_fin, om_one_ph_fin, om_two_ph_fin, om_three_ph_fin ;

  // Declaring tuple to extract photon energy and direction and acceptance-rejection condition
  tuple<double, double, double, double, double> ener_acc_rej_res ;

  /* --------- Declaring Random Number Generator Information -------- */
  // declare a random number generator
  // here we choose the Mersenne Twister engine (32bit)
  mt19937 rng;

  random_device randomSeed; // this uses hardware parameters (clocks etc) to generate a random seed
  rng.seed(randomSeed()); // set the seed from the random device

  /* ---------------------------------------------------------------- */

  // Drawing random numbers; 2 for photon direction after scattering,
  // 1 to determine if scattering will happen, i.e. compare random number to sigma/sigma_T
  // 1 for acceptance-rejection for photon energy and direction after scattering
  uniform_real_distribution<double> rn_mu_dir(0.0,1.0) ;
  uniform_real_distribution<double> rn_phi_dir(0.0,1.0) ;
  uniform_real_distribution<double> rn_scatt_prob_acc_rej(0.0,1.0) ;
  uniform_real_distribution<double> rn_ener_acc_rej_dir(0.0,1.0) ;

  // mu is the dot product of initial electron and initial photon direction
  double mu = ( v10elec_init*om_one_ph_init +
                v20elec_init*om_two_ph_init +
                v30elec_init*om_three_ph_init ) ;

  // Dimensionless photon energy in electron rest frame
  double x_dim_ph_ener = ( 2.0*gamma_e_init*E_phot_comv_init*(1.0-mu*beta_e_init) )/(m_e*c*c) ;

  // Calculating dimensionless cross-section (sigma_hat) for electron photon interaction
  // Just as a note, x_dim_ph_ener is always x_dim_ph_ener>0 (positive)
  if (x_dim_ph_ener<=0.5)
  {
      sigma_hat = (  1.0/3.0 + 0.141*x_dim_ph_ener - 0.12*x_dim_ph_ener*x_dim_ph_ener +
                  (1.0 + 0.5*x_dim_ph_ener)/((1.0 + x_dim_ph_ener)*(1.0 + x_dim_ph_ener))  ) ;
  }
  else if ( (0.5<x_dim_ph_ener) && (x_dim_ph_ener<=3.5) )
  {
      sigma_hat = (log(1.0 + x_dim_ph_ener) + 0.06)/x_dim_ph_ener ;
  }
  else
  {
      sigma_hat = ( (log(1.0 + x_dim_ph_ener) + 0.5 -
                      1.0/(2.0 + 0.076*x_dim_ph_ener) )/x_dim_ph_ener ) ;
  }

  // Draw random number to determine if scattering event will happen
  scatt_prob_acc_rej_rnum = rn_scatt_prob_acc_rej(rng) ;

  // Scattering event acceptance-rejection, comparing random number to ratio
  // of sigma/sigma_T
  scatt_acc_rej_cond  = scatt_prob_acc_rej_rnum - (2.0*PI*r_e*r_e*sigma_hat)/sigma_T ;

  // If scatt_acc_rej_cond<0, scattering event will happen
  if (scatt_acc_rej_cond<0)
  {
    do
    {
      // Drawing and assigning 3 random numbers for photon scattering. 2 for photon direction
      // after scattering and 1 for acceptance-rejection condition for photon direction
      mu_pr_rnum = rn_mu_dir(rng), phi_pr_rnum = rn_phi_dir(rng), ener_acc_rej_rnum = rn_ener_acc_rej_dir(rng) ;

      // funct_scatt_acc_rej calculates energy and direction of photon after scattering and
      // phot_ener_acc_rej. which is used to determine if scattering is accepted or rejected
      ener_acc_rej_res = funct_scatt_acc_rej(mu_pr_rnum, phi_pr_rnum,
                                             ener_acc_rej_rnum,
                                             beta_e_init, gamma_e_init,
                                             v10elec_init, v20elec_init,
                                             v30elec_init, E_phot_comv_init,
                                             x_dim_ph_ener, om_one_ph_init,
                                             om_two_ph_init, om_three_ph_init) ;

      // Assigning values returned from tuple ener_acc_rej_res from funct_scatt_acc_rej func
      phot_ener_acc_rej = get<0>(ener_acc_rej_res) ;
      E_phot_comv_fin = get<1>(ener_acc_rej_res) ;
      om_one_ph_fin = get<2>(ener_acc_rej_res) ;
      om_two_ph_fin = get<3>(ener_acc_rej_res) ;
      om_three_ph_fin = get<4>(ener_acc_rej_res) ;

    }
    // If phot_ener_acc_rej<0, scattering accepted
    // Else, draw new random numbers until phot_ener_acc_rej<0 is satisfied.
    while (phot_ener_acc_rej>=0) ;
  }

  // Else, scatt_acc_rej_cond>0, scattering event does not happen.
  // Set initial photon energy and direction equal to final photon energy and direction
  else
  {
      E_phot_comv_fin = E_phot_comv_init ;
      om_one_ph_fin = om_one_ph_init ;
      om_two_ph_fin = om_two_ph_init ;
      om_three_ph_fin = om_three_ph_init ;
  }

  //Acceptance-Rejection finished. Return final photon energy and direction after scattering
  return make_tuple(E_phot_comv_fin, om_one_ph_fin, om_two_ph_fin, om_three_ph_fin) ;
}


// Function which updates electron energy and direction after photon energy and direction accepted
tuple<double, double, double,
      double, double> update_elec_func(double gamma_e_init, double beta_e_init,
                                       double v10elec_init, double v20elec_init,
                                       double v30elec_init, double E_phot_comv_init,
                                       double E_phot_comv_fin, double om_one_ph_init,
                                       double om_two_ph_init, double om_three_ph_init,
                                       double om_one_ph_fin, double om_two_ph_fin,
                                       double om_three_ph_fin)
{

   double gamma_e_fin, beta_e_fin, v10elec_fin, v20elec_fin, v30elec_fin ;

  // Updating gamma_e and beta_e of electron after scattering with energy conservation
  if( (E_phot_comv_init - E_phot_comv_fin + m_e*c*c*gamma_e_init)/(m_e*c*c) > 1.0 )
  {
    gamma_e_fin = (E_phot_comv_init - E_phot_comv_fin + m_e*c*c*gamma_e_init)/(m_e*c*c) ;
    beta_e_fin = sqrt(1.0 - 1.0/(gamma_e_fin*gamma_e_fin)) ;
  }

  else
  {
    gamma_e_fin = 1.00001 ;
    beta_e_fin = sqrt(1.0 - 1.0/(gamma_e_fin*gamma_e_fin)) ;
  }

  // Updating v10elec, v20elec, and v30elec after scattering with momentum conservation
  v10elec_fin = (  ( (E_phot_comv_init/c)*om_one_ph_init +
                               m_e*c*beta_e_init*gamma_e_init*v10elec_init -
                               (E_phot_comv_fin/c)*om_one_ph_fin )/
                               (m_e*c*beta_e_fin*gamma_e_fin)  ) ;


  v20elec_fin = (  ( (E_phot_comv_init/c)*om_two_ph_init +
                             m_e*c*beta_e_init*gamma_e_init*v20elec_init -
                               (E_phot_comv_fin/c)*om_two_ph_fin )/
                               (m_e*c*beta_e_fin*gamma_e_fin)  ) ;


  // Decrease v10elec_fin,v20elec_fin if round-off causes normalization condition to be imaginary
  if( v10elec_fin*v10elec_fin + v20elec_fin*v20elec_fin >= 1.0)
  {
    v10elec_fin = 0.9999*v10elec_fin ;
    v20elec_fin = 0.9999*v20elec_fin ;
  }

  // Using the normalization condition to update v30elec_fin. Since denominator of
  // v30elec_fin is positive, only using the numerator of v30elec_fin for copysign
  v30elec_fin = copysign( sqrt( 1.0 - v10elec_fin*v10elec_fin - v20elec_fin*v20elec_fin) ,
                              ( (E_phot_comv_init/c)*om_three_ph_init +
                               m_e*c*beta_e_init*gamma_e_init*v30elec_init -
                               (E_phot_comv_fin/c)*om_three_ph_fin ) )  ;

  return make_tuple(gamma_e_fin, beta_e_fin, v10elec_fin, v20elec_fin, v30elec_fin) ;

}






