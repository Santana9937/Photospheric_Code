// --------------------- Photon Variable Declarations -------------------------

// Specifying Energy of Photons in Jet Comoving Frame ------------

// Array for energy of photons in jet comoving frame
vector<double> E_phot_comv(Nphotons) ;

// Array needed to draw BB photons
double jsum[50] ;

// Declaring 4 random numbers needed for each BB photon
double rnum_BB1, rnum_BB2, rnum_BB3, rnum_BB4 ;

// Declaring 1 random number need for each PL distribution photon
double rnum_phot_PL_dist ;


// Specifying Momentum of Photons in Jet Comoving Frame -----------

// 3 arrays for components of photon in x, y, and z directions
vector<double> om_one_ph(Nphotons), om_two_ph(Nphotons), om_three_ph(Nphotons) ;

// Declaring 2 random for the direction of each BB photon
double rnum_phot_dir_one, rnum_phot_dir_two ;



// Specifying position of photons in observer frame at start of simulation ---------

// Initial radius and photospheric radius
const double r_initial = (L*sigma_T)/(8.0*PI*tau_initial*m_p*c*c*c*bulk_gamma*bulk_gamma*bulk_gamma) ;
const double r_photosphere = (L*sigma_T)/(8.0*PI*tau_photosphere*m_p*c*c*c*bulk_gamma*bulk_gamma*bulk_gamma) ;

// 3 arrays for location of each photon in x, y, and z direction
vector<double> r_ph_pos_x(Nphotons), r_ph_pos_y(Nphotons), r_ph_pos_z(Nphotons) ;

// Array to store the propagation time for each photon in the observer frame
vector<double> t_tot_each_ph_obs_fr(Nphotons) ;

// Placing each photon in a random theta and phi in the jet
double theta_ph_pos, phi_ph_pos ;

// 2 random numbers for photon location
double rnum_phot_pos_theta, rnum_phot_pos_phi ;

// array for distance traveled by photon in each step by sampling mean-free path
vector<double> s_disp_pr(Nphotons) ;

// random number for drawin s_disp_pr array
double rnum_l_mfp_pr ;



// Quantities keeping track of for photons ---------------------

// Counter for the number of photons that have currently escaped the photosphere
int N_ph_esc_counter = 0 ;

// Array holding the number of scatterings each photon experiences
vector<int> N_scatt_each_ph(Nphotons) ;

// Final scatterings for collected photons
vector<int> fin_num_scatt_ph_esc(N_photon_collect) ;
// Final photon energies in eV Lorentz transformed to the observer frame
vector<double> fin_ph_ener_obs_eV(N_photon_collect) ;



// This void function initializes photon momentum, energies, and photon position in observer frame
void initialize_photons()
{

  /* --------- Declaring Random Number Generator Information -------- */

  // declare a random number generator
  // here we choose the Mersenne Twister engine (32bit)
  mt19937 rng;

  random_device randomSeed; // this uses hardware parameters (clocks etc) to generate a random seed
  rng.seed(randomSeed()); // set the seed from the random device


  /* --------------------- Drawing Photon Momentum Direction --------------------- */

  // Drawing 2 random numbers for photon direction
  uniform_real_distribution<double> rn1_phot_dir(0.0,1.0) ;
  uniform_real_distribution<double> rn2_phot_dir(0.0,1.0) ;

  // Filling in direction of photons with for loop
  for(int i=0 ; i<Nphotons ; i++)
  {
    // Getting 2 random numbers
    rnum_phot_dir_one = rn1_phot_dir(rng), rnum_phot_dir_two = rn2_phot_dir(rng) ;

    // Filling in direction of photons
    om_three_ph[i] = 2.0*rnum_phot_dir_one - 1.0 ;
    om_two_ph[i] = ( sqrt(1.0 - om_three_ph[i]*om_three_ph[i]) )*sin(2.0*PI*rnum_phot_dir_two) ;
    om_one_ph[i] = ( sqrt(1.0 - om_three_ph[i]*om_three_ph[i]) )*cos(2.0*PI*rnum_phot_dir_two) ;
  }

  /* --------------------------- Drawing BB energies -------------------------- */

  // Filling in the first entry of jsum
  jsum[0] = 1.0 ;

  // Filling in the subsequent entries for jsum
  for(int i=1 ; i<50 ; i++)
  {
    jsum[i] = jsum[i-1] + 1.0/((i+1)*(i+1)*(i+1)) ;
  }

  // Inputing jsum entries into a vector
  vector<double> jsum_vect(jsum,jsum+50) ;
  // Defining iterator for vector
  vector<double>::iterator low ;

  // Drawing 4 random numbers for each BB photon
  uniform_real_distribution<double> BB_rn1(0.0,1.0) ;
  uniform_real_distribution<double> BB_rn2(0.0,1.0) ;
  uniform_real_distribution<double> BB_rn3(0.0,1.0) ;
  uniform_real_distribution<double> BB_rn4(0.0,1.0) ;

  // Declaring integer alpha needed for BB energies
  int alpha ;

  // Drawing 1 random number for power-law distribution of photons
  uniform_real_distribution<double> rn_phot_PL_dist(0.0,1.0) ;

  // Filling in energy of photons to follow BB or PL distribution
  for(int i=0 ; i<Nphotons ; i++)
  {
    switch(photon_distr_type)
    {
      case 1: // Case 1 is BB photons
      {
        // Getting 4 random numbers
        rnum_BB1 = BB_rn1(rng), rnum_BB2 = BB_rn2(rng), rnum_BB3 = BB_rn3(rng), rnum_BB4 = BB_rn4(rng) ;

        // Using lower bound to calculate alpha
        low =lower_bound (jsum_vect.begin(), jsum_vect.end(), 1.202*rnum_BB1);
        alpha = low - jsum_vect.begin() + 1 ;

        // Assining the energy of each BB photon
        E_phot_comv[i] = -((kB*T_phot_comv)/alpha)*log(rnum_BB2*rnum_BB3*rnum_BB4) ;

        break ;
      }
      case 2: // Case 2 is PL photons
      {
        // Getting 1 random numbers for photon energy
        rnum_phot_PL_dist = rn_phot_PL_dist(rng) ;

        // Assining the energy of each BB photon
        E_phot_comv[i] = ( pow( ( rnum_phot_PL_dist*( pow(E_2_phot_PL, (1.0-p_phot)) - pow(E_1_phot_PL, (1.0-p_phot)) ) + pow(E_1_phot_PL, (1.0-p_phot)) ), (1.0/(1.0-p_phot)) ) ) ;

        break ;
      }
    }
  }


  switch(photon_distr_type)
  {
    case 1: // Case 1 is BB photons
    {
      cout << "BB photons drawn" << endl ;
      break ;
    }
    case 2: // Case 2 is PL photons
    {
      cout << "PL photons drawn" << endl ;
      break ;
    }
  }

  // ---------- Initializing photons in jet and performing first propagation -----------

  // Placing observer in +z-direction. Photons escape when r>r_photosphere in hemisphere in +z

  // Photons initialized in r=r_initial in a cone with half-opening angle 1/bulk_gamma with respect
  // to the +z-axis and radius r_initial. (phi, theta) location for each photon chosen randomly

  // Since all photons initialized at r_initial, initial lmfp is the same for all photons
  const double lmfp_initial = (4.0*PI*r_initial*r_initial*m_p*c*c*c*bulk_gamma*bulk_gamma)/(L*sigma_T) ;


  // Drawing 2 random number for photon location in the jet
  uniform_real_distribution<double> rn_phot_pos_theta(0.0,1.0) ;
  uniform_real_distribution<double> rn_phot_pos_phi(0.0,1.0) ;

  // Drawing random number for s_disp
  uniform_real_distribution<double> rn_l_mfp_pr(0.0,1.0) ;

  // This for loop:
  // 1. initializes location of each photon in jet
  // 2. calculates initial energy and momentum in each direction for all photons
  for(int i=0 ; i<Nphotons ; i++)
  {
    // Assigning 2 random numbers for photon location
    rnum_phot_pos_theta = rn_phot_pos_theta(rng), rnum_phot_pos_phi = rn_phot_pos_phi(rng) ;

    // Theta and Phi from random number
    theta_ph_pos = (1.0/bulk_gamma)*rnum_phot_pos_theta, phi_ph_pos = 2.0*PI*rnum_phot_pos_phi ;

    // Initial Position of the Photons
    r_ph_pos_x[i] = r_initial*sin(theta_ph_pos)*cos(phi_ph_pos) ;
    r_ph_pos_y[i] = r_initial*sin(theta_ph_pos)*sin(phi_ph_pos) ;
    r_ph_pos_z[i] = r_initial*cos(theta_ph_pos) ;

    // Drawing random number to calculate s_disp_pr
    rnum_l_mfp_pr = rn_l_mfp_pr(rng) ;
    s_disp_pr[i] = -lmfp_initial*log(rnum_l_mfp_pr) ;
  }

  // Filling in all entries in N_scatt_each_ph with a value 1
  fill(N_scatt_each_ph.begin(), N_scatt_each_ph.end(), 0) ;

  // Filling in all entries of t_tot_each_ph_obs_fr to 0.0
  fill(t_tot_each_ph_obs_fr.begin(), t_tot_each_ph_obs_fr.end(), 0.0) ;
}




tuple<double, double,
      double, double> phot_propag_func(double beta_jet, double bulk_gamma,
                                       double r_ph_pos_x_init, double r_ph_pos_y_init,
                                       double r_ph_pos_z_init,
                                       double s_disp_pr_fin, double om_one_ph_fin,
                                       double om_two_ph_fin, double om_three_ph_fin)
{
  // Declaring the variables that are going to be returned by the function
  double r_ph_pos_x_fin, r_ph_pos_y_fin, r_ph_pos_z_fin, t_delta_fin ;

  // Initial spherical theta and phi angles
  double theta_ph_pos_init = acos(r_ph_pos_z_init/sqrt(r_ph_pos_x_init*r_ph_pos_x_init +
                                                       r_ph_pos_y_init*r_ph_pos_y_init +
                                                       r_ph_pos_z_init*r_ph_pos_z_init ) ) ;

  double phi_ph_pos_init = atan2(r_ph_pos_y_init , r_ph_pos_x_init) ;

  //----- Propogating Photons now that s_disp_prime has been drawn for each photon -----

  // Ratio of beta_jet_x/beta_jet, beta_jet_y/beta_jet, and beta_jet_z/beta_jet
  double betajx_div_betaj = sin(theta_ph_pos_init)*cos(phi_ph_pos_init) ;
  double betajy_div_betaj = sin(theta_ph_pos_init)*sin(phi_ph_pos_init) ;
  double betajz_div_betaj = cos(theta_ph_pos_init) ;

  // Propagation time of photon in the jet comoving frame (t_pr_step)
  // and progation of photon in the x (x_pr_step), y (y_pr_step), and
  // z (z_pr_step) direction in the jet comoving frame
  double t_pr_step = s_disp_pr_fin/c ;
  double x_pr_step = s_disp_pr_fin*om_one_ph_fin ;
  double y_pr_step = s_disp_pr_fin*om_two_ph_fin ;
  double z_pr_step = s_disp_pr_fin*om_three_ph_fin ;

  // Lorentz transformation of photon propagation in x-direction from jet-comoving frame to observer frame
  r_ph_pos_x_fin = (  r_ph_pos_x_init +
                      bulk_gamma*beta_jet*betajx_div_betaj*c*t_pr_step +
                      (1.0 + (bulk_gamma-1.0)*betajx_div_betaj*betajx_div_betaj)*x_pr_step +
                      (bulk_gamma-1.0)*betajx_div_betaj*betajy_div_betaj*y_pr_step +
                      (bulk_gamma-1.0)*betajx_div_betaj*betajz_div_betaj*z_pr_step ) ;

  // Lorentz transformation of photon propagation in y-direction from jet-comoving frame to observer frame
  r_ph_pos_y_fin = (  r_ph_pos_y_init +
                      bulk_gamma*beta_jet*betajy_div_betaj*c*t_pr_step +
                      (bulk_gamma-1.0)*betajy_div_betaj*betajx_div_betaj*x_pr_step +
                      (1.0 + (bulk_gamma-1.0)*betajy_div_betaj*betajy_div_betaj)*y_pr_step +
                      (bulk_gamma-1.0)*betajy_div_betaj*betajz_div_betaj*z_pr_step ) ;

  // Lorentz transformation of photon propagation in z-direction from jet-comoving frame to observer frame
  r_ph_pos_z_fin = (  r_ph_pos_z_init +
                      bulk_gamma*beta_jet*betajz_div_betaj*c*t_pr_step +
                      (bulk_gamma-1.0)*betajz_div_betaj*betajx_div_betaj*x_pr_step +
                      (bulk_gamma-1.0)*betajz_div_betaj*betajy_div_betaj*y_pr_step +
                      (1.0 + (bulk_gamma-1)*betajz_div_betaj*betajz_div_betaj)*z_pr_step ) ;

  // Lorentz transforming time to the observer frame
  t_delta_fin = bulk_gamma*t_pr_step +
                bulk_gamma*beta_jet*betajx_div_betaj*(x_pr_step/c) +
                bulk_gamma*beta_jet*betajy_div_betaj*(y_pr_step/c) +
                bulk_gamma*beta_jet*betajz_div_betaj*(z_pr_step/c) ;

  return make_tuple(r_ph_pos_x_fin, r_ph_pos_y_fin, r_ph_pos_z_fin, t_delta_fin) ;
}


