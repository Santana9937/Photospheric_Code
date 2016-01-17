// -------- Variable Declarations for all Elecron Distributions -------

// Arrays for gamma_e and beta_e for electrons
vector<double> gamma_e(Nelectrons), beta_e(Nelectrons) ;

// Declare array to hold elapsed time for each electron
vector<double> t_tot_each_elec_obs(Nelectrons) ;

// 3 arrays for components of electron in x, y, and z directions
vector<double> v10elec(Nelectrons), v20elec(Nelectrons), v30elec(Nelectrons) ;

// Declaring 2 random for the direction of each electron
double rnum_elec_dir_one, rnum_elec_dir_two ;

// ------ Additional Variable Declarations for MB Distribution Electrons -------

// Defining the electron adabatic index and electron temperature. Temp is in Kelvin
// m_e*c*c = 5.11e5 eV and kB = 8.6173e-5 eV/K
const double gam_elec_ad_ind_init = (4.0*gamma_e_init_MB + 1.0)/(3.0*gamma_e_init_MB) ;
const double T_e_comv = (5.11e5*(gam_elec_ad_ind_init-1.0)*(gamma_e_init_MB-1.0))/8.6173e-5 ;
// Dimensionless electron energy in jet comoving frame
double const n_dim_elec_ener = (kB*T_e_comv)/(m_e*c*c) ;
// Dimensionless momentum in jet comoving frame
double eta ;

// ----- Declaring Vectors for the gamma_e values of the electrons we will track

// Tracking gamma_e of some electrons after scattering
vector<double> vec_gam_e_track_1 ;
vector<double> vec_gam_e_track_2 ;
vector<double> vec_gam_e_track_3 ;


// This void function fills in the electron gamma_e, beta_e, and electron direction
// for each electron from the distribution specified
void initialize_electrons()
{

  // --------- Declaring Random Number Generator Information --------

  // declare a random number generator
  // here we choose the Mersenne Twister engine (32bit)
  mt19937 rng;

  random_device randomSeed; // this uses hardware parameters (clocks etc) to generate a random seed
  rng.seed(randomSeed()); // set the seed from the random device

  // ----------------------------------------------------------------


  // **********************************************************************************
  // ------ Drawing Electron gamma_e, beta_e from Electron distribution specified -----
  // **********************************************************************************

  switch(electron_distr_type)
  {
    case 1: // Case 1 is mono electrons
    {
      // This for loop fills in gamma_e and beta_e for each electron
      for(int i=0 ; i<Nelectrons; i++)
      {
        // Filling in gamma_e and beta_e of each electron
        gamma_e[i] = gamma_e_init_mono ;
        beta_e[i] = sqrt( 1.0 - 1.0/(gamma_e[i]*gamma_e[i]) ) ;
      }

      break ;

    }

    case 2: // Case 2 is MB electrons
    {

      // Drawing electron energies of T_e_comv < 150 keV
      if(T_e_comv < 1.5e5/8.6173e-5)
      {

        /* -------------- Drawing Low-Energy MB Random Numbers ---------------- */

        // Drawing 2 random numbers for electron direction
        uniform_real_distribution<double> low_en_mb_rn1(0.0,1.0) ;
        uniform_real_distribution<double> low_en_mb_rn2(0.0,1.0) ;

        // Declaring variables for electron energy acceptance rejection
        double rn1_low_en_mb, rn2_low_en_mb, xi_pr, low_en_mb_acc_rej ;

        // This for loop fills in gamma_e and beta_e for each electron
        for(int i=0 ; i<Nelectrons; i++)
        {
          // Drawing 2 random numbers
          rn1_low_en_mb = low_en_mb_rn1(rng), rn2_low_en_mb = low_en_mb_rn2(rng) ;

          // Auxilary random number for acceptance rejection
          xi_pr = (-3.0/2.0)*log(rn1_low_en_mb) ;

          // Dimensionless momentum if acceptance-rejection satisfied
          eta = sqrt( n_dim_elec_ener*xi_pr*(2.0 + n_dim_elec_ener*xi_pr) ) ;

          // Acceptance Rejection for low-energy MB condition
          low_en_mb_acc_rej = rn2_low_en_mb*rn2_low_en_mb -
                              0.151*(1.0 + n_dim_elec_ener*xi_pr)*(1.0 + n_dim_elec_ener*xi_pr)*
                              xi_pr*(2.0 + n_dim_elec_ener*xi_pr)*rn1_low_en_mb ;

          // If acceptance-rejection condition is not satisfied,
          // draw new randon numbers until it is satisfied
          while(low_en_mb_acc_rej>=0)
          {
            // Drawing 2 random numbers
            rn1_low_en_mb = low_en_mb_rn1(rng), rn2_low_en_mb = low_en_mb_rn2(rng) ;

            // Auxilary random number for acceptance rejection
            xi_pr = (-3.0/2.0)*log(rn1_low_en_mb) ;

            // Dimensionless momentum if acceptance-rejection satisfied
            eta = sqrt( n_dim_elec_ener*xi_pr*(2.0 + n_dim_elec_ener*xi_pr) ) ;

            // Acceptance Rejection for low-energy MB condition
            low_en_mb_acc_rej = rn2_low_en_mb*rn2_low_en_mb -
                                0.151*(1.0 + n_dim_elec_ener*xi_pr)*(1.0 + n_dim_elec_ener*xi_pr)*
                                xi_pr*(2.0 + n_dim_elec_ener*xi_pr)*rn1_low_en_mb ;
          }

          // Now that eta has been determined, fill in gamma_e[i] and beta_e[i]
          gamma_e[i] = sqrt(eta*eta + 1.0) ;
          beta_e[i] = sqrt( 1.0 - 1.0/(gamma_e[i]*gamma_e[i]) ) ;

        }

      }

      // Drawing electron energies of T_e_comv >= 150 keV
      else
      {

        // Drawing 4 random numbers for electron direction
        uniform_real_distribution<double> high_en_mb_rn1(0.0,1.0) ;
        uniform_real_distribution<double> high_en_mb_rn2(0.0,1.0) ;
        uniform_real_distribution<double> high_en_mb_rn3(0.0,1.0) ;
        uniform_real_distribution<double> high_en_mb_rn4(0.0,1.0) ;


        // Declaring variables for electron energy acceptance rejection
        double rn1_high_en_mb, rn2_high_en_mb, rn3_high_en_mb, rn4_high_en_mb ;
        double eta_pr, eta_pr_pr ;
        double high_en_mb_acc_rej ;


        for(int i=0 ; i<Nelectrons; i++)
        {

          // Drawing 4 random numbers
          rn1_high_en_mb = high_en_mb_rn1(rng), rn2_high_en_mb = high_en_mb_rn2(rng) ;
          rn3_high_en_mb = high_en_mb_rn3(rng), rn4_high_en_mb = high_en_mb_rn4(rng) ;

          // Computing 2 quantities for acceptance-rejection
          eta_pr = -n_dim_elec_ener*log(rn1_high_en_mb*rn2_high_en_mb*rn3_high_en_mb) ;
          eta_pr_pr = -n_dim_elec_ener*log(rn1_high_en_mb*rn2_high_en_mb*rn3_high_en_mb*rn4_high_en_mb) ;

          // Dimensionless momentum if acceptance-rejection satisfied
          eta = eta_pr ;

          // Acceptance Rejection for high-energy MB condition
          high_en_mb_acc_rej = eta_pr_pr*eta_pr_pr - eta_pr*eta_pr ;

          while(high_en_mb_acc_rej<=1.0)
          {

            // Drawing 4 random numbers
            rn1_high_en_mb = high_en_mb_rn1(rng), rn2_high_en_mb = high_en_mb_rn2(rng) ;
            rn3_high_en_mb = high_en_mb_rn3(rng), rn4_high_en_mb = high_en_mb_rn4(rng) ;

            // Computing 2 quantities for acceptance-rejection
            eta_pr = -n_dim_elec_ener*log(rn1_high_en_mb*rn2_high_en_mb*rn3_high_en_mb) ;
            eta_pr_pr = -n_dim_elec_ener*log(rn1_high_en_mb*rn2_high_en_mb*rn3_high_en_mb*rn4_high_en_mb) ;

            // Dimensionless momentum if acceptance-rejection satisfied
            eta = eta_pr ;

            // Acceptance Rejection for high-energy MB condition
            high_en_mb_acc_rej = eta_pr_pr*eta_pr_pr - eta_pr*eta_pr ;
          }

          // Now that eta has been determined, fill in gamma_e[i] and beta_e[i]
          gamma_e[i] = sqrt(eta*eta + 1.0) ;
          beta_e[i] = sqrt( 1.0 - 1.0/(gamma_e[i]*gamma_e[i]) ) ;
        }

      }

      break ;

    }

    case 3: // Case 3 is PL electrons
    {
      // Drawing 1 random number for power-law distribution of electrons
      uniform_real_distribution<double> rn_elec_PL_dist(0.0,1.0) ;

      // Declaring variable for extrancting random number
      double rnum_elec_PL_dist ;

      // This for loop fills in gamma_e and beta_e for each PL distr. elec
      for(int i=0 ; i<Nelectrons; i++)
      {
        // Drawing 1 random number for electron power-law distribution
        rnum_elec_PL_dist = rn_elec_PL_dist(rng) ;

        // Filling in gamma_e and beta_e of each electron
        gamma_e[i] = ( pow( ( rnum_elec_PL_dist*( pow(E_2_PL, (1.0-p_elec)) - pow(E_1_PL, (1.0-p_elec)) ) + pow(E_1_PL, (1.0-p_elec)) ), (1.0/(1.0-p_elec)) ) )/(m_e*c*c) ;
        beta_e[i] = sqrt( 1.0 - 1.0/(gamma_e[i]*gamma_e[i]) ) ;
      }

      break ;
    }

  }


  // **********************************************************************************
  // ------------------ Drawing Random Direction for each electron --------------------
  // **********************************************************************************

  // Drawing 2 random numbers for electron direction
  uniform_real_distribution<double> rn1_elec_dir(0.0,1.0) ;
  uniform_real_distribution<double> rn2_elec_dir(0.0,1.0) ;

  // This for loop fills the direction of each electron from 2 random numbers
  // Initial energy and momentum of the electrons are also calculated in this loop
  for(int i=0 ; i<Nelectrons; i++)
  {
    // Assigning the random numbers for electron direction
    rnum_elec_dir_one = rn1_elec_dir(rng), rnum_elec_dir_two = rn2_elec_dir(rng) ;

    // Filling in direction of electrons
    v30elec[i] = 2.0*rnum_elec_dir_one - 1.0 ;
    v20elec[i] = sqrt(1.0 - v30elec[i]*v30elec[i])*sin(2.0*PI*rnum_elec_dir_two) ;
    v10elec[i] = sqrt(1.0 - v30elec[i]*v30elec[i])*cos(2.0*PI*rnum_elec_dir_two) ;

  }


  // Filling in all entries of the total time elapsed for each electron to 0.0
  fill(t_tot_each_elec_obs.begin(), t_tot_each_elec_obs.end(), 0.0) ;

  // If electron tracking knock is on, initialize the
  // gamma_e of the 3 electrons we are tracking

  switch(elec_tracking_knob)
  {
    case 0: // Case 0 is electron tracking OFF
    {
      break ;
    }
    case 1: // Case 1 is electron tracking ON
    {
      vec_gam_e_track_1.push_back (gamma_e[elec_track_1]) ;
      vec_gam_e_track_2.push_back (gamma_e[elec_track_2]) ;
      vec_gam_e_track_3.push_back (gamma_e[elec_track_3]) ;
      break ;
    }
  }

}


void elec_track_switch_push(int curr_elec_ind)
{
  // If electron-tracking is on, update the gamma_e of the 3 electrons we are tracking
  switch(elec_tracking_knob)
  {
    case 0: // Case 0 is electron tracking OFF
    {
      break ;
    }
    case 1: // Case 1 is electron tracking ON
    {
      if(curr_elec_ind==elec_track_1)
      {
        vec_gam_e_track_1.push_back (gamma_e[elec_track_1]) ;
      }
      else if(curr_elec_ind==elec_track_2)
      {
        vec_gam_e_track_2.push_back (gamma_e[elec_track_2]) ;
      }
      else if(curr_elec_ind==elec_track_3)
      {
        vec_gam_e_track_3.push_back (gamma_e[elec_track_3]) ;
      }
      break ;
    }
  }
}

