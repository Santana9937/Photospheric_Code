#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <tuple>
#include <random>
#include <queue>
#include <utility>
#include <vector>
#include <algorithm>
#include <string>

// Using namespace std to not call standard library with std:: , standard C++ practice
using namespace std ;

// Declaring and assigning fundamental constants
#include "fund_const.h"
// Declaring and assigning input parameters for simulation
#include "input_parameters.h"
// Variable declarations for electrons and void function which initializes
// the direction and gamma_e of the electrons. There is another function in
// this header file, which updates the gamma_e of the electrons we are tracking
#include "electron_properties.h"
// Variable declarations for photons and void function which initializes
// the direction and energy of the photons. There is another function which
// calculates the new location of the photon and the elapsed time of the photon
// (both in the observer frame) after a s_disp drawn in the jet-comoving frame.
#include "photon_properties.h"
// This header file contains the function which draws new directions and gamma_e
// for the electrons after a re-heating episode. The variable declarations for
// electron re-heating are also in this header file
#include "electron_reheat_properties.h"
// Contains a function which updates photon energy and electron
// gam_e, beta_e due to adiabatic cooling
#include "adiab_cool_ph_and_el.h"
// This header file contains 3 function related to the IC/Compton scattering.
// The first function calculates the photon energy and direction after scattering
// and calculates an acceptance rejection quantity. The second function checks to
// see if the energy and direction of the photon after scattering are accepted. The
// last function updates the energy and direction of the electrons after scattering
#include "IC_Compton_scatt.h"
// Function which converts simulation parameters to strings for file-saving.
// A function which bins the photon and electron energies and finds the f_nu
// spectrum is also included. Thus function also writes the spectrum to file.
// Lastly, this header file also contains a function which creates a log-file
// which saves the simulation properties.
#include "filename_bin_save_logfile.h"
// This header-file includes the function prototypes of all the functions in the program
// Tuples to extract values from function are also declared in this file
#include "func_proto_tuple_decl.h"
// This header file includes a function which prints the simulation parameters
// and properties to the command prompt
#include "print_sim_prop.h"


int main(void)
{
  // This function prints simulation parameters and properties to the screen
  // It also returns an integer to determine if proper input parameters were specified
  bool proper_sim_parm = print_sim_param_and_props() ;

  // proper_sim_parm is 0 (i.e. false), simulation is terminated
  if(proper_sim_parm==false)
  {
    cout << "PROGRAM HAS BEEN TERMINATED" << endl << endl ;  ;
    return 0 ;
  }


  // Start timing code
  clock_t start_time = clock() ;

  // Initializing photons and electrons
  initialize_electrons() ;
  initialize_photons() ;


  /* ------ Declaring Random Number Generator Information For Drawing Random Electrons -------- */

  // declare a random number generator
  // here we choose the Mersenne Twister engine (32bit)
  mt19937 rng;

  random_device randomSeed; // this uses hardware parameters (clocks etc) to generate a random seed
  rng.seed(randomSeed()); // set the seed from the random device

  /* ---------------------------------------------------------------- */

  // ****************************************************************************************
  // ****************************************************************************************
  // --- Checking if a Photon Escaped Photosphere With First Propagation. If it did, -----
  // --- collect photon energy for spectrum. Otherwise, push photon into priority queue -----
  // ****************************************************************************************
  // ****************************************************************************************

  // Declaring t_obs_each_ph_pq priority queue. Each entry is a pair (t_tot_each_ph_obs_fr, photon_index)
  // The first entry is the total elapsed time for this particular photon and the second entry is
  // the photon_index, which serves as the tag (or label) for each photon
  priority_queue<pair<double, int> > t_obs_each_ph_pq ;

  // Temporary r_ph_pos_x, r_ph_pos_y, r_ph_pos_z to determine if photon escaped photosphere
  // Also, temporary t_delta_ph_obs_temp to fill the priority queue
  // Also, cos_theta for latitude at which photon escapes for Doopler Boosting
  double r_ph_pos_x_temp, r_ph_pos_y_temp, r_ph_pos_z_temp ;
  double t_delta_ph_obs_temp ;
  double cos_theta ;

  for(int i=0 ; i<Nphotons ; i++)
  {
    // Performing first propagation of photon
    phot_propag_step = phot_propag_func(beta_jet, bulk_gamma,
                                        r_ph_pos_x[i], r_ph_pos_y[i],
                                        r_ph_pos_z[i], s_disp_pr[i],
                                        om_one_ph[i], om_two_ph[i],
                                        om_three_ph[i]) ;

    // Extracting new x, y, and x position of photon in observer frame. Also, extracting
    // time elapsed in between photon scatterings in the observer frame.
    r_ph_pos_x_temp = get<0>(phot_propag_step) ;
    r_ph_pos_y_temp = get<1>(phot_propag_step) ;
    r_ph_pos_z_temp = get<2>(phot_propag_step) ;
    t_delta_ph_obs_temp = get<3>(phot_propag_step) ;

    // If radial position of photon > R_photosphere, collect energy of photon and DO NOT PUSH INTO QUEUE
    // If N_photon_collect photons escaped photosphere (or more), stop the loop
    if( sqrt( r_ph_pos_x_temp*r_ph_pos_x_temp + r_ph_pos_y_temp*r_ph_pos_y_temp +
              r_ph_pos_z_temp*r_ph_pos_z_temp ) > r_photosphere &&  N_ph_esc_counter<N_photon_collect )
    {
      // Add one to photon escape counter
      N_ph_esc_counter ++ ;

      // Cosine of latitude where photon escapes for Doppler boosting
      cos_theta = r_ph_pos_z[i]/sqrt( r_ph_pos_x[i]*r_ph_pos_x[i] +
                                      r_ph_pos_y[i]*r_ph_pos_y[i] +
                                      r_ph_pos_z[i]*r_ph_pos_z[i] ) ;

      // Add photon to array storring final photon energies and final number of scatterings, respectively
      fin_ph_ener_obs_eV[N_ph_esc_counter-1] = (E_phot_comv[i]/(1.6022e-12))*(1.0/bulk_gamma)*
                                                 (1.0/(1.0 - beta_jet*cos_theta))  ;

      fin_num_scatt_ph_esc[N_ph_esc_counter-1] = N_scatt_each_ph[i] ;

      // DO NOT PUSH PHOTON INTO PRIORITY QUEUE!!!
    }

    // Else, photon did not escape, put it in the priority queue
    // Filling the priority with negative values since the top value is the max value
    else
    {
      t_obs_each_ph_pq.push(make_pair(-t_delta_ph_obs_temp, i));
    }

  }

  cout << "Check for photon escape in 1st propagation done and priority queue built" << endl ;

  /* ------- Declaring variables before while loop to avoid re-declaring ------- */

  // Storing elapsed time step for current s_prime of photon
  double t_delta_ph_temp ;

  // Index for photon (first in priority queue) and electron (randon electron) for scattering
  int curr_ph_ind, curr_elec_ind ;

  // Drawing random integer [0,Nelectrons-1] to choose random electron
  uniform_int_distribution<int> rand_elec(0,Nelectrons-1) ;

  // Saving properties of photon prior to scattering
  double E_phot_comv_bef_scatt, om_one_ph_bef_scatt, om_two_ph_bef_scatt, om_three_ph_bef_scatt ;

  // Drawing uniform random numbers in (0.0,1.0) for new s_disp
  uniform_real_distribution<double> rn_l_mfp_pr(0.0,1.0) ;

  // After photon propagation, use current radial position of photon location to calculate new s_disp_pr
  double r_ph_disp_aft_scatt ;

  cout << "Simulation has now begun" << endl ;


  // *************************************************************************************
  /* --------------------------------- Main Program ----------------------------------- */
  // *************************************************************************************

  // If N_photon_collect photons (or more) escaped R_photosphere, stop the simulation
  if(N_ph_esc_counter==N_photon_collect)
  {
    cout << endl ;
    cout << ">= N_photon_collect PHOTONS ESCAPED DURING FIRST PROPAGATION, SIMULATION IS DONE!" << endl ;
    cout << endl ;
  }

  // Else, start the simulation and continue until N_photon_collect have escaped R_photosphere
  else
  {
    // This while loop scatters photon with smallest s_disp_pr with a random electron.
    // Then, a new s_disp_pr is drawn and photon is propagated forward. Then, it is
    // checked if photon has escaped photosphere. Simulation continues until N_photon_collect of photons collected.
    while(N_ph_esc_counter<N_photon_collect)
    {
      // Checking electron re-heating knob to see if it was activated
      switch(reheating_knob)
      {
        case 0: // Case 0 is electron re-heating OFF
        {
          break ;
        }

        case 1: // Case 1 is electron re-heating ON
        {
          // Checking for initializing a electron re-heating episode
          if(N_scatt_counter==N_scatt_for_diss_event && N_curr_reheat_event<=N_tot_reheat_events)
          {
            cout << endl << "Re-reating event number: " << N_curr_reheat_event << "/" << N_tot_reheat_events << endl ;
            cout << "gamma_e[0] before reheating: " << gamma_e[0] << endl ;

            // This function re-heats electrons and draws new random directions for electron
            // If electron tracking is on, gamma_e also updated in vector tracking gamma_e values
            reheat_electrons() ;

            cout << "gamma_e[0] after reheating: " << gamma_e[0] << endl ;

            // Moving on to the next re-heating event
            N_curr_reheat_event ++ ;

            // Re-defining the number of photons collected for the next re-heating event
            N_scatt_for_diss_event = N_curr_reheat_event*((2.0*tau_initial*Nphotons)/(N_tot_reheat_events+1)) ;
          }

          break ;
        }

      }

      // Extract index of photon with smallest s_disp_pr from priority queue
      curr_ph_ind = t_obs_each_ph_pq.top().second ;

      // Now that curr_ph_ind has been saved, pop entry with smallest s_disp_pr from priority queue
      t_obs_each_ph_pq.pop() ;

      // Add one more scatter to photon with index curr_ph_ind
      N_scatt_each_ph[curr_ph_ind] ++ ;

      // Add one more scatter to the scattering counter
      N_scatt_counter ++ ;

      // Draw random electron to scatter with photon curr_ph_ind
      curr_elec_ind = rand_elec(rng) ;

      // Now that photon has been selected, propagate this photon forward
      phot_propag_step = phot_propag_func(beta_jet, bulk_gamma,
                                          r_ph_pos_x[curr_ph_ind], r_ph_pos_y[curr_ph_ind],
                                          r_ph_pos_z[curr_ph_ind], s_disp_pr[curr_ph_ind],
                                          om_one_ph[curr_ph_ind], om_two_ph[curr_ph_ind],
                                          om_three_ph[curr_ph_ind]) ;

      // Updating the x, y, and z location of the photon in the observer frame
      // Storring the propagation time with this s_prime as a temporary value
      // Note, not updating t_tot_each_ph_obs_fr[curr_ph_ind] to know the time
      // elapsed of the photon before propagation with s_prime. This allows us
      // to know the location of the photon before the propagation with s_prime
      r_ph_pos_x[curr_ph_ind] = get<0>(phot_propag_step) ;
      r_ph_pos_y[curr_ph_ind] = get<1>(phot_propag_step) ;
      r_ph_pos_z[curr_ph_ind] = get<2>(phot_propag_step) ;
      t_delta_ph_temp = get<3>(phot_propag_step) ;

      // ***********************************************************************
      // ------------------- Performing Adiabatic Cooling ----------------------
      // ***********************************************************************

      switch(adiab_cool_knob)
      {
        case 0: // Case 0 is adiabatic cooling OFF
        {
         break ;
        }
        case 1: // Case 1 is adiabatic cooling ON
        {
          // Extracting photon energy and electron gamma_e, beta_e from adiabatic cooling function
          abiab_cool_ph_el_return = func_adiab_cool_ph_el(t_tot_each_ph_obs_fr[curr_ph_ind],
                                                          t_delta_ph_temp, t_tot_each_elec_obs[curr_elec_ind],
                                                          E_phot_comv[curr_ph_ind], gamma_e[curr_elec_ind]) ;

          // Updating photon energy and electron gamma_e, beta_e due to adiabatic cooling
          E_phot_comv[curr_ph_ind] = get<0>(abiab_cool_ph_el_return) ;
          gamma_e[curr_elec_ind] =  get<1>(abiab_cool_ph_el_return) ;
          beta_e[curr_elec_ind] = get<2>(abiab_cool_ph_el_return) ;

          // If only adiabatic cooling is on, track electron after adiabatic cooling done
          // If adiabatic cooling + IC/Comp scattering on, track electron after both adiabatic
          // cooling and IC/Comp scattering performed
          if(scatt_cool_knob == 0)
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
                // If curr_elec_ind == index of one of the electrons we are tracking, update
                // its gamma_e value after both adiabatic cooling and IC/Comp scattering
                elec_track_switch_push(curr_elec_ind) ;

                break ;
              }
            }
          }

          break ;
        }
      }

      // ***********************************************************************
      // ---------------------- Adiabatic Cooling Done -------------------------
      // ***********************************************************************

      // Update the propagation time for the photon and electron
      // Even if adiabatic cooling is off, still important to update propagation
      // time for photons since photon priority_queue based on photon propagation time
      t_tot_each_ph_obs_fr[curr_ph_ind] = t_tot_each_ph_obs_fr[curr_ph_ind] + t_delta_ph_temp ;
      t_tot_each_elec_obs[curr_elec_ind] = t_tot_each_ph_obs_fr[curr_ph_ind] ;

      // ***********************************************************************
      // ------------ Performing Photon-Electron Scattering --------------------
      // ***********************************************************************

      switch(scatt_cool_knob)
      {
        case 0: // Case 0 is IC/Comp Scattering OFF
        {
          break ;
        }
        case 1: // Case 1 is IC/Comp Scattering ON
        {
          // Saving properties of photon prior to scattering
          E_phot_comv_bef_scatt = E_phot_comv[curr_ph_ind] ;
          om_one_ph_bef_scatt = om_one_ph[curr_ph_ind] ;
          om_two_ph_bef_scatt = om_two_ph[curr_ph_ind] ;
          om_three_ph_bef_scatt = om_three_ph[curr_ph_ind] ;

          // Scatter photon with index curr_ph_ind with electron with index curr_elec_ind
          phot_IC_scatt_return = funct_IC_Comp_scattering(beta_e[curr_elec_ind],
                                                          gamma_e[curr_elec_ind], E_phot_comv_bef_scatt,
                                                          v10elec[curr_elec_ind], v20elec[curr_elec_ind],
                                                          v30elec[curr_elec_ind], om_one_ph_bef_scatt,
                                                          om_two_ph_bef_scatt, om_three_ph_bef_scatt) ;

          // Updating photon energy and direction after scattering
          E_phot_comv[curr_ph_ind] = get<0>(phot_IC_scatt_return) ;
          om_one_ph[curr_ph_ind] = get<1>(phot_IC_scatt_return) ;
          om_two_ph[curr_ph_ind] = get<2>(phot_IC_scatt_return) ;
          om_three_ph[curr_ph_ind] = get<3>(phot_IC_scatt_return) ;

          // Now that scattering is done, update electron energy and direction
          elec_IC_scatt_return = update_elec_func(gamma_e[curr_elec_ind], beta_e[curr_elec_ind],
                                                  v10elec[curr_elec_ind], v20elec[curr_elec_ind],
                                                  v30elec[curr_elec_ind], E_phot_comv_bef_scatt,
                                                  E_phot_comv[curr_ph_ind], om_one_ph_bef_scatt,
                                                  om_two_ph_bef_scatt, om_three_ph_bef_scatt,
                                                  om_one_ph[curr_ph_ind], om_two_ph[curr_ph_ind],
                                                  om_three_ph[curr_ph_ind]) ;

          // Updating electron gamma_e, beta_e, and direction after scattering
          gamma_e[curr_elec_ind] = get<0>(elec_IC_scatt_return) ;
          beta_e[curr_elec_ind] = get<1>(elec_IC_scatt_return) ;
          v10elec[curr_elec_ind] = get<2>(elec_IC_scatt_return) ;
          v20elec[curr_elec_ind] = get<3>(elec_IC_scatt_return) ;
          v30elec[curr_elec_ind] = get<4>(elec_IC_scatt_return) ;

          // If electron-tracking is on, update the gamma_e of the 3 electrons we are tracking
          switch(elec_tracking_knob)
          {
            case 0: // Case 0 is electron tracking OFF
            {
              break ;
            }
            case 1: // Case 1 is electron tracking ON
            {
              // If curr_elec_ind == index of one of the electrons we are tracking, update
              // its gamma_e value after both adiabatic cooling and IC/Comp scattering
              elec_track_switch_push(curr_elec_ind) ;

              break ;
            }
          }

          break ;
        }

      }

      // ***********************************************************************
      // --------------- Photon-Electron Scattering Done -----------------------
      // ***********************************************************************

      // ***********************************************************************
      // -----------------------------------------------------------------------
      // - Checking to see if Photon is able to escape in the next propagation -
      // -----------------------------------------------------------------------
      // ***********************************************************************

      // Now that scattering is done, use current location of photon in obsever
      // frame to calculate the new s_disp_pr and to propagate photon forward
      r_ph_disp_aft_scatt = sqrt( r_ph_pos_x[curr_ph_ind]*r_ph_pos_x[curr_ph_ind] +
                                  r_ph_pos_y[curr_ph_ind]*r_ph_pos_y[curr_ph_ind] +
                                  r_ph_pos_z[curr_ph_ind]*r_ph_pos_z[curr_ph_ind] ) ;

      // Drawing new random number for s_disp_pr
      rnum_l_mfp_pr = rn_l_mfp_pr(rng) ;

      // Updating s_disp_pr of photon curr_ph_ind with new location in observer frame
      s_disp_pr[curr_ph_ind] = -((4.0*PI*r_ph_disp_aft_scatt*r_ph_disp_aft_scatt*
                                  m_p*c*c*c*bulk_gamma*bulk_gamma)/(L*sigma_T))*log(rnum_l_mfp_pr) ;

      // Now, propagating photon with index curr_ph_ind in observer frame with new s_disp_pr
      phot_propag_step = phot_propag_func(beta_jet, bulk_gamma,
                                          r_ph_pos_x[curr_ph_ind], r_ph_pos_y[curr_ph_ind],
                                          r_ph_pos_z[curr_ph_ind],
                                          s_disp_pr[curr_ph_ind], om_one_ph[curr_ph_ind],
                                          om_two_ph[curr_ph_ind], om_three_ph[curr_ph_ind]) ;

      // Extracting new x, y, and x position of photon in observer frame. Also, extracting
      // time elapsed in between photon scatterings in the observer frame.
      // Temporary r_ph_pos_x, r_ph_pos_y, r_ph_pos_z to determine if photon escaped photosphere
      // Also, temporary t_delta_ph_obs_temp to fill the priority queue
      r_ph_pos_x_temp = get<0>(phot_propag_step) ;
      r_ph_pos_y_temp = get<1>(phot_propag_step) ;
      r_ph_pos_z_temp = get<2>(phot_propag_step) ;
      t_delta_ph_obs_temp = get<3>(phot_propag_step) ;

      // Check to see if photon with curr_ph_ind escaped with new s_disp_pr
      if( sqrt( r_ph_pos_x_temp*r_ph_pos_x_temp +
                r_ph_pos_y_temp*r_ph_pos_y_temp +
                r_ph_pos_z_temp*r_ph_pos_z_temp ) > r_photosphere )
      {
        // If so, update escape counter and store energy and scatterings for this photon
        N_ph_esc_counter ++ ;

        // Cosine of latitude where photon escapes for Doppler boosting
        cos_theta = r_ph_pos_z[curr_ph_ind]/sqrt( r_ph_pos_x[curr_ph_ind]*r_ph_pos_x[curr_ph_ind] +
                                                  r_ph_pos_y[curr_ph_ind]*r_ph_pos_y[curr_ph_ind] +
                                                  r_ph_pos_z[curr_ph_ind]*r_ph_pos_z[curr_ph_ind] ) ;

        fin_ph_ener_obs_eV[N_ph_esc_counter-1] = (E_phot_comv[curr_ph_ind]/(1.6022e-12))*(1.0/bulk_gamma)*
                                                 (1.0/(1.0 - beta_jet*cos_theta))  ;

        fin_num_scatt_ph_esc[N_ph_esc_counter-1] = N_scatt_each_ph[curr_ph_ind] ;

        // DO NOT PUSH PHOTON INTO PRIORITY QUEUE!!!
      }

      // Else, photon did not escape, put it in the priority queue
      // Filling the priority with negative values since the top value is the max value
      else
      {
        // Pushing the sum of the total previous elapsed time for this photon plus
        // the new elapsed time with the new s_disp
        t_obs_each_ph_pq.push(make_pair(-t_tot_each_ph_obs_fr[curr_ph_ind]-t_delta_ph_obs_temp, curr_ph_ind));
      }

    }

  }



  // Stop timing code
  clock_t end_time = clock() ;
  // Find the time the code took to run in sec
  double timeSec = (end_time - start_time) / static_cast<double>( CLOCKS_PER_SEC );
  // Output the time the code took to run in sec
  cout << endl ;
  cout << "Code took: " << timeSec << " sec" << endl ;


  // ****************************************************************************
  // ------ Formating Variable Names to Strings for File Saving and -------------
  // --------------- Storring Simulation Results to Binary File -----------------
  // ****************************************************************************

  // This function saves the photons and electrons and returns a string
  // for the log file filename
  string sim_properties = str_sim_properties_func() ;

  fnu_bin_and_file_save(sim_properties) ;



  // ------- Average # Scatterings of all photons and photons that escaped ----

  // Declaring variables used to hold sums
  double N_scatt_each_ph_sum = 0, fin_num_scatt_ph_esc_sum = 0 ;

  // Calculating the sum of number of scatterings for all the photons
  for(int i=0 ; i<Nphotons ; i++)
  {
    N_scatt_each_ph_sum += N_scatt_each_ph[i] ;
  }

  // Calculating the sum of number of scatterings for photons that escaped photosphere
  for(int i=0 ; i<N_photon_collect ; i++)
  {
    fin_num_scatt_ph_esc_sum += fin_num_scatt_ph_esc[i] ;
  }

  cout << endl ;
  cout << "Ave. Num Scatt for Escape Phot: " << fin_num_scatt_ph_esc_sum/N_photon_collect << endl ;
  cout << "Ave. Scatt All Photons: " << N_scatt_each_ph_sum/Nphotons << endl ;
  cout << endl ;

  //for(int i=0 ; i<Nelectrons ; i++)
  //{
  //  cout << "gamma_e[" << i << "]: " << gamma_e[i] << endl ;
  //}
  //cout << endl ;

  //for(int i=0 ; i<N_photon_collect ; i++)
  //{
    //cout << "fin_num_scatt_ph_esc[" << i << "]: " << fin_ph_ener_obs_eV[i] << endl ;
  //}
  //cout << endl ;

  // ****************************************************************************
  // ---------------- Saving Logfile With Simulation Properties -----------------
  // ****************************************************************************
  logfile_func(sim_properties, fin_num_scatt_ph_esc_sum, N_scatt_each_ph_sum, timeSec) ;

  return 0 ;

}
