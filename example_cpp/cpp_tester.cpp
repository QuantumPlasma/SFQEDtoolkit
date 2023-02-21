/*
! This program shows how to call the main routines provided by SFQEDtoolkit.
! The code at hand mimics what happens inside a generic particle pusher after
! the fields acting on a particle and its new momentum have been established.
! These quantities, which will not be specified in the following, can easily
! be used to determine the particle's quantum parameter chi and gamma, that
! represent our starting point. Therefore before actually moving the particle,
! electron/photon, the user could resort to this suite of functiuons to
! ascertain whether the electron/photon emits/decays /into a synchrophoton/pair.
! Of course, in order to implement these stochastic mechanism we will need some 
! random numbers, which in this case will be fixed.
*/

#include "SFQEDtoolkit_Interface.hpp"

#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <cmath>

double scalar_prod(const double * const arr1, const double * const arr2){
        return arr1[0] * arr2[0] + arr1[1] * arr2[1] + arr1[2] * arr2[2];
}

int main(int argc, char** argv) {

        double PI = 3.141592653589793238462643383279503;
        double c = 299792458.; //[m/s]
        double electron_mass = 9.10938356e-31; // (*kg*)
        double electron_charge = 1.60217662e-19; //(*coulombs*)
        double coulombConst = 8.987551792313e9;
        double HBar = 1.0545718e-34; // (*m^2kg/s*)
        double eps_0 = 8.85418782e-12;
        double Schwinger_E_field= (electron_mass*electron_mass*c*c*c)/(HBar*electron_charge);
    
        // this is the 'random" number we should use to ascertain whether
        // a SFQED event occurs or not. It has been fixed to zero in such a way
        // that the emission is forced. However, this should be a proper random number
        // comprised between 0 and 1
        double rnd_zero = 0.;

        // wave length associated 
        double wave_length; //[m]
        // simulation reference length
        double reference_length;
        // reference angular frequency
        double w_r; //[1/s]
        // time step of the simulation
        double time_step;
        // PIC electric field normalization
        double Enorm;

        // LCFA variables
        double gamma, chi, rnd;
        bool toolkit_init;
        // photon emission
        double emission_rate, phtn_energy;
        // pair production
        double production_rate, electron_energy;
        // values used for theoretical inspections
        double expected_chi;
        
        // LCFA fields and momenta
        double BB[3], EE[3];
        double parts_momenta[2][3];
        double product_momenta[2][3];
        //index used insde while cycling over particles
        int part_num;
        //tmp helper variables
        double tmp_parent, tmp_child;

        // additional variables for BLCFA
        BLCFA_Object * blcfa_entity[2];
        double rnd2;
        double momentum[2][3];
        double pushed_momentum[2][3];
        double Lorentz_F_t_der[2][3];
        double Lorentz_F_tt_der[2][3];
        double delta_blcfa[2];
        double LCFA_threshold;
        bool keep_going_blcfa;

        ///////////////////
        // SET VARIABLES //
        ///////////////////

        // set main simualtions quantities
        time_step = 0.1;
        rnd = 2.1432683788011112e-01;
        rnd2 = 4.6667110607838297e-01;
        wave_length = 0.8e-6; //1e-6; //
        reference_length = wave_length / (2. * PI);
        w_r = c / reference_length;
        Enorm = electron_mass*w_r*c/electron_charge;

        //prepare e.m. fields
        EE[0] = 0.;
        EE[1] = 0.;
        EE[2] = 0.;

        BB[0] = 0.;
        BB[1] = 270.;
        BB[2] = 0.;

        // parts momenta
        parts_momenta[0][0] = -1.2137280955271491e+01;
        parts_momenta[0][1] = 0.0;
        parts_momenta[0][2] = -1.9570511118175094e+04;

        parts_momenta[1][0] = 5.9499124858705343e+00;
        parts_momenta[1][1] = 0.0;
        parts_momenta[1][2] = -1.0000512547779021e+04;

        std::cout << std::scientific << std::setprecision(16);

        std::cout <<"Starting c++ TEST of SFQEDtoolkit!\n";

        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!! SFQEDtoolkit INITIALIZATION !!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        // init SFQEDtoolkit by providing the reference length and the timestep of the simulation
        // char *log_message;
        toolkit_init = SFQED_INIT_ALL_ref_len(reference_length, time_step);
        // toolkit_init = SFQED_INIT_ALL_ref_freq(w_r, time_step);
        std::cout << "SFQEDtoolkit initialization = " << toolkit_init << "\n";
        // std::cout << log_message;

        // std::cout << "photon emission time: " << (1. / (SFQED_INV_COMPTON_rate(1526.4757363249616, 1.0) * 1883651567308853.2)) << '\n'; //
        // std::cout << "pair decay time: " << (1. / (SFQED_BREIT_WHEELER_rate(8242.968976154792, 20.0) * 1883651567308853.2)) << "\n"; //


        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //! LCFA synchrophoton emission !
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        // please notice that in both the photon emission and pair production sections
        // we are cycling over the exact same array of momenta (i.e. the parts_momenta
        // variable defined above). Of course you are free to redefine its components
        // as you please. The same holds true for the magnetic and electric fields.
        // If you do so, mind that the expected_chi variablle is generally no longer
        // reliable, as it was meant to represent the theoretical value of the quantum
        // parameter chi when the particles' momenta are perpendicular to the magnetic
        // field and the electric one is turned off
        std::cout << "\nLCFA SECTION:" << "\n";

        for(part_num = 0; part_num < 2; part_num++){
                
                //compute chi and gamma
                gamma = sqrt(1. + scalar_prod(parts_momenta[part_num], parts_momenta[part_num]));
                chi = compute_chi_with_vectors(gamma, parts_momenta[part_num], EE, BB);

                //theoretical chi expected when the electric field is zero
                //and the magnetic field is perpendicular to the momenta
                //(this must be compared with the chi obtained from the routine)
                expected_chi = gamma * Enorm * sqrt(scalar_prod(BB,BB)) / Schwinger_E_field;
                
                std::cout << "#particle=" << part_num << ", gamma=" << gamma << ", chi=" << chi << " (" << expected_chi << ")\n";
                
                //actual SFQED routines
                emission_rate = SFQED_INV_COMPTON_rate(gamma, chi) * w_r;
                std::cout << "Emission rate = " << emission_rate << '\n';

                if (emission_rate * time_step >= rnd_zero) {
                        phtn_energy = SFQED_LCFA_INV_COMPTON_PHOTON_energy(gamma, chi, rnd);
                        std::cout << "Emitted photon energy = " << phtn_energy << '\n';

                        //compute new particle's momentum
                        SFQED_build_collinear_momentum(phtn_energy, parts_momenta[part_num], product_momenta[part_num]);

                        //output momenta for checking
                        tmp_parent = sqrt(scalar_prod(parts_momenta[part_num],parts_momenta[part_num]));
                        tmp_child = sqrt(scalar_prod(product_momenta[part_num],product_momenta[part_num]));
                        std::cout << "Parent momentum = " 
                                        << parts_momenta[part_num][0] << ' ' << parts_momenta[part_num][1] << ' ' << parts_momenta[part_num][2] << '\n'
                                << "Parent direction = "
                                        << parts_momenta[part_num][0] / tmp_parent << ' ' << parts_momenta[part_num][1] / tmp_parent << ' ' << parts_momenta[part_num][2] / tmp_parent << '\n'
                                << "Child momentum = "
                                        << product_momenta[part_num][0] << ' ' << product_momenta[part_num][1] << ' ' << product_momenta[part_num][2] << '\n'
                                << "Child direction = "
                                        << product_momenta[part_num][0] / tmp_child << ' ' << product_momenta[part_num][1] / tmp_child << ' ' << product_momenta[part_num][2] / tmp_child << '\n';
                }

                std::cout << '\n';
        }

   
        //!!!!!!!!!!!!!!!!!!!
        //! PAIR production !
        //!!!!!!!!!!!!!!!!!!! 
        
        std::cout << "\nPAIRS SECTION:" << "\n";

        // SFQED_BREIT_WHEELER_rate(1000, 20);

        // std::cout << "***********************\n";

        for(part_num = 0; part_num < 2; part_num++){
                //compute chi and gamma
                gamma = sqrt(1. + scalar_prod(parts_momenta[part_num], parts_momenta[part_num]));
                chi = compute_chi_with_vectors(gamma, parts_momenta[part_num], EE, BB);

                //theoretical chi expected when the electric field is zero
                //and the magnetic field is perpendicular to the momenta
                //(this must be compared with the chi obtained from the routine)
                expected_chi = gamma * Enorm * sqrt(scalar_prod(BB,BB)) / Schwinger_E_field;
                
                std::cout << "#particle=" << part_num << ", gamma=" << gamma << ", chi=" << chi << " (" << expected_chi << ")\n";
                
                //actual SFQED routines
                production_rate = SFQED_BREIT_WHEELER_rate(gamma, chi);
                std::cout << "Production rate = " << production_rate << '\n';

                if (production_rate * time_step >= rnd_zero) {
                        electron_energy = SFQED_BREIT_WHEELER_ELECTRON_energy(gamma, chi, rnd);
                        std::cout << "Created electron energy = " << electron_energy << '\n';

                        //compute new particle's momentum
                        SFQED_build_collinear_momentum(electron_energy, parts_momenta[part_num], product_momenta[part_num]);

                        //output momenta for checking
                        tmp_parent = sqrt(scalar_prod(parts_momenta[part_num],parts_momenta[part_num]));
                        tmp_child = sqrt(scalar_prod(product_momenta[part_num],product_momenta[part_num]));
                        std::cout << "Parent momentum = " 
                                        << parts_momenta[part_num][0] << ' ' << parts_momenta[part_num][1] << ' ' << parts_momenta[part_num][2] << '\n'
                                << "Parent direction = "
                                        << parts_momenta[part_num][0] / tmp_parent << ' ' << parts_momenta[part_num][1] / tmp_parent << ' ' << parts_momenta[part_num][2] / tmp_parent << '\n'
                                << "Child momentum = "
                                        << product_momenta[part_num][0] << ' ' << product_momenta[part_num][1] << ' ' << product_momenta[part_num][2] << '\n'
                                << "Child direction = "
                                        << product_momenta[part_num][0] / tmp_child << ' ' << product_momenta[part_num][1] / tmp_child << ' ' << product_momenta[part_num][2] / tmp_child << '\n';
                }

                std::cout << '\n';
      
        }

        //!!!!!!!!!!!!!!!!!!!!!!!
        //! BLCFA synchrophoton !
        //!!!!!!!!!!!!!!!!!!!!!!!

        //! set BLCFA quantities
        blcfa_entity[0] = SFQED_CREATE_BLCFA_OBJECT();
        blcfa_entity[1] = SFQED_CREATE_BLCFA_OBJECT();

        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!! PRELIMINARY BLCFA ACTIONS !!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        // when using the blcfa approach we need some preliminary 
        // actions to properly set the lorentz force derivatives
        // and chi and gamma. We will perform the same action 3
        // times, since we are interested in getting the second time
        // derivatives of the force

        ///////// STEP 1 /////////

        //old momentum
        momentum[0][0] = 2.0923765116618465e+00;
        momentum[0][1] = 0.0000000000000000e+00;
        momentum[0][2] = -1.9570512944081507e+04;
        //new pushed momentum
        pushed_momentum[0][0] = 5.9499124858705343e+00;
        pushed_momentum[0][1] =  0.0000000000000000e+00;
        pushed_momentum[0][2] = -1.9570512547779021e+04;
        //update blcfa obj
        keep_going_blcfa = SFQED_BLCFA_OBJECT_update(blcfa_entity[0], pushed_momentum[0], momentum[0], delta_blcfa[0], gamma, chi);

        //old momentum
        momentum[1][0] = -1.7718145323355337e+01;
        momentum[1][1] = 0.0000000000000000e+00;
        momentum[1][2] = -1.9570508989729224e+04;
        //new pushed momentum
        pushed_momentum[1][0] = -1.6542717780968854e+01;
        pushed_momentum[1][1] =  0.0000000000000000e+00;
        pushed_momentum[1][2] = -1.9570509504165970e+04;
        //update blcfa obj
        keep_going_blcfa = SFQED_BLCFA_OBJECT_update(blcfa_entity[1], pushed_momentum[1], momentum[1], delta_blcfa[1], gamma, chi);

        ///////// STEP 2 /////////

        //old momentum                           
        momentum[0][0] = 5.9499124858705343e+00;
        momentum[0][1] = 0.0000000000000000e+00;
        momentum[0][2] = -1.9570512547779021e+04;
        //new pushed momentum
        pushed_momentum[0][0] = 9.6034048452081251e+00;
        pushed_momentum[0][1] = 0.0000000000000000e+00;
        pushed_momentum[0][2] = -1.9570511821892545e+04;
        //update blcfa obj
        keep_going_blcfa = SFQED_BLCFA_OBJECT_update(blcfa_entity[0], pushed_momentum[0], momentum[0], delta_blcfa[0], gamma, chi);

        //old momentum
        momentum[1][0] = -1.6542717780968854e+01;
        momentum[1][1] = 0.0000000000000000e+00;
        momentum[1][2] = -1.9570509504165970e+04;
        //new pushed momentum
        pushed_momentum[1][0] = -1.4660387531605418e+01;
        pushed_momentum[1][1] =  0.0000000000000000e+00;
        pushed_momentum[1][2] = -1.9570510254460005e+04;
        //update blcfa obj
        keep_going_blcfa = SFQED_BLCFA_OBJECT_update(blcfa_entity[1], pushed_momentum[1], momentum[1], delta_blcfa[1], gamma, chi);

        ///////// STEP 3 /////////
        //in this last step we do not update the blcfa objects,
        //as this is done inside the particle loop

        //old momentum 
        momentum[0][0] = 9.6034048452081251e+00;
        momentum[0][1] = 0.0000000000000000e+00;
        momentum[0][2] = -1.9570511821892545e+04;
        //new pushed momentum
        pushed_momentum[0][0] = 1.2895899733826127e+01;
        pushed_momentum[0][1] = 0.0000000000000000e+00;
        pushed_momentum[0][2] = -1.9570510875586107e+04;

        //old momentum
        momentum[1][0] = -1.4660387531605418e+01;
        momentum[1][1] =  0.0000000000000000e+00;
        momentum[1][2] = -1.9570510254460005e+04;
        //new pushed momentum
        pushed_momentum[1][0] = -1.2137280955271491e+01;
        pushed_momentum[1][1] =  0.0000000000000000e+00;
        pushed_momentum[1][2] = -1.9570511118175094e+04;


        ////////////  ACTUAL BLCFA LOOP  ////////////
        // at this point everything should be fine and the actual blcfa section can now begin
        std::cout << "\nBLCFA SECTION:\n";
        
        for(part_num = 0; part_num < 2; part_num++){

                // update the blcfa objects. Remember that this update operation must be 
                // done at each timestep of the simulation for every particle. The update
                // requires you to know the initial momentum possessed by the given particle,
                // together with the pushed momentum, which would result by normally applying
                // the fields and forces on the particle. The determination of these quantities,
                // which should be performed in this same loop has been here omitted, and replaced
                // with the 3 preliminary operations above
                keep_going_blcfa = SFQED_BLCFA_OBJECT_update(blcfa_entity[part_num], pushed_momentum[part_num], momentum[part_num],
                                                                delta_blcfa[part_num], gamma, chi);

                std::cout << "#particle = " << part_num << '\n';
                std::cout << "chi = " << chi << '\n';
                std::cout << "gamma = " << gamma << '\n';
                // std::cout << "first derivative = " << Lorentz_F_t_der[part_num][0] << ' ' << Lorentz_F_t_der[part_num][1] << ' ' << Lorentz_F_t_der[part_num][2] << '\n';
                // std::cout << "second derivative = " << Lorentz_F_tt_der[part_num][0] << ' ' << Lorentz_F_tt_der[part_num][1] << ' ' << Lorentz_F_tt_der[part_num][2] << '\n';

                emission_rate = SFQED_INV_COMPTON_rate(gamma, chi);
                std::cout << "Emission rate = " << emission_rate << '\n';

                if (keep_going_blcfa && emission_rate * time_step >= rnd_zero) {
                        
                        LCFA_threshold = SFQED_BLCFA_INV_COMPTON_PHOTON_threshold(blcfa_entity[part_num], delta_blcfa[part_num], gamma, chi);
                        std::cout << "LCFA threshold = " << LCFA_threshold << '\n';

                        phtn_energy = SFQED_BLCFA_INV_COMPTON_PHOTON_energy(LCFA_threshold, gamma, chi, rnd, rnd2);
                        std::cout << "Emitted photon energy = " << phtn_energy << '\n';

                        //now you should check if phtn_energy > 0. In that case the emission occurs, otherwise
                        // nothing happens

                        //compute new particle's momentum
                        SFQED_build_collinear_momentum(phtn_energy, pushed_momentum[part_num], product_momenta[part_num]);

                        //output momenta for checking
                        tmp_parent = sqrt(scalar_prod(pushed_momentum[part_num],pushed_momentum[part_num]));
                        tmp_child = sqrt(scalar_prod(product_momenta[part_num],product_momenta[part_num]));
                        std::cout << "Parent momentum = " 
                                        << pushed_momentum[part_num][0] << ' ' << pushed_momentum[part_num][1] << ' ' << pushed_momentum[part_num][2] << '\n'
                                << "Parent direction = "
                                        << pushed_momentum[part_num][0] / tmp_parent << ' ' << pushed_momentum[part_num][1] / tmp_parent << ' ' << pushed_momentum[part_num][2] / tmp_parent << '\n'
                                << "Child momentum = "
                                        << product_momenta[part_num][0] << ' ' << product_momenta[part_num][1] << ' ' << product_momenta[part_num][2] << '\n'
                                << "Child direction = "
                                        << product_momenta[part_num][0] / tmp_child << ' ' << product_momenta[part_num][1] / tmp_child << ' ' << product_momenta[part_num][2] / tmp_child << '\n';
                }

                std::cout << '\n';

        }
   
   
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   //!! SFQEDtoolkit FINALIZATION !!
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SFQED_FINALIZE_ALL();

}