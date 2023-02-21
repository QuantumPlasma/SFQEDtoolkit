/***********************************************************************
    Copyright (c) 2000-2021,
    QuantumPlasma team, Max Planck Institut fur Kernphysik
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:

    * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.  

    * Redistributions in binary form must reproduce the above
        copyright notice, this list of conditions and the following
        disclaimer in the documentation and/or other materials provided
        with the distribution.  

    * Neither the name of the copyright holder nor the names of its
        contributors may be used to endorse or promote products derived
        from this software without specific prior written permission.

    * The software using and/or implementing this library must include
        a citation to [] in the official documentation.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
    COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
    INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
    OF THE POSSIBILITY OF SUCH DAMAGE.
*************************************************************************/

#include "SFQEDtoolkit_Interface.hpp"
#include "SFQED_Processes.h"
#include "BLCFA_Object.h"
#include "Breit_Wheeler_Pair_Production.h"
#include "Nonlinear_Inverse_Compton_Scattering.h"

#include "SFQEDtoolkit_BaseInterface_Functions.h"

#include <iostream>
#include <stdio.h>
#include <string.h>


/**************************/
/* SFQED process instance */
/**************************/
//here we define the (pointers to the) objects that
// a simulation needs in order to apply the Nonlinear Inverse
// Compton Scattering and Breit-Wheeler pair production mechanisms
//The coefficients stored inside must be initialized through
// the init functions
SFQED_Processes* procs;
Breit_Wheeler_Pair_Production* bw_pp;
Nonlinear_Inverse_Compton_Scattering* nics;

//these two variables represent the conversion constants 
// which the rates returned by the toolkit have to be multiplied by
// in order to properly express them in code units.
double phtn_emission_rate_coefficient, pair_creation_rate_coefficient;

void SFQED_INIT_PROCS(){
        //initialization of the pointers
        procs = new SFQED_Processes();
        bw_pp = new Breit_Wheeler_Pair_Production();
        nics = new Nonlinear_Inverse_Compton_Scattering();
}

/**********************************/
/* SIMULATION INITIALIZER SECTION */
/**********************************/

//files to check to load the phtn emission coefficients
const int phtn_file_names_size = 20;
const char *phtn_file_names[20] = {
        //photon emission rate
        "phtn_rate/phtn_emission_rate_0-2.txt",
        "phtn_rate/phtn_emission_rate_2-20.txt",
        "phtn_rate/phtn_emission_rate_20-80.txt",
        "phtn_rate/phtn_emission_rate_80-600.txt",
        "phtn_rate/phtn_emission_rate_600-2000.txt",
        //photon partial emission rate
        "phtn_prtl_rate/phtn_partial_rate_0-2.txt",
        "phtn_prtl_rate/phtn_partial_rate_2-20.txt",
        "phtn_prtl_rate/phtn_partial_rate_20-80.txt",
        "phtn_prtl_rate/phtn_partial_rate_80-600.txt",
        "phtn_prtl_rate/phtn_partial_rate_600-2000.txt",
        //photon energy in w
        "w/phtn_w_nrg_0-2.txt",
        "w/phtn_w_nrg_2-20.txt",
        "w/phtn_w_nrg_20-80.txt",
        "w/phtn_w_nrg_80-600.txt",
        "w/phtn_w_nrg_600-2000.txt",
        //photon energy in 1/w
        "1_over_w/phtn_1_over_w_nrg_0-2.txt",
        "1_over_w/phtn_1_over_w_nrg_2-20.txt",
        "1_over_w/phtn_1_over_w_nrg_20-80.txt",
        "1_over_w/phtn_1_over_w_nrg_80-600.txt",
        "1_over_w/phtn_1_over_w_nrg_600-2000.txt"
};

//files to check to load the pair production coefficients
const int pair_file_names_size = 23;
const char *pair_file_names[23] = {
        //pair creation rate
        "pair_rate/pair_production_rate_024-04.txt",
        "pair_rate/pair_production_rate_0-2.txt",
        "pair_rate/pair_production_rate_2-20.txt",
        "pair_rate/pair_production_rate_20-80.txt",
        "pair_rate/pair_production_rate_80-600.txt",
        "pair_rate/pair_production_rate_600-2000.txt",
        //pair creation partial rate
        "pair_prtl_rate/pair_partial_rate_0-2.txt",
        "pair_prtl_rate/pair_partial_rate_2-20.txt",
        "pair_prtl_rate/pair_partial_rate_20-80.txt",
        "pair_prtl_rate/pair_partial_rate_80-600.txt",
        "pair_prtl_rate/pair_partial_rate_600-2000.txt",
        //photon energy in v
        "v/pair_nrgs_v_001-03.txt",
        "v/pair_nrgs_v_0-2.txt",
        "v/pair_nrgs_v_2-20.txt",
        "v/pair_nrgs_v_20-80.txt",
        "v/pair_nrgs_v_80-600.txt",
        "v/pair_nrgs_v_600-2000.txt",
        //photon energy in v high part
        "v/pair_nrgs_v_high_001-03.txt",
        "v/pair_nrgs_v_high_0-2.txt",
        "v/pair_nrgs_v_high_2-20.txt",
        "v/pair_nrgs_v_high_20-80.txt",
        "v/pair_nrgs_v_high_80-600.txt",
        "v/pair_nrgs_v_high_600-2000.txt"
};

const int blcfa_file_names_size = 5;
const char *blcfa_file_names[5] = {
        "dP_over_dt_dw/phtn_diff_crss_sctn_0-2.txt",
        "dP_over_dt_dw/phtn_diff_crss_sctn_2-20.txt",
        "dP_over_dt_dw/phtn_diff_crss_sctn_20-80.txt",
        "dP_over_dt_dw/phtn_diff_crss_sctn_80-600.txt",
        "dP_over_dt_dw/phtn_diff_crss_sctn_600-2000.txt"
};



/****  COMPLETE SETTERS AND FINALIZERS  ****/

// the log_message argument will be reinitialized inside the function
bool SFQED_INIT_ALL_ref_len(const double& ref_len, const double& ts){
        //initialization of the pointers
        procs = new SFQED_Processes();
        bw_pp = new Breit_Wheeler_Pair_Production();
        nics = new Nonlinear_Inverse_Compton_Scattering();

        /////////////////////////////////////////////////////////////
        // //old//old way with environment variable
        // //retrieve the path to coefficients
        // const char *variable_name = "SFQED_TOOLKIT_USER"; // use string literals to initialize a const pointer to char
        // const char * val = std::getenv(variable_name);
        // if ( val == nullptr ) { // invalid to assign nullptr to std::string
        //         std::cout << "Environment variable " << variable_name
        //                 << "not found! Not able to load the coefficients.\n";
        //         // return false;
        // }
        // std::string str_path = std::string(val) + std::string("/coefficients/");
        // std::cout << "Loading SFQED coefficients from " << str_path << "\n";

        // // //new way (the coefficients folder must be copied to where the executable is launched)
        // // std::cout << "Loading SFQED coefficients from root directory\n";
        // // std::string str_path("./coefficients/");
        /////////////////////////////////////////////////////////////

        //retireve path to coefficients
        char *log_message;
        std::string message_to_add;
        std::string str_path = find_path(message_to_add);

        //check the existence of the coefficients' files
        if( !(file_checker(str_path, phtn_file_names, phtn_file_names_size, message_to_add)
                && file_checker(str_path, pair_file_names, pair_file_names_size, message_to_add)
                && file_checker(str_path, blcfa_file_names, blcfa_file_names_size, message_to_add)) ){
                
                log_message = new char[message_to_add.length() + 1];
                strcpy(log_message, message_to_add.c_str());
                std::cout << log_message;
                return false;
        }

        //copy the message to the newly created
        //[keep this part for when the log message will be passed outside]
        // log_message = new char[message_to_add.length() + 1];
        // strcpy(log_message, message_to_add.c_str());
        
        //load coefficients
        bw_pp->SFQED_init_PAIR_creation(str_path);
        nics->SFQED_init_PHTN_emission(str_path);
        nics->SFQED_init_PHTN_emission_BLCFA(str_path);
        
        //initialize reference length and time
        procs->SFQED_set_reference_length(ref_len);
        //this functions must be called only after SFQED_set_reference_length or SFQED_set_reference_angular_frequency
        BLCFA_Object::SFQED_BLCFA_set_time_constants(ts, procs->get_normalized_compton_time());

        //initialize phtn creation and pair production rates' coefficients
        phtn_emission_rate_coefficient = procs->get_phtn_rate_coefficient();
        pair_creation_rate_coefficient = procs->get_pair_rate_coefficient();
        
        //return true when everythin is done
        return true;
}

bool SFQED_INIT_ALL_ref_freq(const double& ref_freq, const double& ts){
        //initialization of the pointers
        procs = new SFQED_Processes();
        bw_pp = new Breit_Wheeler_Pair_Production();
        nics = new Nonlinear_Inverse_Compton_Scattering();
        
        //////////////////////////////////////////////////////////////////
        // //old way with environment variable
        // //retrieve the path to coefficients
        // const char *variable_name = "SFQED_TOOLKIT_USER"; // use string literals to initialize a const pointer to char
        // const char * val = std::getenv(variable_name);
        // if ( val == nullptr ) { // invalid to assign nullptr to std::string
        //         std::cout << "Environment variable " << variable_name
        //                 << "not found! Not able to load the coefficients.\n";
        //         // return false;
        // }
        // std::string str_path = std::string(val) + std::string("/coefficients/");
        // std::cout << "Loading SFQED coefficients from " << str_path << "\n";

        // // //new way (copy coefficients folder to where the executable is launched)
        // // std::cout << "Loading SFQED coefficients from root directory\n";
        // // std::string str_path("./coefficients/");
        ///////////////////////////////////////////////////////////////
        
        //retireve path to coefficients
        char *log_message;
        std::string message_to_add;
        std::string str_path = find_path(message_to_add);

        //check the existence of the coefficients' files
        if( !(file_checker(str_path, phtn_file_names, phtn_file_names_size, message_to_add)
                && file_checker(str_path, pair_file_names, pair_file_names_size, message_to_add)
                && file_checker(str_path, blcfa_file_names, blcfa_file_names_size, message_to_add)) ){
                
                log_message = new char[message_to_add.length() + 1];
                strcpy(log_message, message_to_add.c_str());
                std::cout << log_message;
                return false;
        }

        //copy the message to the newly created
        //[keep this part for when the log message will be passed outside]
        // log_message = new char[message_to_add.length() + 1];
        // strcpy(log_message, message_to_add.c_str());
        
        //load coefficients
        bw_pp->SFQED_init_PAIR_creation(str_path);
        nics->SFQED_init_PHTN_emission(str_path);
        nics->SFQED_init_PHTN_emission_BLCFA(str_path);
        
        //initialize reference freq and time
        procs->SFQED_set_reference_angular_frequency(ref_freq);
        //this functions must be called only after SFQED_set_reference_length or SFQED_set_reference_angular_frequency
        BLCFA_Object::SFQED_BLCFA_set_time_constants(ts, procs->get_normalized_compton_time());

        //initialize phtn creation and pair production rates' coefficients
        phtn_emission_rate_coefficient = procs->get_phtn_rate_coefficient();
        pair_creation_rate_coefficient = procs->get_pair_rate_coefficient();
        
        //return true when everythin is done
        return true;
}

void SFQED_INIT_ALL_debug(){
        //initialization of the pointers
        procs = new SFQED_Processes();
        bw_pp = new Breit_Wheeler_Pair_Production();
        nics = new Nonlinear_Inverse_Compton_Scattering();

        //////////////////////////////////////////////////////////////////// 
        // //old way with environment variable
        // //retrieve the path to coefficients
        // const char *variable_name = "SFQED_TOOLKIT_USER"; // use string literals to initialize a const pointer to char
        // const char * val = std::getenv(variable_name);
        // if ( val == nullptr ) { // invalid to assign nullptr to std::string
        //         std::cout << "Environment variable " << variable_name
        //                 << "not found! Not able to load the coefficients.\n";
        //         // return false;
        // }
        // std::string str_path = std::string(val) + std::string("/coefficients/");
        // std::cout << "Loading SFQED coefficients from " << str_path << "\n";

        // // //new way (copy coefficients folder to where the executable is launched)
        // // std::cout << "Loading SFQED coefficients from root directory\n";
        // // std::string str_path("./coefficients/");
        ////////////////////////////////////////////////////////////////////

        //retireve path to coefficients
        char *log_message;
        std::string message_to_add;
        std::string str_path = find_path(message_to_add);

        //check the existence of the coefficients' files
        if( !(file_checker(str_path, phtn_file_names, phtn_file_names_size, message_to_add)
                && file_checker(str_path, pair_file_names, pair_file_names_size, message_to_add)
                && file_checker(str_path, blcfa_file_names, blcfa_file_names_size, message_to_add)) ){
                
                log_message = new char[message_to_add.length() + 1];
                strcpy(log_message, message_to_add.c_str());
                std::cout << log_message;
                return;
        }

        //copy the message to the newly created
        //[keep this part for when the log message will be passed outside]
        // log_message = new char[message_to_add.length() + 1];
        // strcpy(log_message, message_to_add.c_str());
        
        //load coefficients
        bw_pp->SFQED_init_PAIR_creation(str_path);
        nics->SFQED_init_PHTN_emission(str_path);
        nics->SFQED_init_PHTN_emission_BLCFA(str_path);
        
        //initialize reference freq and time
        procs->SFQED_set_all_to_one();

        //initialize phtn creation and pair production rates' coefficients
        phtn_emission_rate_coefficient = procs->get_phtn_rate_coefficient();
        pair_creation_rate_coefficient = procs->get_pair_rate_coefficient();
}

void SFQED_FINALIZE_ALL(){
        bw_pp->SFQED_finalize_PAIR_creation();
        nics->SFQED_finalize_PHTN_emission();
        nics->SFQED_finalize_PHTN_emission_BLCFA();

        delete bw_pp;
        delete nics;
}

/****  PARTIAL SETTERS  ****/

//set all the simulation's quantities based on
//the value of the given reference length
void SFQED_set_ref_len(const double& ref_len){
        procs->SFQED_set_reference_length(ref_len);

        //initialize phtn creation and pair production rates' coefficients
        phtn_emission_rate_coefficient = procs->get_phtn_rate_coefficient();
        pair_creation_rate_coefficient = procs->get_pair_rate_coefficient();
}
        
void SFQED_set_ref_freq(const double& ref_freq){
        procs->SFQED_set_reference_angular_frequency(ref_freq);

        //initialize phtn creation and pair production rates' coefficients
        phtn_emission_rate_coefficient = procs->get_phtn_rate_coefficient();
        pair_creation_rate_coefficient = procs->get_pair_rate_coefficient();
}

//useful for debug
void SFQED_set_for_debug(){
        procs->SFQED_set_all_to_one();

        //initialize phtn creation and pair production rates' coefficients
        phtn_emission_rate_coefficient = procs->get_phtn_rate_coefficient();
        pair_creation_rate_coefficient = procs->get_pair_rate_coefficient();
}

void SFQED_set_sim_tstep(const double& ts){
        //this functions must be called only after SFQED_set_reference_length or SFQED_set_reference_angular_frequency
        BLCFA_Object::SFQED_BLCFA_set_time_constants(ts, procs->get_normalized_compton_time());
}


/*****************************************************/
/* OBSOLETE coefficients initializers and finalizers */
/*****************************************************/
// /* PHOTON */
// void SFQED_init_INV_COMPTON(const char* const path){
//         std::string str_path(path);
//         nics->SFQED_init_PHTN_emission(str_path);
// }

// //obsolete functions: the nics pointer is not deleted!!!
// void SFQED_finalize_INV_COMPTON(){
//         nics->SFQED_finalize_PHTN_emission();
// }

// /* PHOTON BLCFA */
// void SFQED_init_INV_COMPTON_BLCFA(const char* const path){
//         std::string str_path(path);
//         nics->SFQED_init_PHTN_emission_BLCFA(str_path);
// }

// //obsolete functions: the nics pointer is not deleted!!!
// void SFQED_finalize_INV_COMPTON_BLCFA(){
//         nics->SFQED_finalize_PHTN_emission_BLCFA();
// }

// /* PAIRS */
// void SFQED_init_BREIT_WHEELER(const char* const path){
//         std::string str_path(path);
//         bw_pp->SFQED_init_PAIR_creation(str_path);
// }

// //obsolete functions: the bw_pp pointer is not deleted!!!
// void SFQED_finalize_BREIT_WHEELER(){
//         bw_pp->SFQED_finalize_PAIR_creation();
// }

/*********************************/
/* QUANTUM PARAMETER CALCULATION */
/*********************************/
double compute_chi_with_vectors(const double& gamma, const double p_in[3], const double EE[3], const double BB[3]){
        return procs->compute_quantum_param(gamma, p_in, EE, BB);
}

double compute_chi_with_components(const double &gamma,
                                                                const double &p_in_x,
                                                                const double &p_in_y,
                                                                const double &p_in_z,
                                                                const double &EE_x,
                                                                const double &EE_y,
                                                                const double &EE_z,
                                                                const double &BB_x,
                                                                const double &BB_y,
                                                                const double &BB_z){

        return procs->compute_quantum_param(gamma, p_in_x, p_in_y, p_in_z, EE_x, EE_y, EE_z, BB_x, BB_y, BB_z);                                                                            
}


/*****************/
/* BUILD MOMENTA */
/*****************/
void SFQED_build_collinear_momentum(const double &gamma_out, const double p_in[3], double p_out[3]){
        SFQED_collinear_momentum(gamma_out, p_in, p_out);
}

/********************************/
/* PHOTON EMISSION SECTION LCFA */
/********************************/
double SFQED_INV_COMPTON_rate(const double &gamma, const double &chi){
        return phtn_emission_rate_coefficient * nics->SFQED_PHTN_emission_rate(gamma, chi);
}

double SFQED_LCFA_INV_COMPTON_PHOTON_energy(const double &gamma, const double &chi, const double &rnd){
        return nics->SFQED_LCFA_emitted_photon_energy(gamma, chi, rnd);
}

/*********************************/
/* PHOTON EMISSION SECTION BLCFA */
/*********************************/

BLCFA_Object* SFQED_CREATE_BLCFA_OBJECT(){
        return new BLCFA_Object();
}

void SFQED_FINALIZE_BLCFA_OBJECT(BLCFA_Object* entity){
        delete entity;
        entity = NULL;
}

//you should add a function that creates an array (a std::vector is better) of BLCFA_Object

bool SFQED_BLCFA_OBJECT_update(BLCFA_Object* entity,
                                                const double* const pushed_momentum,
                                                const double* const momentum,
                                                double& delta,
                                                double& part_gamma,
                                                double& part_chi){
                                                        
        return entity->SFQED_BLCFA_update_entities_quantities(pushed_momentum, momentum,
                                                delta, part_gamma, part_chi);
}


double SFQED_BLCFA_INV_COMPTON_PHOTON_threshold(const BLCFA_Object* const entity,
                                                        const double& delta,
                                                        const double& part_gamma,
                                                        const double& part_chi){

        return entity->SFQED_BLCFA_find_energy_threshold(delta, part_gamma, part_chi);
}


double SFQED_BLCFA_INV_COMPTON_PHOTON_energy(const double& LCFA_limit,
                                                        const double& gamma, 
                                                        const double& chi, 
                                                        const double& rnd, 
                                                        const double& rnd2){

        return nics->SFQED_BLCFA_emitted_photon_energy(LCFA_limit, gamma, chi, rnd, rnd2);
}

//debug (this functions could be deleted)*********************
void SFQED_BLCFA_printALL(BLCFA_Object* entity){
        std::cout << entity->is_just_created() << "\n";
        double *tmp = entity->FL_old();
        std::cout << tmp[0] << ' ' << tmp[1] << ' ' << tmp[2] << "\n";
        tmp = entity->deltaFL_old();
        std::cout << tmp[0] << ' ' << tmp[1] << ' ' << tmp[2] << "\n";
}

void SFQED_BLCFA_get_F_old(BLCFA_Object* entity, double* F_old){
        double *tmp = entity->FL_old();
        F_old[0] = tmp[0];
        F_old[1] = tmp[1];
        F_old[2] = tmp[2];
        std::cout << F_old[0] << ' ' << F_old[1] << ' ' << F_old[2] << "\n";
}

void SFQED_BLCFA_update_F_old(BLCFA_Object* entity, const double* const F_old_new){
        std::cout << F_old_new[0] << ' ' << F_old_new[1] << ' ' << F_old_new[2] << "\n";
        entity->update_force(F_old_new);
        double *arr = entity->FL_old();
        std::cout << arr[0] << ' ' << arr[1] << ' ' << arr[2] << "\n";
}

void SFQED_BLCFA_get_delta_F_old(BLCFA_Object* entity, double* delta_F_old){
        double *tmp = entity->deltaFL_old();
        // double *tmp = entity->FL_old();
        delta_F_old[0] = tmp[0];
        delta_F_old[1] = tmp[1];
        delta_F_old[2] = tmp[2];
        std::cout << delta_F_old[0] << ' ' << delta_F_old[1] << ' ' << delta_F_old[2] << "\n";
}
//**************************************************************************************

/***************************/
/* PAIR PRODUCTION SECTION */
/***************************/
double SFQED_BREIT_WHEELER_rate(const double &gamma, const double &chi){
        return pair_creation_rate_coefficient * bw_pp->SFQED_PAIR_creation_rate(gamma, chi);
}

double SFQED_BREIT_WHEELER_rate_fast(const double &gamma, const double &chi){
        return pair_creation_rate_coefficient * bw_pp->SFQED_PAIR_creation_rate_fast(gamma, chi);
}

double SFQED_BREIT_WHEELER_ELECTRON_energy(const double &nrg, const double &chi, const double &rnd){
        return bw_pp->SFQED_PAIR_created_electron_energy(nrg, chi, rnd);
}

double SFQED_BREIT_WHEELER_ELECTRON_energy_fast(const double &nrg, const double &chi, const double &rnd){
        return bw_pp->SFQED_PAIR_created_electron_energy_fast(nrg, chi, rnd);
}
