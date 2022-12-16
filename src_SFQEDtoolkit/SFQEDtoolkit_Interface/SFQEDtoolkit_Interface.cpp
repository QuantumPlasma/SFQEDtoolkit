#include "SFQEDtoolkit_Interface.hpp"
#include "SFQED_Processes.h"
#include "BLCFA_Object.h"
#include "Breit_Wheeler_Pair_Production.h"
#include "Nonlinear_Inverse_Compton_Scattering.h"

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

/****  COMPLETE SETTERS AND FINALIZERS  ****/

//remember to check the presence of the files you want to read!!!

bool SFQED_INIT_ALL_ref_len(const double& ref_len, const double& ts){
        //initialization of the pointers
        procs = new SFQED_Processes();
        bw_pp = new Breit_Wheeler_Pair_Production();
        nics = new Nonlinear_Inverse_Compton_Scattering();

        //old//old way with environment variable
        //retrieve the path to coefficients
        const char *variable_name = "SFQED_TOOLKIT_USER"; // use string literals to initialize a const pointer to char
        const char * val = std::getenv(variable_name);
        if ( val == nullptr ) { // invalid to assign nullptr to std::string
                std::cout << "Environment variable " << variable_name
                        << "not found! Not able to load the coefficients.\n";
                // return false;
        }
        std::string str_path = std::string(val) + std::string("/coefficients/");
        std::cout << "Loading SFQED coefficients from " << str_path << "\n";

        // //new way (the coefficients folder must be copied to where the executable is launched)
        // std::cout << "Loading SFQED coefficients from root directory\n";
        // std::string str_path("./coefficients/");
        
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
        
        //old way with environment variable
        //retrieve the path to coefficients
        const char *variable_name = "SFQED_TOOLKIT_USER"; // use string literals to initialize a const pointer to char
        const char * val = std::getenv(variable_name);
        if ( val == nullptr ) { // invalid to assign nullptr to std::string
                std::cout << "Environment variable " << variable_name
                        << "not found! Not able to load the coefficients.\n";
                // return false;
        }
        std::string str_path = std::string(val) + std::string("/coefficients/");
        std::cout << "Loading SFQED coefficients from " << str_path << "\n";

        // //new way (copy coefficients folder to where the executable is launched)
        // std::cout << "Loading SFQED coefficients from root directory\n";
        // std::string str_path("./coefficients/");
        
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
        
        //old way with environment variable
        //retrieve the path to coefficients
        const char *variable_name = "SFQED_TOOLKIT_USER"; // use string literals to initialize a const pointer to char
        const char * val = std::getenv(variable_name);
        if ( val == nullptr ) { // invalid to assign nullptr to std::string
                std::cout << "Environment variable " << variable_name
                        << "not found! Not able to load the coefficients.\n";
                // return false;
        }
        std::string str_path = std::string(val) + std::string("/coefficients/");
        std::cout << "Loading SFQED coefficients from " << str_path << "\n";

        // //new way (copy coefficients folder to where the executable is launched)
        // std::cout << "Loading SFQED coefficients from root directory\n";
        // std::string str_path("./coefficients/");
        
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
/* PHOTON */
void SFQED_init_INV_COMPTON(const char* const path){
        std::string str_path(path);
        nics->SFQED_init_PHTN_emission(str_path);
}

//obsolete functions: the nics pointer is not deleted!!!
void SFQED_finalize_INV_COMPTON(){
        nics->SFQED_finalize_PHTN_emission();
}

/* PHOTON BLCFA */
void SFQED_init_INV_COMPTON_BLCFA(const char* const path){
        std::string str_path(path);
        nics->SFQED_init_PHTN_emission_BLCFA(str_path);
}

//obsolete functions: the nics pointer is not deleted!!!
void SFQED_finalize_INV_COMPTON_BLCFA(){
        nics->SFQED_finalize_PHTN_emission_BLCFA();
}

/* PAIRS */
void SFQED_init_BREIT_WHEELER(const char* const path){
        std::string str_path(path);
        bw_pp->SFQED_init_PAIR_creation(str_path);
}

//obsolete functions: the bw_pp pointer is not deleted!!!
void SFQED_finalize_BREIT_WHEELER(){
        bw_pp->SFQED_finalize_PAIR_creation();
}

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
