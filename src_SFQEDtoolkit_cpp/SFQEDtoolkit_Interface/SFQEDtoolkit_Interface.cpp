#include "SFQEDtoolkit_Interface.hpp"
#include "SFQED_Processes.h"

#include <iostream>
#include <stdio.h>
#include <string.h>


/**************************/
/* SFQED process instance */
/**************************/
SFQED_Processes procs;

void SFQED_INIT_PROCS(){
        procs = SFQED_Processes();
}

/**********************************/
/* SIMULATION INITIALIZER SECTION */
/**********************************/

/****  COMPLETE SETTERS AND FINALIZERS  ****/

bool SFQED_INIT_ALL_ref_len(const double& ref_len, const double& ts){
        //initialize process instance
        procs = SFQED_Processes();

        //retrieve the path to coefficients
        const char *variable_name = "SFQED_TOOLKIT_USER"; // use string literals to initialize a const pointer to char
        const char * val = std::getenv(variable_name);
        if ( val == nullptr ) { // invalid to assign nullptr to std::string
                std::cout << "Environment variable " << variable_name
                        << "not found! Not able to load the coefficients.\n";
                return false;
        }
        std::string str_path = std::string(val) + std::string("/coefficients/");
        std::cout << "Loading SFQED coefficients from " << str_path << "\n";
        
        //load coefficients
        procs.SFQED_init_PHTN_emission(str_path);
        procs.SFQED_init_PAIR_creation(str_path);
        
        //initialize reference length and time
        procs.SFQED_set_reference_length(ref_len);
        procs.SFQED_set_time_step(ts);
        
        //return true when everythin is done
        return true;
}

bool SFQED_INIT_ALL_ref_freq(const double& ref_freq, const double& ts){
        //initialize process instance
        procs = SFQED_Processes();
        
        //retrieve the path to coefficients
        const char *variable_name = "SFQED_TOOLKIT_USER"; // use string literals to initialize a const pointer to char
        const char * val = std::getenv(variable_name);
        if ( val == nullptr ) { // invalid to assign nullptr to std::string
                std::cout << "Environment variable " << variable_name
                        << "not found! Not able to load the coefficients.\n";
                return false;
        }
        std::string str_path = std::string(val) + std::string("/coefficients/");
        std::cout << "Loading SFQED coefficients from " << str_path << "\n";
        
        //load coefficients
        procs.SFQED_init_PHTN_emission(str_path);
        procs.SFQED_init_PAIR_creation(str_path);
        
        //initialize reference freq and time
        procs.SFQED_set_reference_angular_frequency(ref_freq);
        procs.SFQED_set_time_step(ts);
        
        //return true when everythin is done
        return true;
}

void SFQED_FINALIZE_ALL(){
        procs.SFQED_finalize_PHTN_emission();
        procs.SFQED_finalize_PAIR_creation();
}

/****  PARTIAL SETTERS  ****/

//set all the simulation's quantities based on
//the value of the given reference length
void SFQED_set_ref_len(const double& ref_len){
        procs.SFQED_set_reference_length(ref_len);
}
        
void SFQED_set_ref_freq(const double& ref_freq){
        procs.SFQED_set_reference_angular_frequency(ref_freq);
}

//useful for debug
void SFQED_set_for_debug(){
        procs.SFQED_set_all_to_one();
}
        
void SFQED_set_sim_tstep(const double& ts){
        procs.SFQED_set_time_step(ts);
}


/********************************************/
/* coefficients initializers and finalizers */
/********************************************/
/* PHOTON */
void SFQED_init_INV_COMPTON(const char* const path){
        std::string str_path(path);
        procs.SFQED_init_PHTN_emission(str_path);
}

void SFQED_finalize_INV_COMPTON(){
        procs.SFQED_finalize_PHTN_emission();
}

/* PAIRS */
void SFQED_init_BREIT_WHEELER(const char* const path){
        std::string str_path(path);
        procs.SFQED_init_PAIR_creation(str_path);
}

void SFQED_finalize_BREIT_WHEELER(){
        procs.SFQED_finalize_PAIR_creation();
}

/*********************************/
/* QUANTUM PARAMETER CALCULATION */
/*********************************/
double compute_chi_with_vectors(const double& gamma, const double p_in[3], const double EE[3], const double BB[3]){
        return procs.compute_quantum_param(gamma, p_in, EE, BB);
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

        return procs.compute_quantum_param(gamma, p_in_x, p_in_y, p_in_z, EE_x, EE_y, EE_z, BB_x, BB_y, BB_z);                                                                            
}


/*****************/
/* BUILD MOMENTA */
/*****************/
void SFQED_build_collinear_momentum(const double &gamma_out, const double p_in[3], double p_out[3]){
        procs.SFQED_collinear_momentum(gamma_out, p_in, p_out);
}

/********************************/
/* PHOTON EMISSION SECTION LCFA */
/********************************/
double SFQED_INV_COMPTON_rate(const double &gamma, const double &chi){
        return procs.SFQED_PHTN_emission_rate(gamma, chi);
}

double SFQED_LCFA_INV_COMPTON_PHOTON_energy(const double &gamma, const double &chi, const double &rnd){
        return procs.SFQED_LCFA_emitted_photon_energy(gamma, chi, rnd);
}

/*********************************/
/* PHOTON EMISSION SECTION BLCFA */
/*********************************/

BLCFA_Object* SFQED_CREATE_BLCFA_OBJECT(){
        return new BLCFA_Object();
}

void SFQED_FINALIZE_BLCFA_OBJECT(BLCFA_Object* entity){
        delete entity;
}

bool SFQED_BLCFA_OBJECT_update(BLCFA_Object* entity,
                                                const double* const pushed_momentum,
                                                const double* const momentum,
                                                double* Lorentz_F_t_der,
                                                double* Lorentz_F_tt_der,
                                                double& part_gamma,
                                                double& part_chi){
                                                        
        return entity->SFQED_BLCFA_update_entities_quantities(procs, pushed_momentum, momentum,
                                                Lorentz_F_t_der, Lorentz_F_tt_der, part_gamma, part_chi);
}

double SFQED_BLCFA_INV_COMPTON_PHOTON_threshold(const BLCFA_Object* const entity,
                                                        const double* const Lorentz_F_t_der,
                                                        const double* const Lorentz_F_tt_der,
                                                        const double& part_gamma,
                                                        const double& part_chi){

        return entity->SFQED_BLCFA_find_energy_threshold(procs, Lorentz_F_t_der, Lorentz_F_tt_der, part_gamma, part_chi);
}

double SFQED_BLCFA_INV_COMPTON_PHOTON_energy(const double& LCFA_limit,
                                                        const double& gamma, 
                                                        const double& chi, 
                                                        const double& rnd, 
                                                        const double& rnd2){

        return procs.SFQED_BLCFA_emitted_photon_energy(LCFA_limit, gamma, chi, rnd, rnd2);
}


//debug
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

/***************************/
/* PAIR PRODUCTION SECTION */
/***************************/
double SFQED_BREIT_WHEELER_rate(const double &gamma, const double &chi){
        return procs.SFQED_PAIR_creation_rate(gamma, chi);
}

double SFQED_BREIT_WHEELER_rate_fast(const double &gamma, const double &chi){
        return procs.SFQED_PAIR_creation_rate_fast(gamma, chi);
}

double SFQED_BREIT_WHEELER_ELECTRON_energy(const double &nrg, const double &chi, const double &rnd){
        return procs.SFQED_PAIR_created_electron_energy(nrg, chi, rnd);
}

double SFQED_BREIT_WHEELER_ELECTRON_energy_fast(const double &nrg, const double &chi, const double &rnd){
        return procs.SFQED_PAIR_created_electron_energy_fast(nrg, chi, rnd);
}
