#pragma once
#ifndef SFQEDtoolkit_INTERFACE
#define SFQEDtoolkit_INTERFACE

class BLCFA_Object;

extern "C" {

        /**************************/
        /* SFQED process instance */
        /**************************/
        void SFQED_INIT_PROCS();


        /**********************************/
        /* SIMULATION INITIALIZER SECTION */
        /**********************************/
        /****  COMPLETE SETTERS AND FINALIZERS  ****/
        bool SFQED_INIT_ALL_ref_len(const double& ref_len, const double& ts);

        bool SFQED_INIT_ALL_ref_freq(const double& ref_freq, const double& ts);

        void SFQED_FINALIZE_ALL();
        
        /****  PARTIAL SETTERS  ****/

        //set all the simulation's quantities based on
        //the value of the given reference length
        void SFQED_set_ref_len(const double& ref_len);

        void SFQED_set_ref_freq(const double& ref_freq);

        void SFQED_set_sim_tstep(const double& ts);

        //useful for debug
        void SFQED_set_for_debug();


        /********************************************/
        /* coefficients initializers and finalizers */
        /********************************************/
        /* PHOTON */
        void SFQED_init_INV_COMPTON(const char* const path);

        void SFQED_finalize_INV_COMPTON();

        /* PAIRS */
        void SFQED_init_BREIT_WHEELER(const char* const path);

        void SFQED_finalize_BREIT_WHEELER();


        /*********************************/
        /* QUANTUM PARAMETER CALCULATION */
        /*********************************/
        double compute_chi_with_vectors(const double& gamma, const double p_in[3], const double EE[3], const double BB[3]);

        double compute_chi_with_components(const double &gamma,
                                                const double &p_in_x, const double &p_in_y, const double &p_in_z,
                                                const double &EE_x, const double &EE_y, const double &EE_z,
                                                const double &BB_x, const double &BB_y, const double &BB_z);


        /*****************/
        /* BUILD MOMENTA */
        /*****************/
        void SFQED_build_collinear_momentum(const double &gamma_out, const double p_in[3], double p_out[3]);


        /********************************/
        /* PHOTON EMISSION SECTION LCFA */
        /********************************/
        double SFQED_INV_COMPTON_rate(const double &gamma, const double &chi);

        double SFQED_LCFA_INV_COMPTON_PHOTON_energy(const double &gamma, const double &chi, const double &rnd);


        /*********************************/
        /* PHOTON EMISSION SECTION BLCFA */
        /*********************************/
        BLCFA_Object* SFQED_CREATE_BLCFA_OBJECT();

        void SFQED_FINALIZE_BLCFA_OBJECT(BLCFA_Object* entity);

        bool SFQED_BLCFA_OBJECT_update(BLCFA_Object* entity,
                                                        const double* const pushed_momentum,
                                                        const double* const momentum,
                                                        double* Lorentz_F_t_der,
                                                        double* Lorentz_F_tt_der,
                                                        double& part_gamma,
                                                        double& part_chi);

        double SFQED_BLCFA_INV_COMPTON_PHOTON_threshold(const BLCFA_Object* const entity,
                                                                const double* const Lorentz_F_t_der,
                                                                const double* const Lorentz_F_tt_der,
                                                                const double& part_gamma,
                                                                const double& part_chi);

        double SFQED_BLCFA_INV_COMPTON_PHOTON_energy(const double& LCFA_limit,
                                                                const double& gamma, 
                                                                const double& chi, 
                                                                const double& rnd, 
                                                                const double& rnd2);

        //debug
        void SFQED_BLCFA_printALL(BLCFA_Object* entity);

        void SFQED_BLCFA_get_F_old(BLCFA_Object* entity, double* F_old);

        void SFQED_BLCFA_update_F_old(BLCFA_Object* entity, const double* const F_old_new);

        void SFQED_BLCFA_get_delta_F_old(BLCFA_Object* entity, double* delta_F_old);


        /***************************/
        /* PAIR PRODUCTION SECTION */
        /***************************/
        double SFQED_BREIT_WHEELER_rate(const double &gamma, const double &chi);

        double SFQED_BREIT_WHEELER_rate_fast(const double &gamma, const double &chi);

        double SFQED_BREIT_WHEELER_ELECTRON_energy(const double &nrg, const double &chi, const double &rnd);

        double SFQED_BREIT_WHEELER_ELECTRON_energy_fast(const double &nrg, const double &chi, const double &rnd);

}

#endif