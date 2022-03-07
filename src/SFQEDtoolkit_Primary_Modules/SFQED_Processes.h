#pragma once
#ifndef SFQED_PROCS
#define SFQED_PROCS

#include "Chebyshev_User_1D_PH.h"
//#include "Chebyshev_User_1D.h"
#include "Chebyshev_User_2D_PH.h"
// #include "Chebyshev_User_2D.h"
#include "Chebyshev_User_2D_PH.h"
#include "BLCFA_Object.h"

// #include "Chebyshev_Coefficients.h"
// #include "SFQED_Functions.h"

struct SFQED_Processes{
private:
        /*********************/
        /* OVERALL VARIABLES */
        /*********************/
        //reference length and angular frequency
        double Lambda, omega_r;

        //simple constants initialized during
        //the reference length initialization
        double twopiComptonDivLambda,
                ComptonDivLambda,
                LambdaDivtwopiCompton,
                LambdaDivCompton,
                coef_rate_Wrad,
                coef_rate_Wpair;
        double norm_Compton_time,
                norm_Compton_time_2;

        double one_over_dt,
                one_over_dt2;

        /***************************/
        /* PHOTON EMISSION SECTION */
        /***************************/

        /*******************/
        /* LCFA subsection */
        /*******************/
        Chebyshev_User_1D phtn_mx_rt_0_2,
                        phtn_mx_rt_2_20,
                        phtn_mx_rt_20_80,
                        phtn_mx_rt_80_600,
                        phtn_mx_rt_600_2000;

        //this look up table will be employed to mitigate the
        //aftermath of a branch misprediction
        //the size is due to (0+1)+1+2+4+8 = 16 = 2^(5-1)
        //notice that objects with module scope should be static by default
        Chebyshev_User_1D *look_up_table_phtn_mx_rt[16];

        Chebyshev_User_2D phtn_prtl_rt_0_2,
                        phtn_prtl_rt_2_20,
                        phtn_prtl_rt_20_80,
                        phtn_prtl_rt_80_600,
                        phtn_prtl_rt_600_2000;

        Chebyshev_User_2D *look_up_table_phtn_prtl_rt[16];

        Chebyshev_User_2D phtn_momw_0_2,
                        phtn_momw_2_20,
                        phtn_momw_20_80,
                        phtn_momw_80_600,
                        phtn_momw_600_2000;

        Chebyshev_User_2D *look_up_table_phtn_momw[16];

        Chebyshev_User_2D phtn_1_over_w_0_2,
                        phtn_1_over_w_2_20,
                        phtn_1_over_w_20_80,
                        phtn_1_over_w_80_600,
                        phtn_1_over_w_600_2000;

        Chebyshev_User_2D *look_up_table_1_over_w[16];

        Chebyshev_User_1D phtn_1_over_w_proj_0_2,
                        phtn_1_over_w_proj_2_20,
                        phtn_1_over_w_proj_20_80,
                        phtn_1_over_w_proj_80_600,
                        phtn_1_over_w_proj_600_2000;

        Chebyshev_User_1D *look_up_table_phtn_1_o_w_proj[16];

        //auxiliary inner function
        //it computes the w-value energy of the emitted photon
        double SFQED_LCFA_phtn_nrg_aux(const double&, const double&, const int&);

        /********************/
        /* BLCFA subsection */
        /********************/

        Chebyshev_User_2D phtn_diff_crss_sctn_0_2,
                        phtn_diff_crss_sctn_2_20,
                        phtn_diff_crss_sctn_20_80,
                        phtn_diff_crss_sctn_80_600,
                        phtn_diff_crss_sctn_600_2000;

        Chebyshev_User_2D *look_up_table_phtn_diff_crss_sctn[16];

        /***************************/
        /* PAIR PRODUCTION SECTION */
        /***************************/
        
        Chebyshev_User_1D_PH pair_prd_rt_001_03, pair_prd_rt_001_03_fast;

        //these will host the rate of pair production
        Chebyshev_User_1D pair_prd_rt_0_2,
                        pair_prd_rt_2_20,
                        pair_prd_rt_20_80,
                        pair_prd_rt_80_600,
                        pair_prd_rt_600_2000;

        Chebyshev_User_1D *look_up_table_pair_prd_rt[32];

        Chebyshev_User_1D *look_up_table_pair_prd_rt_fast[32];

        //these will store the coefficients of the inverse function W[k,v] - rW[k,0] = 0.

        Chebyshev_User_2D pair_v_nrgs_001_03,
                        pair_v_nrgs_0_2,
                        pair_v_nrgs_2_20,
                        pair_v_nrgs_20_80,
                        pair_v_nrgs_80_600,
                        pair_v_nrgs_600_2000;

        Chebyshev_User_2D *look_up_table_pair_v_nrgs[32];

        Chebyshev_User_2D pair_v_nrgs_001_03_high,
                        pair_v_nrgs_0_2_high,
                        pair_v_nrgs_2_20_high,
                        pair_v_nrgs_20_80_high,
                        pair_v_nrgs_80_600_high,
                        pair_v_nrgs_600_2000_high;

        Chebyshev_User_2D *look_up_table_pair_v_nrgs_high[32];

        //projections of the above brent inverse coeffs
        //upon the corresponding r low limit value
        Chebyshev_User_1D pair_v_nrg_proj_001_03,
                        pair_v_nrg_proj_0_2,
                        pair_v_nrg_proj_2_20,
                        pair_v_nrg_proj_20_80,
                        pair_v_nrg_proj_80_600,
                        pair_v_nrg_proj_600_2000;

        Chebyshev_User_1D *look_up_table_pair_v_nrgs_proj[32];

        //these will approximate the differential probability integrated
        //from a certain v to 1. Be careful, the v domain is not complete!
        //it does not need to be

        Chebyshev_User_2D_PH pair_prtl_rt_001_03;

        Chebyshev_User_2D pair_prtl_rt_0_2,
                        pair_prtl_rt_2_20,
                        pair_prtl_rt_20_80,
                        pair_prtl_rt_80_600,
                        pair_prtl_rt_600_2000;

        Chebyshev_User_2D *look_up_table_pair_prtl_rt[32];

        //pair production auxiliary function
        double SFQED_PAIR_emitted_electron_energy_aux(const int &, const double &, const double &);

public:

        //debug purposes
        // double norm_Compton_time,
        //         norm_Compton_time_2;

        SFQED_Processes(){
                //do nothing
        }

        //overall class init routines
        void SFQED_set_reference_length(const double&);

        void SFQED_set_reference_angular_frequency(const double&);

        void SFQED_set_time_step(const double&);

        // COMMONLY USED ROUTINES 

        // this function is used for the calculation of \chi (for electron)
        // or \kappa (photon).
        // this form provides the most accurate result with also very good performance
        double SFQED_compute_quantum_param(double, const double[3], const double[3], const double[3]);

        double SFQED_compute_quantum_param(const double &gamma,
                                                const double &p_in_x,
                                                const double &p_in_y,
                                                const double &p_in_z,
                                                const double &EE_x,
                                                const double &EE_y,
                                                const double &EE_z,
                                                const double &BB_x,
                                                const double &BB_y,
                                                const double &BB_z);

        /***************************/
        /* PHOTON EMISSION SECTION */
        /***************************/

        //init and finalize routines
        void SFQED_init_PHTN_emission(std::string);

        void SFQED_finalize_PHTN_emission();

        //debug no longer needed
        // double* intdWrad_0chi2;
        // void initialize_intdWrad(std::string);

        /*******************/
        /* LCFA subsection */
        /*******************/

        double SFQED_PHTN_emission_rate(const double&, const double&);

        double SFQED_LCFA_emitted_photon_energy(const double&, const double&, const double&);

        void SFQED_collinear_momentum(const double&, const double[3], double[3]);

        /********************/
        /* BLCFA subsection */
        /********************/

        //define friend functions
        friend bool BLCFA_Object::SFQED_BLCFA_update_entities_quantities(const SFQED_Processes&, const double* const, const double* const, double*, double*, double&, double&);

        friend double BLCFA_Object::SFQED_BLCFA_find_energy_threshold(const SFQED_Processes&, const double* const, const double* const, const double&, const double&);

        double SFQED_BLCFA_emitted_photon_energy(const double&, const double&, const double&, const double&, const double&);

        //these are no longer used
        // double SFQED_BLCFA_emitted_photon_energy_derivative(const double&, const double&, const double&, const double&, const double&);

        // double SFQED_BLCFA_emitted_photon_energy_no_approx(const double&, const double&, const double&, const double&, const double&);

        // double SFQED_BLCFA_matteo(const double&, const double&, const double&, const double&, const double&);

        // double SFQED_BLCFA_DEBUG(const double&, const double&, const double&, const double&, const double&, std::ofstream&);


        /***************************/
        /* PAIR PRODUCTION SECTION */
        /***************************/

        void SFQED_init_PAIR_creation(std::string);

        double SFQED_PAIR_creation_rate(const double&, const double&);

        double SFQED_PAIR_creation_rate_fast(const double&, const double&);

        double SFQED_PAIR_created_electron_energy(const double&, const double&, const double&);

        double SFQED_PAIR_created_electron_energy_fast(const double&, const double&, const double&);

        void SFQED_finalize_PAIR_creation();
       
};

#endif