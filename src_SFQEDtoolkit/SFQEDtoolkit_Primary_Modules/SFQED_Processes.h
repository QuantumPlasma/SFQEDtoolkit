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

#pragma once
#ifndef SFQED_PROCS
#define SFQED_PROCS

#include "Chebyshev_User_1D.h"
#include "Chebyshev_User_2D.h"
#include "BLCFA_Object.h"

#include <cmath>

// #include "SFQED_Functions.h"

//upper bounds of the \chi parameter in each interval
const double bound_chi_1st = 2.0,
        bound_chi_2nd = 20.0,
        bound_chi_3rd = 80.0,
        bound_chi_4th = 600.0,
        bound_chi_5th = 2000.0;

const double
        bound_kappa_0th = 0.3,
        bound_kappa_1st = 2.0,
        bound_kappa_2nd = 20.0,
        bound_kappa_3rd = 80.0,
        bound_kappa_4th = 600.0,
        bound_kappa_5th = 2000.0;

const double pigreek = 3.141592653589793238462643383279503;

// //the rates lookup tables are defined in the cpp file
// extern Chebyshev_User_1D look_up_table_pair_prd_rt[6], look_up_table_phtn_mx_rt[5];


struct SFQED_Processes{
private:
        /*********************/
        /* GENEAL VARIABLES */
        /*********************/
        //reference length and angular frequency
        double Lambda, omega_r;

        // static double zbrent(double (*)(double,const double*), const double&, const double&, const double&, const double* const);

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


        /*********************************/
        /* LCFA subsection lookup tables */
        /*********************************/
        //photon emission rate
        Chebyshev_User_1D look_up_table_phtn_mx_rt[5];
        //photon partial emission rate (used in the exponential tale asymptotic solution)
        Chebyshev_User_2D look_up_table_phtn_prtl_rt[5];
        //photon energy in w values
        Chebyshev_User_2D look_up_table_phtn_momw[5];
        //photon energy in 1/w values
        Chebyshev_User_2D look_up_table_1_over_w[5];
        //1D coefficient projection of the previous 2D coefficients
        //(used in the exponential tale asymptotic solution)
        Chebyshev_User_1D look_up_table_phtn_1_o_w_proj[5];

        /**********************************/
        /* BLCFA subsection lookup tables */
        /**********************************/
        //photon differential cross section (used in the BLCFA mechanism)
        Chebyshev_User_2D look_up_table_phtn_diff_crss_sctn[5];

        /**********************************/
        /* PAIRS subsection lookup tables */
        /**********************************/
        //pair production rate
        Chebyshev_User_1D look_up_table_pair_prd_rt[6];
        //pair's electron energy in v values
        //(coefficients of the inverse function W[k,v] - rW[k,0] = 0.)
        Chebyshev_User_2D look_up_table_pair_v_nrgs[6];
        //pair's electron energy in v values for the intermediate region of v
        Chebyshev_User_2D look_up_table_pair_v_nrgs_high[6];
        //1D projections of the above brent inverse 2D coeffs upon the corresponding r high limit value
        //(used in the exponential tale asymptotic solution)
        Chebyshev_User_1D look_up_table_pair_v_nrgs_proj[6];
        //pair production partial rate 
        //(used in the exponential tale asymptotic solution)
        Chebyshev_User_2D look_up_table_pair_prtl_rt[6];


        /***************************/
        /* PHOTON EMISSION SECTION */
        /***************************/    

        //auxiliary inner function
        //it computes the w-value energy of the emitted photon
        double SFQED_LCFA_phtn_nrg_aux(const double&, const double&, const int&) const;
        
        /***************************/
        /* PAIR PRODUCTION SECTION */
        /***************************/

        //pair production auxiliary function
        double SFQED_PAIR_emitted_electron_energy_aux(const int &, const double &, const double &) const;

public:

        SFQED_Processes(){
                //do nothing
        }

        /********************************************/
        /* coefficients initializers and finalizers */
        /********************************************/
        //photon
        void SFQED_init_PHTN_emission(std::string);

        void SFQED_finalize_PHTN_emission();

        //pairs
        void SFQED_init_PAIR_creation(std::string);

        void SFQED_finalize_PAIR_creation();

        /**********************************/
        /* SIMULATION INITIALIZER SECTION */
        /**********************************/
        //set all the simulation's quantities based on
        //the value of the given reference length
        void SFQED_set_reference_length(const double&);

        void SFQED_set_reference_angular_frequency(const double&);

        void SFQED_set_all_to_one();

        //needed in the BLCFA calculation 
        void SFQED_set_time_step(const double&);

        /************************/
        /* QUANTUM CHI ROUTINES */
        /************************/

        // this function is used for the calculation of \chi (for electron)
        // or \kappa (photon).
        // this form provides the most accurate result with also very good performance

       
        //MIND THE FIELDS NORMALIZED UNITS
        //remember that, if the momentum of the fermion p_in is given in terms of "m_e c", than
        //gamma = sqrt(1 + |p_in|^2)
        // this function should be used with p_in corresponding to the e- or e+ momentum at half-step
        // calculate the \chi parameter. Momentum and fields MUST BE in units of
        // "m_e c" and of "m_e c \omega / |e|" for EE and "m_e \omega / |e|" for BB, respectively
        inline double __attribute__((always_inline)) compute_quantum_param(const double& gamma,
                                                                                const double p_in[3],
                                                                                const double EE[3],
                                                                                const double BB[3]) const
        {
                double cross[3], ff[3];
    
                // the calculation of \chi or \kappa in this form provides
                // the most accurate result with also very good performance
                cross[0] = p_in[1]*BB[2] - p_in[2]*BB[1];
                cross[1] = p_in[2]*BB[0] - p_in[0]*BB[2];
                cross[2] = p_in[0]*BB[1] - p_in[1]*BB[0];
                
                ff[0] = gamma * EE[0] + cross[0];
                ff[1] = gamma * EE[1] + cross[1];
                ff[2] = gamma * EE[2] + cross[2];
                
                double arg = p_in[0]*EE[0] + p_in[1]*EE[1] + p_in[2]*EE[2];
                arg = ff[0]*ff[0] + ff[1]*ff[1] + ff[2]*ff[2] - (arg*arg);
                
                //if ff is parallel to p_in and large, then round-off errors
                //(or an invalid input gam factor) may render arg negative
                //=> NaN is returned instead of zero
                //we avoid this by checking the sign of arg
                arg = arg < 0. ? 0. : arg;
                
                //pay attention to the UNITS
                //return twopiComptonDivLambda * sqrt(arg);
                return ComptonDivLambda * sqrt(arg);
        }

        inline double __attribute__((always_inline)) compute_quantum_param(const double &gamma,
                                                                                const double &p_in_x,
                                                                                const double &p_in_y,
                                                                                const double &p_in_z,
                                                                                const double &EE_x,
                                                                                const double &EE_y,
                                                                                const double &EE_z,
                                                                                const double &BB_x,
                                                                                const double &BB_y,
                                                                                const double &BB_z) const
        {
                double cross[3], ff[3];
    
                // the calculation of \chi or \kappa in this form provides
                // the most accurate result with also very good performance
                cross[0] = p_in_y*BB_z - p_in_z*BB_y;
                cross[1] = p_in_z*BB_x - p_in_x*BB_z;
                cross[2] = p_in_x*BB_y - p_in_y*BB_x;
                
                ff[0] = gamma * EE_x + cross[0];
                ff[1] = gamma * EE_y + cross[1];
                ff[2] = gamma * EE_z + cross[2];
                
                double arg = p_in_x * EE_x + p_in_y * EE_y + p_in_z * EE_z;
                arg = ff[0]*ff[0] + ff[1]*ff[1] + ff[2]*ff[2] - (arg*arg);
                
                //if ff is parallel to p_in and large, then round-off errors
                //(or an invalid input gam factor) may render arg negative
                //=> NaN is returned instead of zero
                //we avoid this by checking the sign of arg
                arg = arg < 0. ? 0. : arg;
                
                //pay attention to the UNITS
                //return twopiComptonDivLambda * sqrt(arg);
                return ComptonDivLambda * sqrt(arg);
        }


        /***************************/
        /* PHOTON EMISSION SECTION */
        /***************************/
        
        /*******************/
        /* LCFA subsection */
        /*******************/

        //the function below returns the rate of photon emission for an electron
        //having  a certain chi and gamma (gamma is the relativistic factor of 
        //the electron, but it can also be interpreted as its energy normalized in 
        //units of m_e c^2)
        inline double __attribute__((always_inline)) SFQED_PHTN_emission_rate(const double &gamma, const double &chi) const {
        
                //coefficient for the rate of emission
                //multiply by \tilde{W_rad} to get the rate of emission
                double coefW_rad = coef_rate_Wrad * chi / gamma;


                int lookup_index;

                if(chi <= bound_chi_1st){
                        lookup_index = 4;
                }else if(chi <= bound_chi_2nd){
                        lookup_index = 3;
                }else if(chi <= bound_chi_3rd){
                        lookup_index = 2;
                }else if(chi <= bound_chi_4th){
                        lookup_index = 1;
                }else{
                        lookup_index = 0;
                }

                
                return coefW_rad * (look_up_table_phtn_mx_rt[lookup_index]).evaluate(chi);
        }

        double SFQED_LCFA_emitted_photon_energy(const double&, const double&, const double&) const;

        //this function will be moved outside in the next versions
        void SFQED_collinear_momentum(const double&, const double[3], double[3]);

        /********************/
        /* BLCFA subsection */
        /********************/

        //define friend functions
        friend bool BLCFA_Object::SFQED_BLCFA_update_entities_quantities(const SFQED_Processes&, const double* const, const double* const, double*, double*, double&, double&);

        friend double BLCFA_Object::SFQED_BLCFA_find_energy_threshold(const SFQED_Processes&, const double* const, const double* const, const double&, const double&) const;

        double SFQED_BLCFA_emitted_photon_energy(const double&, const double&, const double&, const double&, const double&) const;


        /***************************/
        /* PAIR PRODUCTION SECTION */
        /***************************/

        inline double __attribute__((always_inline)) SFQED_PAIR_creation_rate(const double &gamma, const double &chi) const {
    
                //coefficient for the rate of emission
                //multiply by \tilde{W_rad} to get the rate of emission
                double coefW_pair = coef_rate_Wpair / gamma;


                int lookup_index;

                if(chi <= bound_kappa_0th){
                        // lookup_index = 5;
                        return coefW_pair * //asympt_pair_emission_rate_low_k(chi);
                                        (27. * pigreek * chi)/(16. * sqrt(2.)) * exp(-(8/(3 * chi))) * (1. - 11./64.*chi + 7585./73728.*chi*chi);
                }else if(chi <= bound_kappa_1st){
                        lookup_index = 4;
                }else if(chi <= bound_kappa_2nd){
                        lookup_index = 3;
                }else if(chi <= bound_kappa_3rd){
                        lookup_index = 2;
                }else if(chi <= bound_kappa_4th){
                        lookup_index = 1;
                }else{
                        lookup_index = 0;
                }                   
                
                return coefW_pair * (look_up_table_pair_prd_rt[lookup_index]).evaluate(chi);//
        }

        inline double __attribute__((always_inline)) SFQED_PAIR_creation_rate_fast(const double &gamma, const double &chi) const {
    
                //coefficient for the rate of emission
                //multiply by \tilde{W_rad} to get the rate of emission
                double coefW_pair = coef_rate_Wpair / gamma;
                
                
                int lookup_index;

                if(chi <= bound_kappa_0th){
                        lookup_index = 5;
                }else if(chi <= bound_kappa_1st){
                        lookup_index = 4;
                }else if(chi <= bound_kappa_2nd){
                        lookup_index = 3;
                }else if(chi <= bound_kappa_3rd){
                        lookup_index = 2;
                }else if(chi <= bound_kappa_4th){
                        lookup_index = 1;
                }else{
                        lookup_index = 0;
                }

                
                return coefW_pair * (look_up_table_pair_prd_rt[lookup_index]).evaluate(chi);//
        }

        double SFQED_PAIR_created_electron_energy(const double&, const double&, const double&) const;

        double SFQED_PAIR_created_electron_energy_fast(const double&, const double&, const double&) const;

       
};

#endif