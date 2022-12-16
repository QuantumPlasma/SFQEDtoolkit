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
#ifndef PAIRS_TOOLS
#define PAIRS_TOOLS

#include "Chebyshev_User_1D.h"
#include "Chebyshev_User_2D.h"
#include "SFQED_User_Constants.h"

#include <cmath>

//you should make this a class!!!

class Breit_Wheeler_Pair_Production {
private:
        /***********************/
        /* PAIRS lookup tables */
        /***********************/
        //pair production rate
        Chebyshev_User_1D look_up_table_pair_prd_rt[7];
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

        //auxiliary function used in the energy calculation
        double SFQED_PAIR_emitted_electron_energy_aux(const int &lookup_index, const double &chi, const double &rescaled_rnd) const;

public:
        /********************************************/
        /* coefficients initializers and finalizers */
        /********************************************/
        //pairs
        void SFQED_init_PAIR_creation(std::string);

        void SFQED_finalize_PAIR_creation();

        /***************************/
        /* PAIR PRODUCTION SECTION */
        /***************************/

        inline double __attribute__((always_inline)) SFQED_PAIR_creation_rate(const double &gamma, const double &chi) const {

                int lookup_index;

                if(chi <= bound_kappa_0th_rate){
                        // lookup_index = 6;
                        return (9. * pigreek * chi)/(16. * sqrt(2.)) * exp(-(8/(3 * chi))) * (1. - 11./64.*chi + 7585./73728.*chi*chi) / gamma;
                }else if(chi <= bound_kappa_intermediate_rate){
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

                return (look_up_table_pair_prd_rt[lookup_index]).evaluate(chi) / gamma;//
        }

        inline double __attribute__((always_inline)) SFQED_PAIR_creation_rate_fast(const double &gamma, const double &chi) const {
                        
                int lookup_index;

                //we use bound_kappa_0th_nrgs in such a way that
                // the emission under bound_kappa_0th_nrgs = 0.3 
                // is still skipped!
                if(chi <= bound_kappa_0th_nrgs){
                        lookup_index = 6;
                }else if(chi <= bound_kappa_intermediate_rate){
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

                        
                return (look_up_table_pair_prd_rt[lookup_index]).evaluate(chi) / gamma;//
        }

        double SFQED_PAIR_created_electron_energy(const double&, const double&, const double&) const;

        double SFQED_PAIR_created_electron_energy_fast(const double&, const double&, const double&) const;

};

#endif