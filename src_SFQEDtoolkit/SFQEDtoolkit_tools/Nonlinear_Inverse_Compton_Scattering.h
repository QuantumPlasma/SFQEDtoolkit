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
#ifndef LCFA_and_PAIRS_TOOLS
#define LCFA_and_PAIRS_TOOLS

#include "Chebyshev_User_1D.h"
#include "Chebyshev_User_2D.h"
#include "SFQED_User_Constants.h"

class Nonlinear_Inverse_Compton_Scattering{
private:
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

    //auxiliary function used in the energy calculation
    double SFQED_LCFA_phtn_nrg_aux(const double &chi, const double &rnd, const int& lookup_index) const;

    double SFQED_BLCFA_phtn_nrg_aux(const double &chi, const double &rnd, const int& lookup_index) const;

public:
    /********************************************/
    /* coefficients initializers and finalizers */
    /********************************************/
    //photon LCFA
    void SFQED_init_PHTN_emission(std::string);

    void SFQED_finalize_PHTN_emission();


    //photon BLCFA
    void SFQED_init_PHTN_emission_BLCFA(std::string);

    void SFQED_finalize_PHTN_emission_BLCFA();

    /***************************/
    /* PHOTON EMISSION SECTION */
    /***************************/

    /********/
    /* LCFA */
    /********/

    //the function below returns the rate of photon emission for an electron
    //having a certain chi and gamma (gamma is the relativistic factor of 
    //the electron, but it can also be interpreted as its energy normalized in 
    //units of m_e c^2). However the quantity returned is not expressed in the
    // proper code units. and should be multiplied by the coef_rate_Wrad attribute,
    // which specific of each simulation, and gets initilaized in the SFQED_Processes
    // instance
    inline double __attribute__((always_inline)) SFQED_PHTN_emission_rate(const double &gamma, const double &chi) {

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

            return chi / gamma * (look_up_table_phtn_mx_rt[lookup_index]).evaluate(chi);
    }

    double SFQED_LCFA_emitted_photon_energy(const double &gamma, const double &chi, const double &rnd) const;


    /*********/
    /* BLCFA */
    /*********/
    double SFQED_BLCFA_emitted_photon_energy(const double& LCFA_limit,
                                                                const double& gamma, 
                                                                const double& chi, 
                                                                const double& rnd, 
                                                                const double& rnd2) const;

    double SFQED_BLCFA_emitted_photon_energy_2(const double& LCFA_limit,
                                                            const double& gamma_photon,
                                                            const double& gamma, 
                                                            const double& chi, 
                                                            const double& rnd2) const;

};

#endif