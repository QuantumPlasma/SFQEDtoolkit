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

#include <cmath>

//for the moment we keep the SFQED_Processes and the BLCFA_Object
// classes intertwined by forcing the user to compile both (since 
// the 2 classes have friend methods). In the future we would like 
// to have a new class, derived from SFQED_Processes, implementing
// the BLCFA mechanisms. This would require a change of the fortran
// and c++ interfaces.

struct SFQED_Processes{
private:
        /*********************/
        /* GENERAL VARIABLES */
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

        double norm_Compton_time;

public:          

        SFQED_Processes(){
                //do nothing since you have to explicitly
                // call for the setters of the simulation
                // parameters
        }


        /*******************************/
        /* COEFFICIENTS GETTER SECTION */
        /*******************************/

        double get_phtn_rate_coefficient(){
                return coef_rate_Wrad;
        }

        double get_pair_rate_coefficient(){
                return coef_rate_Wpair;
        }

        double get_normalized_compton_time(){
                return norm_Compton_time;
        }


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

};

/***************************/
/* GENERIC ROUTINE SECTION */
/***************************/
void SFQED_collinear_momentum(const double&, const double[3], double[3]);

#endif