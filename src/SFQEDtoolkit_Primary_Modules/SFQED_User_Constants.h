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

#ifndef USER_CONST
#define USER_CONST

//useful constants used throughout the module
const double pi = 3.141592653589793238462643383279503;
const double twopi = 6.283185307179586476925286766559006;
const double four_over_threepi = 4. / (3. * pi);
const double sqrt3 = 1.732050807568877293527446341505872;

//light speed [m/s]
const double light_speed = 299792458.;

///////////////////////////////////////////////////////////////
//physical constants relevant to SFQED processes

//fine structure constant
const double alpha = 7.2973525698e-3;

//Compton length [ \hbar / m_e c ], in meters
const double compton_length = 3.86159268e-13;

//Compton time [(\hbar / m_e c) / c ], in seconds
const double compton_time = 1.2880887e-21;

//squared Compton time 
// const double compton_time_2 = compton_time * compton_time;

//2*\pi*Compton length [ h / m_e c ], in meters
const double twopiCompton = twopi*compton_length;
  
//coef. for the rate of photon emission W_rad
const double coef_W_rad = alpha/(sqrt3*pi);

//coef. for the rate of pair creation W_pair
const double coef_W_pair = alpha/(3.0*sqrt3*pi);

//constant parameter needed in the BLCFA mechanism to
//ascertain whether the emission occurs or not
const double zeta = 2.22e-16;
const double zeta_2 = zeta*zeta;

const double thresh_factor = 0.7;

/************************************/
/*PHOTON EMISSION SPECIFIC CONSTANTS*/
/************************************/

//upper bounds of the \chi parameter in each interval
const double bound_chi_1st = 2.0,
        bound_chi_2nd = 20.0,
        bound_chi_3rd = 80.0,
        bound_chi_4th = 600.0,
        bound_chi_5th = 2000.0;

const double heval_r_0_2 = 0.99, // 0.9999;
        heval_r_2_20 = 0.9999,
        heval_r_20_80 = 0.9999, //0.99959,
        heval_r_80_600 = 0.9998,
        heval_r_600_2000 = 0.99984;

//new lower bounds on r (much simpler and straightforward than the previous hybrid way)
double lu_table_lower_r_bounds[5];

double lu_table_inverse_r_bounds[5];

//the following represents the coefficient of the term in w
//when taking the small w limit (w -> 0) of the function
//g(\chi,w) = \int_0^w{W_rad(\chi, w_{integrated}) dw_{integrated}}:
//(using Series[integral[w],{w,0,2}]//Normal)
//\int_0^w{W_rad(\chi, w) dw} \approx [9 \Gamma(2/3) / 2^(1/3)] * w.
//at that point, since the equation to solve has the form
//g(\chi,w) = rG(\chi) --(w->0)--> [9 \Gamma(2/3) / 2^(1/3)] * w = r*G(\chi),
//with G(\chi) = \tilde{W_rad}
//we can easily invert it to find a value for w
//w \approx [r * G(\chi)] / [9 \Gamma(2/3) / 2^(1/3)]
double soft_ph_coef = 1.03381857436581047459e-1;

//new upper bounds on r (much simpler and straightforward than the previous hybrid way)
double lu_table_upper_r_bounds[5];

//this constants should represent the upper bound
//of the emitted phtn energy, expressed as w value,
//above which the differential probability approximation
//is no longer valid. When this happens the emission
//does not occur
double lu_table_upper_w_bounds[5];

//this function is meant to initialize every photon 
//emission lookup table defined above
void init_phtn_mx_tables(){

    //new r edges init
    lu_table_lower_r_bounds[4] = 0.028; //0.020; //w_min = 0.017
    lu_table_lower_r_bounds[3] = 0.041; //0.031; //w_min = 0.0158
    lu_table_lower_r_bounds[2] = 0.045; //0.035; //w_min = 0.01
    lu_table_lower_r_bounds[1] = 0.040; //0.029; //w_min = 0.0055
    lu_table_lower_r_bounds[0] = 0.050; //0.039; //w_min = 0.0038


    // "maybe these limits could be moved a little backwards?"
    // The answer is NO and, indeed, it is the exact opposite:
    // while the brent inverse (w) is regular basicaclly
    // throughout the whole domain and has problems at the very end,
    // its reciprocal (1/w), on the contrary, is regular at the ending
    // part of the domain and has issues everywhere else. If we tried to
    // extend the 1/w domain backwards, things would start going far worse. 

    lu_table_inverse_r_bounds[4] = 0.986; //0.980905;
    lu_table_inverse_r_bounds[3] = 0.987; //0.980905;
    lu_table_inverse_r_bounds[2] = 0.987; //0.980605;
    lu_table_inverse_r_bounds[1] = 0.987; //0.980705,
    lu_table_inverse_r_bounds[0] = 0.987; //0.980805;

    lu_table_upper_r_bounds[4] = 0.99999;
    lu_table_upper_r_bounds[3] = 0.99999;
    lu_table_upper_r_bounds[2] = 0.99996;
    lu_table_upper_r_bounds[1] = 0.99997;
    lu_table_upper_r_bounds[0] = 0.99998;

    lu_table_upper_w_bounds[4] = 1.6, //2.2999,
    lu_table_upper_w_bounds[3] = 1., //2.,
    lu_table_upper_w_bounds[2] = 0.47, //1.8,
    lu_table_upper_w_bounds[1] = 0.3, //1.7,
    lu_table_upper_w_bounds[0] = 0.15; //1.5;

}

/************************************/
/*PAIR PRODUCTION SPECIFIC CONSTANTS*/
/************************************/

const double
        bound_kappa_0th = 0.3,
        bound_kappa_1st = 2.0,
        bound_kappa_2nd = 20.0,
        bound_kappa_3rd = 80.0,
        bound_kappa_4th = 600.0,
        bound_kappa_5th = 2000.0;

//change this shit
const double pair_r_exp_limit_001_03 = 0.9999,
        pair_r_exp_limit_0_2 = 0.9999,
        pair_r_exp_limit_2_20 = 0.9999,
        pair_r_exp_limit_20_80 = 0.9999,
        pair_r_exp_limit_80_600 = 0.9999,
        pair_r_exp_limit_600_2000 = 0.9999;

double lu_table_pair_lower_r_bounds[6];

double lu_table_pair_upper_r_bounds[6];

void init_pair_crtn_tables(){

    //new r edges init
    lu_table_pair_lower_r_bounds[5] = 0.899;
    lu_table_pair_lower_r_bounds[4] = 0.899; // 0.005;
    lu_table_pair_lower_r_bounds[3] = 0.899; // 0.001;
    lu_table_pair_lower_r_bounds[2] = 0.899; // 0.003;
    lu_table_pair_lower_r_bounds[1] = 0.899; 
    lu_table_pair_lower_r_bounds[0] = 0.899;

    lu_table_pair_upper_r_bounds[5] = 0.9999;
    lu_table_pair_upper_r_bounds[4] = 0.9999;
    lu_table_pair_upper_r_bounds[3] = 0.9999;
    lu_table_pair_upper_r_bounds[2] = 0.9999;
    lu_table_pair_upper_r_bounds[1] = 0.9999;
    lu_table_pair_upper_r_bounds[0] = 0.9999;

}


#endif