#pragma once
#ifndef SFQED_FUNCS
#define SFQED_FUNCS

/*****************/
/*PHOTON EMISSION*/
/*****************/

//1-dim functions

//this function describe the integral contributing to the TOTAL
//photon emission rate W_{rad} (with reference to the notes).
//It depends only on chi and is the one
//that is denoted as \tilde{W}_{rad} in the notes
double tildeWrad(double, void*);

//one dimensional version of photon_w, where r is fixed by the params
double photon_w_fixed_r(double, void*);

//2-dim functions

//photon emission differential cross section
//this basically returns the integrand of Wrad_integral
double differential_cross_section_phtn_emission(double, double, void*);

//rate of photon emission of an electron having a certain amount of energy
double Wrad_integral(double, double, void*); 

//inverse of the g(\chi, w) - r*W_{rad}(\chi) = 0 equation
double photon_w(double, double, void*);

double rec_phtn_w(double, double, void*);

double photon_nrg(double, double);

/*****************/
/*PAIR PRODUCTION*/
/*****************/

double c0(double, void*);

double pairProductionRate(double, void*);

double pairPartialRate_wrapper(double, double, void *);

double pairPartialRateSmallk_wrapper(double, double, void *);

double pair_production_inv_wrapper(double, double, void*);

double pair_production_inv_wrapper_small_k(double, double, void*);


//debug purpose
double brent_inverse_test(double, double, void*);
double low_approx_inverse_test(double, double, void*);
double esponential_tail_test(double, double, void*);


//////////////////////////////////////////
// function needed for accuracy testing //
//////////////////////////////////////////
double photon_nrg(double chi, double rnd);
double electron_nrg(double chi, double rnd);
double electron_nrg_small_k(double chi, double rnd);

//////////////////
//failed attempts
double prova(double, void*);
double provaIntegrata(double, double);

//pair production rate (whole spectrum). it is the equivalent of tildeWrad above
//remember to approximate this in the range k > 3
double pair_production_rate_full_spectrum(double, void*);

double pairFullRate(double, double);

#endif