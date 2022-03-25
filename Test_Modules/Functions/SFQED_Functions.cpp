#include "SFQED_Functions.h"

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>  // error handling, for root finding
#include <gsl/gsl_math.h>  //evaluation of functions, for root finding
#include <gsl/gsl_roots.h>  //root finding routines
#include <gsl/gsl_math.h>

const double max_w = 8.88;
const double max_v = 1.;
const double pigreco = 3.141592653589793;

//this constant will be used to retrieve the value u_star of
//the variable u above which the pair production function
//pp_func(u,k) = 1/sqrt(u*(u-1))*BesselK[2/3,8*u/(3*k)]
//may be considered as zero. It has been chosen in such a
//way that the function in question lays below 1e-300 from
//u_star on.
//Its determination follows from considering BesselK[2/3,8*u/(3*k)]
//as the dominating part of the function. Thus, working on
//the case where k = 2000, it is easy to find that
//pp_func(u, 2000) < 1e-300 for u > 5.1e5. Now, since the
//argument of the dominating BesselK is 8*u/(3*k)
//we establish 5.1e5/2.0e3 to be the proportionality factor
//above which pp_func(u,k) decays to zero
const double pair_production_pre_constant = (5.1e5/2.0e3);
const double pair_production_max_constant = pair_production_pre_constant * (8./3.);

//integrand of tildeWrad
double integrand(double v, void* params){
	double chi = *(double*)params;

	double vchi = v * chi;
	double inv_den = 1.0 / (2.0 + 3. * vchi);

	//BesselK(2/3, x) diverges as x^-(2/3) for x->0
	//and goes in undeflow for x > about 700 because < 4.6712580780e-306
	double result = (20.0 + (42.0 + 45.0 * vchi) * vchi) * (inv_den * inv_den * inv_den)
		* gsl_sf_bessel_Knu(2.0 / 3.0, v);

	return result;
}

double tildeWrad(double X, void* pars){
	//the following line is needed only in case we use a X variable such that
	// -1 <= X <= +1. In that case we would need the conversion to a proper chi =>  a <= chi <= b
	//double chi = X; // ((b - a) * X + a + b) * 0.5;

	double result, error;

	gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);

	gsl_function Fintegrand;
	Fintegrand.function = &integrand;
	Fintegrand.params = &X;

	//BesselK[2/3, x] undeflows for x \gtrsim 700 because < 10^-306
	gsl_integration_qags(&Fintegrand, 0.0, 7.0e2, 1.0e-12, 0.0, 10000, w,
		&result, &error);

	gsl_integration_workspace_free(w);

	return result;
}


/***********************************/
/*      SOME NOTES ON MODIFIED 
       BESSEL FUNCTIONS OF THE
	          2ND KIND             */
/***********************************/
//sinxe we know that that K_{n}(x):
//-for x >> n
//		~ e^{-x} / sqrt(x*2./\pi)
//-for x = 0 and n = 0
//      ~ inf
//-for x = 0 and n != 0
//      ~ complex inf

//this function is the regularized version of w^2 * BesselK[2/3, w^3]
//the limit for small w has been obtained with wolfram mathematica through
//w -> 0 behaviour = Series[w^2 * BesselK[2/3, w^3], {w, 0, 5}, Assumptions -> w \[Element] PositiveReals] // Normal
//while, quite similarly, the limit for large w comes from
//w -> inf behaviour = s = Series[(1/w)^2*BesselK[2/3, (1/w)^3], {w, 0, 5}, Assumptions -> w \[Element] PositiveReals] // Normal; s /. w -> 1/x // FullSimplify
double w2BesselK23w3(double w)
{
	double result;
	double w2 = w * w;
	double w3 = w * w2;

	if (w <= 3.0e-3)
	{
		//no singularity in the function
		//avoid to evaluate BesselK[2/3, w^3] for w->0
		result = 1.074764120767239318355146 - 1.265719144219806941922564 * w2 * w2;
	}
	else if (w >= 8.88)
	{
		//undeflow for w>>1 because of the fast decay
		result = 0.0;
	}
	else
	{ //safe zone
		result = w2 * gsl_sf_bessel_Knu(2.0 / 3.0, w3);
	}

	return result;
}


//same reasoning as in w2BesselK23w3, but this time the function to approximate is BesselK[1/3, z^{3/2}] * \sqrt{z} * 3 / 2
//we used wolfram mathematica in the same way as before
double intBesselK13(double z, void* params)
{
	double result;

	if (z <= 1.0e-6)
	{
		//no singularity in the function
		//avoid to evaluate BesselK[1/3, z^{3/2}] for z->0
		result = 2.531438288439613883845128 - 2.418219271726288466299078 * z;
	}
	else if (z >= 78.84)
	{
		//undeflow for z>>1 because of the fast decay
		result = 0.0;
	}
	else
	{ //safe zone
		double sqrtz = sqrt(z);
		double sqrtz3 = sqrtz * sqrtz * sqrtz;
		result = 1.5 * sqrtz * gsl_sf_bessel_Knu(1.0 / 3.0, sqrtz3);
	}

	return result;
}

//returns: \int_{x}^{700}{BesselK[1/3, y] dy}
//note that: \int_{700}^{\inf}{BesselK[1/3, y] dy} < 4.6*10^-306
//in fact: BesselK[1/3, y] < 4.67*10^-306 for y > 700
double IBesselK13(double x)
{
	//for better convergency change of variable in the integrand
	// y = z^3/2 in order to remove the singularity
	double chi = 0.0;
	//ATTENTION!
	double cbrtx = cbrt(x);
	double z = cbrtx* cbrtx;

	if(z >= 78.84) return 0.;

	int status;

	double result, error;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);

	gsl_function Fintegrand;
	Fintegrand.function = &intBesselK13;
	Fintegrand.params = &chi;

	//BesselK[1/3, x] undeflows for x \gtrsim 700 because < 10^-306
	status = gsl_integration_qag(&Fintegrand, z, 78.84, 1.0e-12, 0.0, 10000,
		GSL_INTEG_GAUSS61, w, &result, &error);//eventually replace x -> z
	//if (status) std::cerr << "IBesselK13 error gsl_integration_qag" << std::endl;

	gsl_integration_workspace_free(w);

	return result;
}

//returns the photon emission differential probability (using the regularized
//version of the BesselK combinations given above)
double intChebinw(double w, void* params){
	
    double chi = *(double*)params;

	double w2 = w * w;
	double w3 = w * w2;
	double coef = 1.0 + 1.5 * chi * w3;
	double coef2 = coef * coef;
	double coef3 = coef * coef2;

	double result;

	result = 4.5 * (w2BesselK23w3(w) * (1.0 + coef2) - w2 * coef * IBesselK13(w3)) / coef3;

	return result;
}

double differential_cross_section_phtn_emission(double chi, double w, void* params){
	return intChebinw(w, &chi);
}

//integral of the photon emission differential probability
//thus giving the photon emission rate
double IChebinw(double w, double chi){

	double result, error;
	int status;

	gsl_integration_workspace* w1 = gsl_integration_workspace_alloc(10000);

	gsl_function Fintegrand;
	Fintegrand.function = &intChebinw;
	Fintegrand.params = &chi;

	status = gsl_integration_qag(&Fintegrand, 0.0, w, 1.0e-12, 0.0, 10000,
		GSL_INTEG_GAUSS61, w1, &result, &error);
	//if (status) std::cerr << "IChebinw error gsl_integration_qag" << std::endl;

	gsl_integration_workspace_free(w1);

	return result;
}

//this simply wraps the IChebinw function above (notice I exchange the variables' order)
double Wrad_integral(double chi, double w, void *params) {

	return IChebinw(w, chi);
}

//function whose inverse will give w and thus the energy of the emitted photon
//being \epsilon_{\gamma} = \frac{3\chi w^3}{2 + 3\chi w^3}\epsilon_e
double equation(double w, void* params)
{
	double *p = static_cast<double*>(params);

	double chi = p[0];
	double r = p[1];

	double result = IChebinw(w, chi) - r * IChebinw(max_w, chi);

	return result;
}

//inverse of equation
double inverse(double chi, double r)
{
	double *params = new double[2]{ chi, r };

	//avoid to allocate and deallocate at each call of the function
	// shifted to => the "global" definitions of the program
	//
	//const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
	//gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);
	//gsl_function F;

	const gsl_root_fsolver_type* T = gsl_root_fsolver_brent;
	gsl_root_fsolver* s = gsl_root_fsolver_alloc(T);
	gsl_function F;

	F.function = &equation;
	F.params = params;

	int status;
	int iter = 0;
	const int max_iter = 100;

	double root;
	double x_lo = 0.0;
	double x_hi = max_w;

	gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	do
	{
		iter++;

		status = gsl_root_fsolver_iterate(s);

		root = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);

		status = gsl_root_test_interval(x_lo, x_hi, 1.0e-10, 0.0);

	} while (status == GSL_CONTINUE && iter < max_iter);

	//avoid to allocate and deallocate at each call of the function
	// shifted to => the "main" of the program
	//
	gsl_root_fsolver_free (s);

	delete[] params;

	return root;
}

//wrapper of inverse
double photon_w(double x, double y, void* params) {
	return inverse(x, y);
}

double rec_phtn_w(double x, double y, void* params){
	return 1 / inverse(x, y);
}

double photon_nrg(double chi, double rnd) {
	double w = inverse(chi, rnd);
	w *= w*w;

	w *= 3.0*chi;
    return w / (2.0 + w);
}

double photon_w_fixed_r(double chi, void *params) {
	double r = *static_cast<double*>(params);
	return inverse(chi, r);
}

//////////////////////////////////////////////////////

//this function is the regularized version of \frac{2(2 u - 1)}{(u (u - 1))^(1/2)} * BesselK[2/3, \frac{8u}{3k}]
//the lower limit for u has been obtained with wolfram mathematica through
//u -> 1 behaviour = Series[2*(2 u - 1)*(u (u - 1))^(-1/2)*BesselK[2/3, 8*u/(3*k) ], {u, 1, 2}, Assumptions -> u \[Element] PositiveReals] // Normal
//while, quite similarly, the limit for large u comes from considering
//


double c0(double k, void* params){
	return 2. * gsl_sf_bessel_Knu(2.0 / 3.0, 8. / (3. * k)) + IBesselK13(8. / (3. * k));
}

//this function has no problems at all when 
double pairProductionRateIntegrand(double v, void* params){
	double k = *(static_cast<double*>(params));

	double cut;
	if(k <= 20.){
		cut = 0.998;
	}
	else if(k <= 80.){
		cut = 0.9998;
	}
	else{
		cut = 0.99998;
	}

	if(v >= cut){
		return exp(-2./(3. * k) + 4./(3. * (-1. + v) * k) + (-1. + v)/(3. * k)) * sqrt(pigreco * 6. * k / (1. - v));
	}

	return (9. - v*v) / (1. - v*v) * gsl_sf_bessel_Knu(2.0 / 3.0, 8. / (3. * k * (1. - v*v)));
	
}

double pairProductionRate(double k, void* params){
	double result, error;

	double upper;
	/*if(k < 0.05){
		// return (27. * pigreco * k)/(16. * sqrt(2.)) * exp(-(8/(3 * k))) * (1. - 11./64.*k + 7585./73728.*k*k);
		upper = 0.7;
	}
	else if(k < 0.1){
		// return (27. * pigreco * k)/(16. * sqrt(2.)) * exp(-(8/(3 * k))) * (1. - 11./64.*k + 7585./73728.*k*k);
		upper = 0.8;
	}
	else*/ if(k < 0.4){
		upper = 0.9;
	}
	else if(k < 2.){
		upper = 0.99;
	}
	else if(k < 20.){
		upper = 0.999;
	}
	else if(k < 80.){
		upper = 0.9999;
	}
	else{
		upper = 0.99999;
	}

	gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);

	gsl_function Fintegrand;
	Fintegrand.function = &pairProductionRateIntegrand;
	Fintegrand.params = &k;
	
	//before we were trying with epsabs = 1.0e-12, epsrel = 0.0
	gsl_integration_qags(&Fintegrand, 0., upper, 1.0e-9, 1.0e-10, 10000, w, &result, &error);
	// gsl_integration_qag(&Fintegrand, 0, upper, 1.0e-9, 1.0e-10, 10000, GSL_INTEG_GAUSS61, w, &result, &error);

	gsl_integration_workspace_free(w);

	return result;
}

// double integrandPairPartialRateSmallk(double v, void* params){
	
//     double k = *(double*)params;
// 	double common_denominator = (1. - v*v);
// 	double eta = 8./(3. * k * common_denominator);

// 	return 0.5 * sqrt(1. / (common_denominator * pigreco * k)) * exp(-eta) * (3 + v*v);
// }

double integrandPairPartialRateSmallk(double v, void* params){
	
    double k = *(double*)params;
	double common_denominator = (1. - v*v);
	double eta = 8./(3. * k * common_denominator);

	double cut;
	
	if(k <= 0.05){
		cut = 0.69;
	}
	else if(k <= 0.1){
		cut = 0.79;
	}
	else{
		cut = 0.89;
	}

	double first_term;
	
	if(v >= cut){
		first_term = exp(-2./(3. * k) + 4./(3. * (-1. + v) * k) + (-1. + v)/(3. * k)) * sqrt(pigreco * 1.5 * k / (1. - v));
	}
	else {
		first_term = 2. * (1. + v*v) / common_denominator * gsl_sf_bessel_Knu(2.0 / 3.0, eta);
	}

	return first_term + IBesselK13(eta);
}

double pairPartialRateSmallk(double v_up, double k){

	if(k < 0.01){
	// if(k == 0.){
		return 0.;
	}

	double abserr = 1.0e-9,
			relerr = 1.0e-10;

	double result, error;
	int status;

	gsl_integration_workspace* w1 = gsl_integration_workspace_alloc(10000);

	gsl_function Fintegrand;
	Fintegrand.function = &integrandPairPartialRateSmallk;
	Fintegrand.params = &k;

	//new
	//these correspond to the transformation \epsilon = \epsilon_{\gamma} / 2 (1 + v)
	status = gsl_integration_qags(&Fintegrand, 0, v_up, abserr, relerr, 10000, w1, &result, &error);

	//if (status) std::cerr << "IChebinw error gsl_integration_qag" << std::endl;

	gsl_integration_workspace_free(w1);

	return result;
}

double pairPartialRateSmallk_wrapper(double chi, double v_low, void *params) {

	return pairPartialRateSmallk(v_low, chi);
}

double pp_equation_small_k(double lower_v, void* params){
	double *p = static_cast<double*>(params);

	double kappa = p[0];
	double r = p[1];

	//adjust this upper bound depending
	//on the value of k
	double maximum_v;

	if(kappa <= 0.05){
		maximum_v = 0.7;
	}
	else if(kappa <= 0.1){
		maximum_v = 0.8;
	}
	else{
		maximum_v = 0.9;
	}

	//new
	double result = pairPartialRateSmallk(lower_v, kappa) - r * pairPartialRateSmallk(maximum_v, kappa);

	return result;

}

double pp_inverted_small_k(double kappa, double r)
{

	double *params = new double[2]{ kappa, r };

	//avoid to allocate and deallocate at each call of the function
	// shifted to => the "global" definitions of the program
	//
	//const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
	//gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);
	//gsl_function F;

	const gsl_root_fsolver_type* T = gsl_root_fsolver_brent;
	gsl_root_fsolver* s = gsl_root_fsolver_alloc(T);
	gsl_function F;

	F.function = &pp_equation_small_k;
	F.params = params;

	int status;
	int iter = 0;
	const int max_iter = 100;

	double root;
	// double x_lo = 0.0;
	//if the extended version goes wrong, please comment the statement below and uncomment the one above
	// double x_lo = -upper;
	double x_lo = -0.01;
	double x_hi = 1.; //upper;

	gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	do
	{
		iter++;

		status = gsl_root_fsolver_iterate(s);

		root = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);

		status = gsl_root_test_interval(x_lo, x_hi, 1.0e-10, 0.0);

	} while (status == GSL_CONTINUE && iter < max_iter);

	//avoid to allocate and deallocate at each call of the function
	// shifted to => the "main" of the program
	//
	gsl_root_fsolver_free (s);

	delete[] params;

	return root;
}

double pair_production_inv_wrapper_small_k(double kappa, double r, void* params) {
	return pp_inverted_small_k(kappa, r);
}

double integrandPairPartialRate(double v, void* params){
	
    double k = *(double*)params;
	double common_denominator = (1. - v*v);
	double eta = 8./(3. * k * common_denominator);

	double cut;
	//move the cut in the range 0.3 -> 2
	// if(k <= 0.4){
	// 	cut = 0.89;

	// 	// if(v > 0.9){
	// 	// 	return 0.;
	// 	// }
	// }
	if(k <= 2.){
		cut = 0.98;

		// if(v > 0.99){
		// 	return 0.;
		// }
	}
	else if(k <= 20.){
		//with cut = 0.999 the approximation is never used
		cut = 0.998;
	}
	else if(k <= 80.){
		cut = 0.9998;
	}
	else if(k <= 600){
		//before we tried with 0.99998
		cut = 0.99993;
	}
	else{
		cut = 0.99998;
	}

	double first_term;
	
	if(v >= cut){
		first_term = exp(-2./(3. * k) + 4./(3. * (-1. + v) * k) + (-1. + v)/(3. * k)) * sqrt(pigreco * 1.5 * k / (1. - v));
	}
	//if the extended version goes wrong, please comment the if below and uncomment the one above
	// if(v >= cut || v <= -cut){
	// 	first_term = exp(-2./(3. * k) + 4./(3. * (-1. + std::abs(v)) * k) + (-1. + std::abs(v))/(3. * k)) * sqrt(pigreco * 1.5 * k / (1. - std::abs(v)));
	// }
	else {
		first_term = 2. * (1. + v*v) / common_denominator * gsl_sf_bessel_Knu(2.0 / 3.0, eta);
	}

	return first_term + IBesselK13(eta);
}

double pairPartialRate(double v_low, double k){

	double upper,
			abserr = 1.0e-9,
			relerr = 1.0e-10;

	if(k < 0.1){
		upper = 0.8;
	}
	else if(k < 0.4){
		upper = 0.9;
	}
	else if(k < 2.){
		upper = 0.99;
	}
	else if(k < 20.){
		upper = 0.999;
	}
	else if(k < 80.){
		upper = 0.9999;
	}
	else{
		upper = 0.99999;

		//we require a less stringent error threshold
		//in the 600 < k < 2000 case
		if(k >= 600){
			abserr = 1.0e-8;
			relerr = 1.0e-9;
		}
		
	}

	double result, error;
	int status;

	gsl_integration_workspace* w1 = gsl_integration_workspace_alloc(10000);

	gsl_function Fintegrand;
	Fintegrand.function = &integrandPairPartialRate;
	Fintegrand.params = &k;

	//old
	//these correspond to the transformation \epsilon = \epsilon_{\gamma} / 2 (1 - v)
	// status = gsl_integration_qag(&Fintegrand, v_low, upper, 1.0e-9, 1.0e-10, 10000, GSL_INTEG_GAUSS61, w1, &result, &error);
	// status = gsl_integration_qags(&Fintegrand, v_low, upper, abserr, relerr, 10000, w1, &result, &error);

	//new
	//these correspond to the transformation \epsilon = \epsilon_{\gamma} / 2 (1 + v)
	status = gsl_integration_qags(&Fintegrand, 0, v_low, abserr, relerr, 10000, w1, &result, &error);

	//if (status) std::cerr << "IChebinw error gsl_integration_qag" << std::endl;

	gsl_integration_workspace_free(w1);

	return result;
}

//this simply wraps the IChebinw function above (notice I exchange the variables' order)
double pairPartialRate_wrapper(double chi, double v_low, void *params) {

	return pairPartialRate(v_low, chi);
}

double pp_equation(double lower_v, void* params){
	double *p = static_cast<double*>(params);

	double kappa = p[0];
	double r = p[1];

	//be careful to the last argument of this pp_diff_cross_section function (and use the full rate for the second func)
	//and remember to include the 2 factor!
	//old
	// double result = 3. * pairPartialRate(lower_v, kappa) - 2. * r * pairProductionRate(kappa, NULL);

	//new
	double result = 3. * pairPartialRate(lower_v, kappa) - r * pairProductionRate(kappa, NULL);

	return result;

}

//inverse of equation
double pp_inverted(double kappa, double r)
{
	
	double upper;
	if(kappa < 0.1){
		upper = 0.8;
	}
	else if(kappa < 0.4){
		upper = 0.9;
	}
	else if(kappa < 2.){
		upper = 0.99;
	}
	else if(kappa < 20.){
		upper = 0.999;
	}
	else if(kappa < 80.){
		upper = 0.9999;
	}
	else{
		upper = 0.99999;
	}

	double *params = new double[2]{ kappa, r };

	//avoid to allocate and deallocate at each call of the function
	// shifted to => the "global" definitions of the program
	//
	//const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
	//gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);
	//gsl_function F;

	const gsl_root_fsolver_type* T = gsl_root_fsolver_brent;
	gsl_root_fsolver* s = gsl_root_fsolver_alloc(T);
	gsl_function F;

	F.function = &pp_equation;
	F.params = params;

	int status;
	int iter = 0;
	const int max_iter = 100;

	double root;
	// double x_lo = 0.0;
	//if the extended version goes wrong, please comment the statement below and uncomment the one above
	// double x_lo = -upper;
	double x_lo = -0.01;
	double x_hi = 1.; //upper;

	gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	do
	{
		iter++;

		status = gsl_root_fsolver_iterate(s);

		root = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);

		status = gsl_root_test_interval(x_lo, x_hi, 1.0e-10, 0.0);

	} while (status == GSL_CONTINUE && iter < max_iter);

	//avoid to allocate and deallocate at each call of the function
	// shifted to => the "main" of the program
	//
	gsl_root_fsolver_free (s);

	delete[] params;

	return root;
}

//wrapper of inverse
double pair_production_inv_wrapper(double kappa, double r, void* params) {
	return pp_inverted(kappa, r);
}

double brent_inverse_test(double kappa, double r, void* params) {
	return 1. * 0.5 * (1. + pp_inverted(kappa, r));
}

double low_approx_inverse_test(double kappa, double r, void* params) {
	return 0.5 + (r / c0(kappa, NULL) * pairProductionRate(kappa, NULL) / 3.);
}

double esponential_tail_test(double kappa, double r, void* params) {
	double limit_r = 0.995;
	double v0_pivot = pp_inverted(kappa, limit_r);

	double partial_pivot_rate = pairPartialRate(v0_pivot, kappa);
	double tildeW = pairProductionRate(kappa, NULL) / 3.;
	double ln = log((1. - r) * tildeW / (tildeW - partial_pivot_rate));

	v0_pivot = 1. / (1. - v0_pivot * v0_pivot);
	double v = sqrt(1. - 1./(v0_pivot - 3. * kappa * ln / 8.));

	return 0.5 * (1. + v);
}

double electron_nrg(double chi, double rnd) {
	
	double rescaled_rnd = 2. * rnd - 1.;
    double sgn = -2. * std::signbit(rescaled_rnd) + 1.;
    //extract the absolute value
    rescaled_rnd *= sgn;
	
	return 0.5 * (1. + sgn * pp_inverted(chi, rescaled_rnd));
}

double electron_nrg_small_k(double chi, double rnd) {
	
	double rescaled_rnd = 2. * rnd - 1.;
    double sgn = -2. * std::signbit(rescaled_rnd) + 1.;
    //extract the absolute value
    rescaled_rnd *= sgn;
	
	return 0.5 * (1. + sgn * pp_inverted_small_k(chi, rescaled_rnd));
}

//////////////////////
// various attempts //
//////////////////////

//attempt 0
double pairCoeffBesselK23Xi(double u, double k){
	double ph = 8./3./k;
	double xi = ph * u;

	//this is the minimum value of u below which
	//the relative error between the actual function
	//and the corresponding polynomial approximation is less
	//than 1e-10 for every sampled value of k (namely for k \in [0.3,2000]) 
	if(u <= 1.00008){
		double K23 = gsl_sf_bessel_Knu(2.0 / 3.0, ph);
		double K53 = gsl_sf_bessel_Knu(5.0 / 3.0, ph);
		double sqrt_u_minus_one = sqrt(u - 1.);

		return K23 / sqrt_u_minus_one +
				(k*K23 - 16.*K53) * sqrt_u_minus_one / 6. / k + 
				(256.*K23 - 5.*k*k*K23 + 192.*K53) * sqrt_u_minus_one*sqrt_u_minus_one*sqrt_u_minus_one / 72. / k / k;
	}
	else if(xi >= pair_production_max_constant){
		//this treats the underflow
		return 0.;
	} else {
		return gsl_sf_bessel_Knu(2.0 / 3.0, xi) / sqrt(u*(u-1.));
	}
}

double wildPairCoeffBesselK23Xi(double u, double k){

	if(k == 0.){
		//this treats the underflow
		return 0.;
	} else {
		return gsl_sf_bessel_Knu(2.0 / 3.0, 8. * u / 3. / k) / sqrt(u*(u-1.));
	}
}

//pair production rate integrand
double pair_production_rate_integrand(double u, void *params){
	double kappa = *static_cast<double*>(params);

	double coeff = 8. + 1./u;
	
	return coeff * wildPairCoeffBesselK23Xi(u, kappa); //pairCoeffBesselK23Xi(u, kappa);
}

//pair production rate (whole spectrum)
//remember to approximate this in the range k > 3
double pair_production_rate_full_spectrum(double kappa, void* params){
	double result, error;

	gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);

	gsl_function Fintegrand;
	Fintegrand.function = &pair_production_rate_integrand;
	Fintegrand.params = &kappa;

	//the upper integration bound is determined using the constant
	//defined at the beginning of the module. We do this as function
	//pairCoeffBesselK23Xi (look at its definition) is zero from a certain u on
	double upper_bound = pair_production_pre_constant * kappa;
	gsl_integration_qags(&Fintegrand, 1., upper_bound, 1.0e-12, 0.0, 10000, w, &result, &error);

	gsl_integration_workspace_free(w);

	return result;
}

double pp_diff_cross_section_integrand(double u, void *params){
	double kappa = *static_cast<double*>(params);

	return kappa == 0. ? 0. : (4.*u - 2.) * wildPairCoeffBesselK23Xi(u, kappa) + IBesselK13(8.*u / (3.*kappa));
}

double pp_diff_cross_section(double lower_u, double kappa, void* params){


	double some_upper_val = *(static_cast<double*>(params));

	// if(kappa == 0.){
	// 	return 0.;
	// }
	// else if(kappa <= 2.){
	// 	some_upper_val = 60.;
	// }
	// else if(kappa <= 20.){
	// 	some_upper_val = 245.;
	// }
	// else if(kappa <= 80.){
	// 	some_upper_val = 1050.;
	// }
	// else if(kappa <= 600.){
	// 	some_upper_val = 8000.;
	// }
	// else{
	// 	//this correspond to an error of 2e-12
	// 	some_upper_val = 25000.;
	// }


	double result, error;
	int status;

	gsl_integration_workspace* w1 = gsl_integration_workspace_alloc(10000);

	gsl_function Fintegrand;
	Fintegrand.function = &pp_diff_cross_section_integrand;
	Fintegrand.params = &kappa;

	// status = gsl_integration_qagiu(&Fintegrand, lower_u,  1.0e-12, 0.0, 10000, w1, &result, &error);
	status = gsl_integration_qag(&Fintegrand, lower_u, some_upper_val, 1.0e-12, 0.0, 10000, GSL_INTEG_GAUSS61, w1, &result, &error);
	//if (status) std::cerr << "IChebinw error gsl_integration_qag" << std::endl;

	gsl_integration_workspace_free(w1);

	return result;
}

double pp_equation_to_invert(double lower_u, void* params){
	double *p = static_cast<double*>(params);

	double kappa = p[0];
	double r = p[1];

	//be careful to the last argument of this pp_diff_cross_section function (and use the full rate for the second func)
	double result = pp_diff_cross_section(lower_u, kappa, NULL) - r * pp_diff_cross_section(1., kappa, NULL);

	return result;

}

//inverse of equation
double pp_inverse(double kappa, double r)
{
	double *params = new double[2]{ kappa, r };

	//avoid to allocate and deallocate at each call of the function
	// shifted to => the "global" definitions of the program
	//
	//const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
	//gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);
	//gsl_function F;

	const gsl_root_fsolver_type* T = gsl_root_fsolver_brent;
	gsl_root_fsolver* s = gsl_root_fsolver_alloc(T);
	gsl_function F;

	F.function = &pp_equation_to_invert;
	F.params = params;

	int status;
	int iter = 0;
	const int max_iter = 100;

	double root;
	double x_lo = 0.0;
	double x_hi = max_w;

	gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	do
	{
		iter++;

		status = gsl_root_fsolver_iterate(s);

		root = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);

		status = gsl_root_test_interval(x_lo, x_hi, 1.0e-10, 0.0);

	} while (status == GSL_CONTINUE && iter < max_iter);

	//avoid to allocate and deallocate at each call of the function
	// shifted to => the "main" of the program
	//
	gsl_root_fsolver_free (s);

	delete[] params;

	return root;
}

//wrapper of inverse
double pair_production_u(double kappa, double r, void* params) {
	return pp_inverse(kappa, r);
}

//////////////////////////////

//attempts 1
double prova(double u, void* params){
	double k = *(static_cast<double*>(params));

	if(k <= 0.17){
		if(k <= 0.01) return 0.;

		return (8. + 1./u) * (
								(-455./98304. * sqrt(k*k*k*k*k * pigreco / (u*u*u*u*u * 3)))
								+ (7./256. * sqrt(k*k*k * pigreco / (3 * u*u*u)))
								+ (0.25 * sqrt(k * 3. * pigreco / u))
							)
				/ (exp(8. * u / (3. * k)) * sqrt(u*(u-1.)));
	}

	return (8. + 1./u) * gsl_sf_bessel_Knu(2.0 / 3.0, 8. * u / 3. / k) / sqrt(u*(u-1.));
	
}

double provaIntegrata(double upper, double k){
	double result, error;

	gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);

	gsl_function Fintegrand;
	Fintegrand.function = &prova;
	Fintegrand.params = &k;

	// gsl_integration_qags(&Fintegrand, 1., upper, 1.0e-12, 0.0, 10000, w, &result, &error);
	gsl_integration_qag(&Fintegrand, 1, upper, 1.0e-12, 0.0, 10000, GSL_INTEG_GAUSS61, w, &result, &error);

	gsl_integration_workspace_free(w);

	return result;
}
////////////////////

//attempts 2
double pairFullRateIntegrand(double xi, void* params){
	double k = *(static_cast<double*>(params));

	double coeff = 8. * (3. * k * xi + 1.) / xi / sqrt(3. * k * xi * (3. * k * xi - 8.));

	return coeff * gsl_sf_bessel_Knu(2.0 / 3.0, xi);
	
}

double pairFullRate(double k, double upper_bound){
	double result, error;

	gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000);

	gsl_function Fintegrand;
	Fintegrand.function = &pairFullRateIntegrand;
	Fintegrand.params = &k;

	double lower_bound = 8. / 3. / k;

	gsl_integration_qags(&Fintegrand, lower_bound, upper_bound, 1.0e-12, 0.0, 10000, w, &result, &error);
	// gsl_integration_qag(&Fintegrand, lower_bound, 1000, 1.0e-12, 0.0, 10000, GSL_INTEG_GAUSS61, w, &result, &error);

	gsl_integration_workspace_free(w);

	return result;
	
}
///////////////////