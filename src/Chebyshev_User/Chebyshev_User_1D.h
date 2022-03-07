#pragma once
#ifndef CHEBYSHEV_USER_1D
#define CHEBYSHEV_USER_1D

#include <fstream>
#include <mpi.h>

/*
 * This object stores the information needed to
 * describe the Chebyshev approximation of a given function,
 * just like the domain of validity ([a,b]) of the approximation
 * and the coefficients. 
 * It also implements the method to evaluate the approx at
 * any point (x \belong [a,b]).
*/
struct Chebyshev_User_1D {

public:

	//default constructor
	Chebyshev_User_1D(){
		//do nothing, this is needed just to ensure the availability of the default constructor
	}

	//constructor expecting the number of coefficients and the domain's extrema
	//(PAY ATTENTION: you are supposed to set the coefficients by hand)
	Chebyshev_User_1D(unsigned int, double, double);

	//this pointer is meant to store the array of chebyshev coeffs
	double* last_coeffs;

	//double precision floating point nums storing the domain's extrema [a,b]
	double a, b;
	//middle point of the domain (a+b)/2
	double domain_middle_point;
	//half range of the domain (b-a)/2
	double domain_half_range;

	//this number is used inside the evaluation function
	//and determines the "effective" order of the approx
	//(basically the number of coefficients kept)
	unsigned int evaluation_order;

    //this method contains the Clenshaw recurrence formula
    //that computes the approximation at a given point x
    //laying in the domain [a,b].
    virtual double evaluate(double);

	//initialize approximation object from txt file
	static Chebyshev_User_1D init_from_txt_file(std::ifstream&);

	//initialize approximation object from binary file (parallelized)
	//Notice this function will alter the file view.
	static Chebyshev_User_1D init_from_bin_file(MPI_File&, MPI_Comm&);

	//if you decide to implement the destructor
	//pay attention to how you use these objs: be sure
	//you always provide them with a last_coeffs
	//dynamically initialized
	// ~Chebyshev_User_1D(){
	// 	delete[] last_coeffs;
	// }

};

#endif