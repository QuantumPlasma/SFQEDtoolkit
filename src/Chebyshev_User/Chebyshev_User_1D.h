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
    double evaluate(double) const;

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