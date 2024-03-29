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
#ifndef CHEBYSHEV_USER_2D
#define CHEBYSHEV_USER_2D

#include <fstream>

struct Chebyshev_User_2D {

public:

    //default constructor setting the SELF communicator as default
    Chebyshev_User_2D(){
        
    }

    Chebyshev_User_2D(unsigned int, double, double, unsigned int, double, double);

	//this pointer is meant to store the array of chebyshev coeffs
	double* last_coeffs;

	//double precision floating point nums storing the domain's extrema [a,b] of the first variable
	double a, b;
    //middle point of the domain (a+b)/2
	double domain_middle_point_N;
    //half range of the domain (b-a)/2
	double domain_half_range_N;

    //same as the block above, but this refers to the second variable
	double c, d;
	double domain_middle_point_M;
	double domain_half_range_M;

    //these numbers are used inside the evaluation function
	//and determine the "effective" orders of the approx
	//(basically the numbers of coefficients that are kept)
	unsigned int evaluation_order_N;
	unsigned int evaluation_order_M;

    //these store the full lengths (along N and M) of last_coeffs array.
    //NOTICE THAT they refer to the ACTUAL approximation order, and not to the number of
    //coefficients that have been kept (like evaluation_order_N/M do). This int is needed
    //to keep track of the approx order while using the various algorithms.
    //As far as the user is concerned, these 2 will coincide with the evaluation orders!
    unsigned int last_computed_order_N;
    unsigned int last_computed_order_M;

    //this method contains the 2D Clenshaw recurrence formula
    //that computes the approximation at a given point (x, y)
    //laying in the domain [a,b] x [c,d].
    double evaluate(double x, double y) const;

	//initialize approximation obeject from txt file
	static Chebyshev_User_2D init_from_txt_file(std::ifstream&);

    //these functions allow to evaluate the Chebyshev approximation 
    //coefficients at a certain x/y only, collapsing the 2D N x M matrix
    //to a 1D array of length N/M, that is then returned by the func
    double* evaluate_x(double) const;
    double* evaluate_y(double) const;

};


#endif