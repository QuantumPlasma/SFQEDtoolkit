#pragma once
#ifndef CHEBYSHEV_USER_3D
#define CHEBYSHEV_USER_3D

#include <mpi.h>
#include <fstream>

struct Chebyshev_User_3D {
	
public:

    Chebyshev_User_3D(){
        MPI_Comm tmpComm;
        MPI_Comm_dup(MPI_COMM_SELF, &tmpComm);
        evalcom = &tmpComm;
    }

    Chebyshev_User_3D(unsigned int, double, double, unsigned int, double, double, unsigned int, double, double);

	//this pointer is meant to store the array of chebyshev coeffs
	double* last_coeffs;

    //these store the full lengths (along N, M and L) of last_coeffs array.
    //NOTICE THAT they refer to the ACTUAL approximation orders, and not to the number of
    //coefficients that have been kept (like evaluation_order_N/M/L do). This int is needed
    //to keep track of the approx order while using the various algorithms.
    //As far as the user is concerned, these 2 will coincide with the evaluation orders!
	unsigned int last_computed_order_N;
	unsigned int last_computed_order_M;
	unsigned int last_computed_order_L;

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

    //same as the block above, but this refers to the third variable
	double e, f;
	double domain_middle_point_L;
	double domain_half_range_L;

    //these numbers are used inside the evaluation function
	//and determine the "effective" orders of the approx
	//(basically the numbers of coefficients that are kept)
	unsigned int evaluation_order_N;
	unsigned int evaluation_order_M;
	unsigned int evaluation_order_L;

	//communicator gathering the procs that share the workload in the evaluate
	MPI_Comm* evalcom;

    //this method contains the 3D Clenshaw recurrence formula
    //that computes the approximation at a given point (x, y, z)
    //laying in the domain [a,b] x [c,d] x [e,f].
    double evaluate(double, double, double);

	//initialize approximation obeject from txt file
	static Chebyshev_User_3D init_from_txt_file(std::ifstream&);

    //initialize approximation object from binary file (parallelized)
	//Notice this function will alter the file view.
	static Chebyshev_User_3D init_from_bin_file(MPI_File&, MPI_Comm&);
};

#endif