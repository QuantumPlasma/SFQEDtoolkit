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

#include <mpi.h>
#include <fstream>

struct Chebyshev_User_2D {

public:

    //default constructor setting the SELF communicator as default
    Chebyshev_User_2D(){
        MPI_Comm tmpComm;
        MPI_Comm_dup(MPI_COMM_SELF, &tmpComm);
        evalcom = &tmpComm;
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

	//communicator gathering the procs that share the workload in the evaluate
	MPI_Comm* evalcom;

    //this method contains the 2D Clenshaw recurrence formula
    //that computes the approximation at a given point (x, y)
    //laying in the domain [a,b] x [c,d].
    inline double __attribute__((always_inline)) evaluate(double x, double y) const{
        //variable conversion to -1 <= X <= +1 and -1 <= Y >= +1

        double X = (x - domain_middle_point_N) / domain_half_range_N;

        double Y = (y - domain_middle_point_M) / domain_half_range_M;

        //Clenshaw recurrence algorithm (easy peasy)
        //this stores the actual value that will be used in the algorithm
        double x_alg = X * 2.;
        double y_alg = Y * 2.;

        //define the first 2 coefficients of the Clenshaw series to be 0.
        //(inside the algorithm loop we will update these two double variables,
        //whose names should be meaningful)
        double tmp, clensh_j_plus_1 = 0., clensh_j = 0.;

        //Clenshaw coefficients
        double* d_i = new double[evaluation_order_N];

        /******************************/
        /* WORKLOAD SPLITTING SESSION */
        /******************************/
        //notice that we are just giving the possibility to evaluate
        //the approx by parallelizing the computation. Of course, the
        //splitting occurs only in the row of the coefficients matrix.
        //If you do not want to parallelize this, just set the MPI_COMM_SELF in the evalcom
    /*
        int rank;
        MPI_Comm_rank(*evalcom, &rank);
        int nProcs;
        MPI_Comm_size(*evalcom, &nProcs);

        //vectors used to store the counts and displacements informations
        int* fcounts = new int[nProcs];
        int* fdispls = new int[nProcs];
        fdispls[0] = 0;
        fcounts[0] = evaluation_order_N / nProcs + (0 < (evaluation_order_N % nProcs));

        //the following for section cannot be parallelized (omp),
        //as each iteration employs the result of the one
        //supposedly ended right before
        for (int i = 1; i < nProcs; i++) {
            fcounts[i] = evaluation_order_N / nProcs + (i < (evaluation_order_N% nProcs));
            fdispls[i] = fcounts[i - 1] + fdispls[i - 1];
        }
    */
        int loop_start = 0; //fdispls[rank];
        int loop_end = evaluation_order_N; //loop_start + fcounts[rank];

        //start computation

        //new coefficients determination loop (that on y)
        #pragma omp parallel for default(shared) private(tmp) firstprivate(clensh_j_plus_1, clensh_j)
        for (int i = loop_start; i < loop_end; i++) {

            //first Clenshaw algorithm
            for (int j = evaluation_order_M - 1; j >= 0; j--) {
                tmp = clensh_j;
                //don't be mistaken, you have to use approx_infos.last_computed_order_M
                clensh_j = last_coeffs[i*last_computed_order_M + j] - clensh_j_plus_1 + y_alg * clensh_j;
                clensh_j_plus_1 = tmp;
            }

            //assign new coefficient (remember to use Y)
            d_i[i] = clensh_j - Y * clensh_j_plus_1;

            //Reset the first 2 Clenshaw coefficient at each outer iteration
            //this must be here in case omp is not available
            clensh_j_plus_1 = 0.;
            clensh_j = 0.;
        }

        //TRANSMIT d_i computation to all procs
        //MPI_Allgatherv(MPI_IN_PLACE, fcounts[rank], MPI_DOUBLE, d_i, fcounts, fdispls, MPI_DOUBLE, *evalcom);

        //second Clenshaw algorithm application (on x)
        //ATTENTION: not parallelizable
        for (int i = evaluation_order_N - 1; i >= 0; i--) {
            tmp = clensh_j;
            clensh_j = d_i[i] - clensh_j_plus_1 + x_alg * clensh_j;
            clensh_j_plus_1 = tmp;
        }

        delete[] d_i;
        //delete[] fcounts;
        //delete[] fdispls;

        //remember to use X
        return clensh_j - X * clensh_j_plus_1;
    }

	//initialize approximation obeject from txt file
	static Chebyshev_User_2D init_from_txt_file(std::ifstream&);

    //initialize approximation object from binary file (parallelized)
	//Notice this function will alter the file view.
	static Chebyshev_User_2D init_from_bin_file(MPI_File&, MPI_Comm&);

    //these functions allow to evaluate the Chebyshev approximation 
    //coefficients at a certain x/y only, collapsing the 2D N x M matrix
    //to a 1D array of length N/M, that is then returned by the func
    double* evaluate_x(double) const;
    double* evaluate_y(double) const;

};


#endif