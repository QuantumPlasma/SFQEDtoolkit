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

#include "Chebyshev_User_2D.h"

#include <omp.h>
#include <iostream>

Chebyshev_User_2D::Chebyshev_User_2D(unsigned int order_N, double a, double b, unsigned int order_M, double c, double d): Chebyshev_User_2D(){
    this->evaluation_order_N = order_N;
    this->last_computed_order_N = order_N;
    this->a = a;
    this->b = b;
    this->domain_half_range_N = (b-a)*0.5;
    this->domain_middle_point_N = (b+a)*0.5;

    this->last_computed_order_M = order_M;
    this->evaluation_order_M = order_M;
    this->c = c;
    this->d = d;
    this->domain_half_range_M = (d-c)*0.5;
    this->domain_middle_point_M = (d+c)*0.5;
}

double Chebyshev_User_2D::evaluate(double x, double y) const{
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

//comments omitted (simplified version of evaluate)
double* Chebyshev_User_2D::evaluate_y(double y) const{

	double Y = (y - domain_middle_point_M) / domain_half_range_M;

	double y_alg = Y * 2.;

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

	//delete[] fcounts;
    //delete[] fdispls;

	//size = evaluation_order_N
	return d_i;
}

//comments omitted (simplified version of evaluate)
double* Chebyshev_User_2D::evaluate_x(double x) const{

	double X = (x - domain_middle_point_N) / domain_half_range_N;

	double x_alg = X * 2.;

	double tmp, clensh_j_plus_1 = 0., clensh_j = 0.;

	//Clenshaw coefficients
	double* d_i = new double[evaluation_order_M];

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
	fcounts[0] = evaluation_order_M / nProcs + (0 < (evaluation_order_M % nProcs));

	//the following for section cannot be parallelized (omp),
	//as each iteration employs the result of the one
	//supposedly ended right before
	for (int i = 1; i < nProcs; i++) {
		fcounts[i] = evaluation_order_M / nProcs + (i < (evaluation_order_M% nProcs));
		fdispls[i] = fcounts[i - 1] + fdispls[i - 1];
	}
*/
	int loop_start = 0; //fdispls[rank];
	int loop_end = evaluation_order_M; //loop_start + fcounts[rank];

    //start computation

	//new coefficients determination loop (that on x) notice that i and j have been exchanged (wrt evaluate_y)
	#pragma omp parallel for default(shared) private(tmp) firstprivate(clensh_j_plus_1, clensh_j)
	for (int j = loop_start; j < loop_end; j++) {

		//first Clenshaw algorithm
		for (int i = evaluation_order_N - 1; i >= 0; i--) {
			tmp = clensh_j;
			//don't be mistaken, you have to use approx_infos.last_computed_order_M
			clensh_j = last_coeffs[i*last_computed_order_M + j] - clensh_j_plus_1 + x_alg * clensh_j;
			clensh_j_plus_1 = tmp;
		}

		//assign new coefficient (remember to use Y)
		d_i[j] = clensh_j - X * clensh_j_plus_1;

		//Reset the first 2 Clenshaw coefficient at each outer iteration
		//this must be here in case omp is not available
		clensh_j_plus_1 = 0.;
		clensh_j = 0.;
	}

    //TRANSMIT d_i computation to all procs
    //MPI_Allgatherv(MPI_IN_PLACE, fcounts[rank], MPI_DOUBLE, d_i, fcounts, fdispls, MPI_DOUBLE, *evalcom);

	//delete[] fcounts;
    //delete[] fdispls;
	
	//size = evaluation_order_M
	return d_i;
}

Chebyshev_User_2D Chebyshev_User_2D::init_from_txt_file(std::ifstream& in_file){
    unsigned int order_N, order_M;
    double a, b, c, d;

    in_file >> order_N >> a >> b
	        >> order_M >> c >> d;

    Chebyshev_User_2D cheby_user_2d(order_N, a, b, order_M, c, d);

	unsigned int order = order_N*order_M;
	cheby_user_2d.last_coeffs = new double[order];

	for (unsigned int i = 0; i < order; i++) {
		
		in_file >> cheby_user_2d.last_coeffs[i];
	}

    return cheby_user_2d;
    
}

Chebyshev_User_2D Chebyshev_User_2D::init_from_bin_file(MPI_File& file, MPI_Comm& readcom){
    int rank;
	MPI_Comm_rank(readcom, &rank);

	int nProcs;
	MPI_Comm_size(readcom, &nProcs);

	if (rank == 0) {
		std::cout << "Reading orders and domain from file...";
	}

    unsigned int order_N, order_M;
    double a, b, c, d;

	//let's read the infos obj related informations we have available on the file.
	//All the procs will do so. THIS WILL BE VERY UNELEGANT...
	MPI_File_read(file, &(order_N), 1, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_read(file, &(a), 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_read(file, &(b), 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_read(file, &(order_M), 1, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_read(file, &(c), 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_read(file, &(d), 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

    Chebyshev_User_2D cheby_user_2d(order_N, a, b, order_M, c, d);

	if (rank == 0) {
		std::cout << " Done!\nSplitting workload...";
	}

	int* counts = new int[nProcs];
	int* displs = new int[nProcs];

	unsigned int coeff_order = order_N * order_M;

	counts[0] = coeff_order / nProcs + (0 < (coeff_order % nProcs));
	displs[0] = 0;

	for (int i = 1; i < nProcs; i++) {
		counts[i] = (coeff_order / nProcs + (i < coeff_order % nProcs));
		displs[i] = (displs[i - 1] + counts[i - 1]);
	}

	int start = displs[rank];

	//we already read something so a pre offset has to be built accordingly
	int pre_offset = sizeof(unsigned int) * 2 + sizeof(double) * 4;

	//start should contain the distance in (double), thus we adjust it to be a proper offset
	MPI_Offset offset = sizeof(double) * start + pre_offset;

	MPI_File_set_view(file, offset, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

	//allocate the memory for our coefficients
	double* coefficients = new double[coeff_order];

	//in case you need to reset the file pointer
	//MPI_File_seek(file, offset, MPI_SEEK_SET);

	if (rank == 0) {
		std::cout << " Done!\nReading coefficients from file...";
	}

	//actual reading step
	MPI_File_read(file, &(coefficients[start]), counts[rank], MPI_DOUBLE, MPI_STATUS_IGNORE);

	//gather what a proc read to all other procs
	MPI_Allgatherv(MPI_IN_PLACE, counts[rank], MPI_DOUBLE, coefficients, counts, displs, MPI_DOUBLE, readcom);

	//pack the coeffs inside the infos obj
	cheby_user_2d.last_coeffs = coefficients;

	/********************* TESTING PART *******************
	//DO NOT DELETE THIS PART (it may come in handy)
	//now we test what we wrote to the file
	//by rewriting the stuff we just read to a txt file, only first proc
	if (rank == 0) {

		std::ofstream out_file("infos_2D_read_test.txt");

		out_file << std::setprecision(16) << std::scientific;

		out_file <<
			infos.last_computed_order_N << ' ' <<
			infos.a << ' ' << infos.b << ' ' <<
			infos.last_computed_order_M << ' ' <<
			infos.c << ' ' << infos.d << '\n';

		out_file << "********************\n****  infos 2D  ****\n********************\n";

		unsigned int absolute_index;

		for (int i = 0; i < infos.last_computed_order_N; i++) {
			for (int j = 0; j < infos.last_computed_order_M; j++) {

				absolute_index = i * infos.last_computed_order_M + j;
				out_file << coefficients[absolute_index] << ' ';
			}

			out_file << '\n';
		}

		out_file.close();
	}*/
	/********************* END OF TESTING ********************/

	if (rank == 0) {
		std::cout << " Done!\n";
	}

	delete[] counts;
	delete[] displs;

	return cheby_user_2d;

}
