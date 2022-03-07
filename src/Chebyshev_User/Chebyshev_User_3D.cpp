#include "Chebyshev_User_3D.h"

#include <omp.h>
#include <iostream>

Chebyshev_User_3D::Chebyshev_User_3D(unsigned int order_N, double a, double b,
                                        unsigned int order_M, double c, double d,
                                        unsigned int order_L, double e, double f): Chebyshev_User_3D(){
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

    this->last_computed_order_L = order_L;
    this->evaluation_order_L = order_L;
    this->e = e;
    this->f = f;
    this->domain_half_range_L = (f-e)*0.5;
    this->domain_middle_point_M = (f+e)*0.5;
}

double Chebyshev_User_3D::evaluate(double x, double y, double z) {
	//variable conversion to -1 <= X <= +1 and -1 <= Y >= +1 and -1 <= Z <= +1

	double X = (x - domain_middle_point_N) / domain_half_range_N;

	double Y = (y - domain_middle_point_M) / domain_half_range_M;

	double Z = (z - domain_middle_point_L) / domain_half_range_L;

	//Clenshaw recurrence algorithm (easy peasy)
	//this stores the actual value that will be used in the algorithm
	double x_alg = X * 2.;
	double y_alg = Y * 2.;
	double z_alg = Z * 2.;

	//define the first 2 coefficients of the Clenshaw series to be 0.
	//(inside the algorithm loop we will update these two double variables,
	//whose names should be meaningful)
	double tmp, clensh_n_plus_1 = 0., clensh_n = 0.;

	unsigned int 
				//mtrx_order = last_computed_order_N * last_computed_order_L,
				mtrx_order = evaluation_order_N * evaluation_order_L,
				plane_number = last_computed_order_N*last_computed_order_M;
	int i, j, k;

	//Clenshaw coefficients
	double *d_ki = new double[mtrx_order],
		*d_k = new double[last_computed_order_L];

    /**********************************/
    /* WORKLOAD SPLITTING 1st SESSION */
    /**********************************/
    //notice that we are just giving the possibility to evaluate
    //the approx by parallelizing the computation.
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
	fcounts[0] = mtrx_order / nProcs + (0 < (mtrx_order % nProcs));

	//the following for section cannot be parallelized (omp),
	//as each iteration employs the result of the one
	//supposedly ended right before
	for (int i = 1; i < nProcs; i++) {
		fcounts[i] = mtrx_order / nProcs + (i < (mtrx_order % nProcs));
		fdispls[i] = fcounts[i - 1] + fdispls[i - 1];
	}
*/

	int loop_start = 0; //fdispls[rank];
	int loop_end = mtrx_order; //loop_start + fcounts[rank];

    //start computation

	//new coefficients determination loop
	#pragma omp parallel for default(shared) private(i,j,k,tmp) firstprivate(clensh_n_plus_1, clensh_n)
	for (unsigned int aux = loop_start; aux < loop_end; aux++) {
		//contrary to the usual way of thinking k represents the row number of the d_ki clenshaw 3d matrix
		k = aux / evaluation_order_N;
		//while i represents the col number
		i = aux % evaluation_order_N;
		//unsigned int pippo;

		//first Clenshaw algorithm (sum over the ys)
		for (j = evaluation_order_M - 1; j >= 0; j--) {
			//pippo = k* plane_number + i * M + j;
			//std::cout << pippo << "(k=" << k << ",i=" << i << ",j=" << j << ")\n";
			tmp = clensh_n;
			clensh_n = last_coeffs[k * plane_number + i * last_computed_order_M + j] - clensh_n_plus_1 + y_alg * clensh_n;
			clensh_n_plus_1 = tmp;
		}

		d_ki[aux] = clensh_n - Y * clensh_n_plus_1;

		//Reset the first 2 Clenshaw coefficient at each outer iteration
		//this must be here in case omp is not available
		clensh_n_plus_1 = 0.;
		clensh_n = 0.;
	}

    //TRANSMIT COEFFICIENTS TO ALL PROCS
    //MPI_Allgatherv(MPI_IN_PLACE, fcounts[rank], MPI_DOUBLE, d_ki, fcounts, fdispls, MPI_DOUBLE, *evalcom);

    /**********************************/
    /* WORKLOAD SPLITTING 2nd SESSION */
    /**********************************/
    //notice that we are just giving the possibility to evaluate
    //the approx by parallelizing the computation.
    //If you do not want to parallelize this, just set the MPI_COMM_SELF in the evalcom

/*
    fdispls[0] = 0;
	fcounts[0] = evaluation_order_L / nProcs + (0 < (evaluation_order_L % nProcs));

	//the following for section cannot be parallelized (omp),
	//as each iteration employs the result of the one
	//supposedly ended right before
	for (int i = 1; i < nProcs; i++) {
		fcounts[i] = evaluation_order_L / nProcs + (i < (evaluation_order_L % nProcs));
		fdispls[i] = fcounts[i - 1] + fdispls[i - 1];
	}
*/

	loop_start = 0; //fdispls[rank];
	loop_end = evaluation_order_L; //loop_start + fcounts[rank];

	//new coefficients determination loop (that on y)
	#pragma omp parallel for default(shared) private(k, i, tmp) firstprivate(clensh_n_plus_1, clensh_n)
	for (k = loop_start; k < loop_end; k++) {

		//second Clenshaw algorithm (sum on the xs)
		for (i = evaluation_order_N - 1; i >= 0; i--) {
			tmp = clensh_n;
			clensh_n = d_ki[k * evaluation_order_N + i] - clensh_n_plus_1 + x_alg * clensh_n;
			clensh_n_plus_1 = tmp;
		}

		//assign new coefficient (remember to use Y)
		d_k[k] = clensh_n - X * clensh_n_plus_1;

		//Reset the first 2 Clenshaw coefficient at each outer iteration
		//this must be here in case omp is not available
		clensh_n_plus_1 = 0.;
		clensh_n = 0.;
	}

    //TRANSMIT COEFFICIENTS TO ALL PROCS
    //MPI_Allgatherv(MPI_IN_PLACE, fcounts[rank], MPI_DOUBLE, d_k, fcounts, fdispls, MPI_DOUBLE, *evalcom);

	//really important part: delete what you're not using
	delete[] d_ki;

	//third Clenshaw algorithm application (sum over zs)
	//ATTENTION: not parallelizable
	for (k = evaluation_order_L - 1; k >= 0; k--) {
		tmp = clensh_n;
		clensh_n = d_k[k] - clensh_n_plus_1 + z_alg * clensh_n;
		clensh_n_plus_1 = tmp;
	}

	//really important part: delete what you're not using
	delete[] d_k;
    // delete[] fcounts;
    // delete[] fdispls;

	//remember to use X
	return clensh_n - Z * clensh_n_plus_1;
}


Chebyshev_User_3D Chebyshev_User_3D::init_from_txt_file(std::ifstream& in_file) {

    unsigned int order_N, order_M, order_L;
    double a, b, c, d, e, f;
	
    in_file >> order_N >> a >> b
	        >> order_M >> c >> d
			>> order_L >> e >> f;

    Chebyshev_User_3D cheby_user_3d(order_N, a, b, order_M, c, d, order_L, e, f);

	unsigned int order = order_N * order_M * order_L;
	cheby_user_3d.last_coeffs = new double[order];

	for (unsigned int i = 0; i < order; i++) {
		
		in_file >> cheby_user_3d.last_coeffs[i];
	}

    return cheby_user_3d;
}


Chebyshev_User_3D Chebyshev_User_3D::init_from_bin_file(MPI_File& file, MPI_Comm& comm) {

	int rank;
	MPI_Comm_rank(comm, &rank);

	int nProcs;
	MPI_Comm_size(comm, &nProcs);

	if (rank == 0) {
		std::cout << "Reading orders and domain from file...";
	}

    unsigned int order_N, order_M, order_L;
    double a, b, c, d, e, f;

	//let's read the infos obj related informations we have available on the file.
	//All the procs will do so. THIS WILL BE VERY UNELEGANT...
	MPI_File_read(file, &(order_N), 1, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_read(file, &(a), 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_read(file, &(b), 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_read(file, &(order_M), 1, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_read(file, &(c), 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_read(file, &(d), 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_read(file, &(order_L), 1, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_read(file, &(e), 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_read(file, &(f), 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

    Chebyshev_User_3D cheby_user_3d(order_N, a, b, order_M, c, d, order_L, e, f);

	if (rank == 0) {
		std::cout << " Done!\nSplitting workload...";
	}

	int* counts = new int[nProcs];
	int* displs = new int[nProcs];

	unsigned int coeff_order = order_N * order_M * order_L;

	counts[0] = coeff_order / nProcs + (0 < (coeff_order % nProcs));
	displs[0] = 0;

	for (int i = 1; i < nProcs; i++) {
		counts[i] = (coeff_order / nProcs + (i < coeff_order % nProcs));
		displs[i] = (displs[i - 1] + counts[i - 1]);
	}

	int start = displs[rank];

	//we already read something so a pre offset has to be built accordingly
	int pre_offset = sizeof(unsigned int) * 3 + sizeof(double) * 6;

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
	MPI_Allgatherv(MPI_IN_PLACE, counts[rank], MPI_DOUBLE, coefficients, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

	//pack the coeffs inside the infos obj
	cheby_user_3d.last_coeffs = coefficients;

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

	return cheby_user_3d;

}