#include "SFQED_Processes.h"
#include "SFQED_Functions.h"

#include <iostream>
#include <omp.h>
#include <cmath>
// #include <chrono>
#include <random>
#include <iomanip> 

using namespace std;

//environment variable reader.
//taken from https://stackoverflow.com/questions/5866134/how-to-read-linux-environment-variables-in-c/5866174
std::string GetEnv( const std::string & var ) {
     const char * val = std::getenv( var.c_str() );
     if ( val == nullptr ) { // invalid to assign nullptr to std::string
         return "";
     }
     else {
         return val;
     }
}

//MAIN FOR ACCURACY TESTING

void write_vec_on_file(MPI_File& file, int& write_count, int& start, double *vec) {

	int pre_offset = sizeof(unsigned int) * 2 + sizeof(double) * 4;

	// //start should contain the distance in (double), thus we adjust it to be a proper offset
	MPI_Offset offset = sizeof(double) * start + pre_offset;

	MPI_File_set_view(file, offset, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

	MPI_File_write(file, vec, write_count, MPI_DOUBLE, MPI_STATUS_IGNORE);
}

void test_photon_emission_rate_accuracy(SFQED_Processes &procs_instance, int &loop_start, int &loop_end, double *diffs,
									double &chi_min, double &chi_max, const int &N_points_chi, MPI_Comm& comm){
	
    double dchi = (chi_max - chi_min) / N_points_chi;

	int i_x_2, index_x_2_plus_1, index;
	double tmp_diff;

	#pragma omp parallel for default(shared) private(i_x_2, index_x_2_plus_1, index, tmp_diff)
	for (int i = loop_start; i < loop_end; i++) {

		index = i - loop_start;

		//new
		i_x_2 = index * 2;
		index_x_2_plus_1 = i_x_2 + 1;

		diffs[i_x_2] = chi_min + i * dchi;

		diffs[index_x_2_plus_1] = tildeWrad(diffs[i_x_2], NULL);

		tmp_diff = procs_instance.SFQED_PHTN_emission_rate(1., diffs[i_x_2]);
		tmp_diff = (diffs[i_x_2] != 0.) ? tmp_diff / diffs[i_x_2] : tmp_diff;

		//NaN check
		if(tmp_diff != tmp_diff){
			#pragma omp critical
			cout << "error at chi=" << diffs[i_x_2] << ". Aborting!\n";
			MPI_Abort(comm, 666);
		}

		//relative difference
		if(diffs[index_x_2_plus_1]  != 0.){
			diffs[index_x_2_plus_1] = std::abs((tmp_diff - diffs[index_x_2_plus_1]) / diffs[index_x_2_plus_1]);
		}
		else if(tmp_diff != 0.) {
			diffs[index_x_2_plus_1] = std::abs((tmp_diff - diffs[index_x_2_plus_1]) / tmp_diff);
		}
		//we arrive in this last case only if both the prevs are null
		else {
			diffs[index_x_2_plus_1] = 0.;
		}

	}
	
}

void test_photon_emission_accuracy(SFQED_Processes &procs_instance, int &loop_start, int &loop_end, double *diffs,
									double &chi_min, double &chi_max, const int &N_points_chi,
									double &r_min, double &r_max, const int &N_points_r, MPI_Comm& comm){
	
    double dchi = (chi_max - chi_min) / N_points_chi;
	double dr = (r_max - r_min) / N_points_r;
    int i, j, index;
    double chi, rnd, tmp_val;
	
	//simulation loop (remember that nested parallelism is turned off by default and you have a parallel region inside SFQED_photon_momentum)
	#pragma omp parallel for default(shared) private(i, j, chi, rnd, index, tmp_val)
	for(int iter = loop_start; iter < loop_end; iter++){
		
        //remember this is a 2D collapsed loop
        i = iter / N_points_r;
		j = iter % N_points_r;

        chi = chi_min + i * dchi;
		rnd = r_min + j * dr;

        index = iter - loop_start;

		diffs[index] = photon_nrg(chi, rnd);
        
		//we use 1 as gamma value to normalize the results
        tmp_val = procs_instance.SFQED_LCFA_emitted_photon_energy(1., chi, rnd);

		//NaN check
		if(tmp_val != tmp_val){
			#pragma omp critical
			cout << "error at chi=" << chi << ", rnd=" << rnd << ". Aborting!\n";
			MPI_Abort(comm, 666);
		}

		//relative difference
		if(diffs[index] != 0.) {
			diffs[index] = std::abs((tmp_val - diffs[index]) / diffs[index]);
		}
		else if(tmp_val != 0.){
			diffs[index] = std::abs((tmp_val - diffs[index]) / tmp_val);
		}
		//we arrive in this last case only if both the prevs are zero
		else {
			diffs[index] = 0.;
		}
	}
}

void test_pair_creation_accuracy(SFQED_Processes &procs_instance, int &loop_start, int &loop_end, double *diffs,
									double &chi_min, double &chi_max, const int &N_points_chi,
									double &r_min, double &r_max, const int &N_points_r, MPI_Comm& comm){
	
    double dchi = (chi_max - chi_min) / N_points_chi;
	double dr = (r_max - r_min) / N_points_r;
    int i, j, index;
    double chi, rnd, tmp_val;
	
	//simulation loop (remember that nested parallelism is turned off by default and you have a parallel region inside SFQED_photon_momentum)
	#pragma omp parallel for default(shared) private(i, j, chi, rnd, index, tmp_val)
	for(int iter = loop_start; iter < loop_end; iter++){
		
        //remember this is a 2D collapsed loop
        i = iter / N_points_r;
		j = iter % N_points_r;

        chi = chi_min + i * dchi;
		rnd = r_min + j * dr;

        index = iter - loop_start;

		diffs[index] = electron_nrg(chi, rnd);
		// diffs[index] = electron_nrg_small_k(chi, rnd);
        
		//we use 1 as gamma value to normalize the results
        tmp_val = procs_instance.SFQED_PAIR_created_electron_energy(1., chi, rnd);

		//NaN check
		if(tmp_val != tmp_val){
			#pragma omp critical
			cout << "error at chi=" << chi << ", rnd=" << rnd << ". Aborting!\n";
			MPI_Abort(comm, 666);
		}

		//relative difference
		if(diffs[index] != 0.) {
			diffs[index] = std::abs((tmp_val - diffs[index]) / diffs[index]);
		}
		else if(tmp_val != 0.){
			diffs[index] = std::abs((tmp_val - diffs[index]) / tmp_val);
		}
		//we arrive in this last case only if both the prevs are zero
		else {
			diffs[index] = 0.;
		}
	}
}


//accuracy tester
int main(int argc, char** argv) {

	int provided_thread_support;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided_thread_support);

	if (provided_thread_support < MPI_THREAD_FUNNELED) {
		std::fprintf(stderr, "Hybrid parallelization requires thread compliant MPI library.\n");
		std::exit(EXIT_FAILURE);
	}

	//create new communicator (we duplicate the old one)
	MPI_Comm newcomm;
	MPI_Comm_dup(MPI_COMM_WORLD, &newcomm);

	int rank;
	MPI_Comm_rank(newcomm, &rank);
	int size;
	MPI_Comm_size(newcomm, &size);

	if (rank == 0) {
		cout << "nProcs = " << size << endl;
		cout << "nThreads = " << omp_get_max_threads() << endl;
	}
    
	SFQED_Processes procs_instance;

    if (rank == 0) {
		std::cout << "Initializing coefficients... " << std::flush;
	}

	string SFQEDtoolkit_location_var_name = "SFQED_TOOLKIT";
	string SFQEDtoolkit_location = GetEnv(SFQEDtoolkit_location_var_name);

	//BEFORE WE WERE USING
	// "/home/smonte/extend/spectrum_tester/new_coefficients/"
	procs_instance.SFQED_init_PHTN_emission(SFQEDtoolkit_location + "/coefficients/");
	procs_instance.SFQED_init_PAIR_creation(SFQEDtoolkit_location + "/coefficients/");

	procs_instance.SFQED_set_all_to_one();
	
	//TOT points
	unsigned int N_points_chi = 500;
	unsigned int N_points_r = 500;

	//////////////////////////////////////////

    unsigned int N_tot = N_points_chi * N_points_r;
	
	// unsigned int N_tot = 100000;

	if (rank == 0) {
		std::cout << " Done!\n" << std::flush;
	}

	///////////////////////////////////////////////
	/*workload splitting section (rank and size have already been determined above)
	* every 'count' or 'displacement' array encountered below
	* in the MPI collective functions will refer to some
	* variable defined within this section
	*/

	//vectors used to store the counts and displacements informations
	//used during process' communication
	int* counts = new int[size];
	int* displs = new int[size];

	//////////////////////////////////////////////////
	//FOR RATE ACCURACY TESTER

	displs[0] = 0;
	counts[0] = ((N_points_chi) / size + (0 < ((N_points_chi) % size))) * 2;

	//the following for section cannot be parallelized,
	//as each iteration employs the result of the one
	//supposedly ended right before
	for (int i = 1; i < size; i++) {
		counts[i] = ((N_points_chi) / size + (i < ((N_points_chi) % size))) * 2;
		displs[i] = counts[i - 1] + displs[i - 1];
	}

	//range on function evals loop
	int loop_start = displs[rank] / 2;
	int loop_end = loop_start + (counts[rank] / 2);

	///////////////////////////////////////////////////
	//FOR ENERGY ACCURACY TESTER

	// //set the first displacement to 0 and manually set the first entry of the function counts
	// displs[0] = 0;
	// counts[0] = (N_tot) / size + (0 < ((N_tot) % size));

	// //the following for section cannot be parallelized,
	// //as each iteration employs the result of the one
	// //supposedly ended right before
	// for (int i = 1; i < size; i++) {
	// 	counts[i] = N_tot / size + (i < (N_tot % size));
	// 	displs[i] = counts[i - 1] + displs[i - 1];
	// }

	// //range on function evals loop
	// int loop_start = displs[rank];
	// int loop_end = loop_start + counts[rank];

	//////////////////////////////////////////////////////////////////

	//this pointer will store the energy values
	double *diffs_high = new double[counts[rank]];

	/////////////////////
	// main variables
	/////////////////////

	// double chi_min = 0.0, chi_max = 2.0;

	// double chi_min = 0.01, chi_max = 0.3;
	// double chi_min = 0.3, chi_max = 2.0;

	// double chi_min = 2.0, chi_max = 20.0;
	double chi_min = 20.0, chi_max = 80.0;
	// double chi_min = 80.0, chi_max = 600.0;
	// double chi_min = 600.0, chi_max = 2000.0;

	double r_min = 0.0, r_max = 1.0;
	// double r_min = 0.5, r_max = 1.0;


	if (rank == 0) {
		std::cout << "Entering testing Region\n" << std::flush;
	}

	// std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();


	////////////////////////////////////////////////////////TEST RATE ACCURACY
	test_photon_emission_rate_accuracy(procs_instance, loop_start, loop_end, diffs_high,
										chi_min, chi_max, N_points_chi, newcomm);

	////////////////////////////////////////////////////////TEST ENERGY ACCURACY
	// test_photon_emission_accuracy(procs_instance, loop_start, loop_end, diffs_high,
    //                     chi_min, chi_max, N_points_chi,
	// 					r_min, r_max, N_points_r, newcomm);

	// test_pair_creation_accuracy(procs_instance, loop_start, loop_end, diffs_high,
    //                     chi_min, chi_max, N_points_chi,
	// 					r_min, r_max, N_points_r, newcomm);

	

	//write on file section
	MPI_File file;
	MPI_File_open(newcomm, "transition.dat", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
	MPI_File_set_errhandler(file, MPI_ERRORS_ARE_FATAL);

	//////////////////////////////////////////////////////////////////
	//WRITE ENERGY ACCURACY ON FILE

	// if (rank == 0) {
	// 	std::cout << "Writing photon energies on file... " << std::flush;
	// }

    // if (rank == 0) {
	// 	//notice this thing could have been done by using a MPI type
	// 	//but this should ensure the file produced to have the smallest
	// 	//size
	// 	MPI_File_write(file, &N_points_chi, 1, MPI_INT, MPI_STATUS_IGNORE);

	// 	double* edges_tmp = new double[2] {chi_min, chi_max};
	// 	MPI_File_write(file, edges_tmp, 2, MPI_DOUBLE, MPI_STATUS_IGNORE);

	// 	delete[] edges_tmp;

	// 	MPI_File_write(file, &N_points_r, 1, MPI_INT, MPI_STATUS_IGNORE);

	// 	edges_tmp = new double[2] {r_min, r_max};
	// 	MPI_File_write(file, edges_tmp, 2, MPI_DOUBLE, MPI_STATUS_IGNORE);

	// 	delete[] edges_tmp;
	// }

	// write_vec_on_file(file, counts[rank], displs[rank], diffs_high);

	// MPI_File_close(&file);

	// if (rank == 0) {
	// 	std::cout << " Done!\n" << std::flush;
	// }

	//////////////////////////////////////////////////////////////////
	//WRITE RATE ACCURACY ON FILE

	//first let's write the infos' obj main feature first
	if (rank == 0) {
		//notice this thing could have been done by using a MPI type
		//but this should ensure the file produced to have the smallest
		//size
		MPI_File_write(file, &(N_points_chi), 1, MPI_INT, MPI_STATUS_IGNORE);
	}

	//we already write something so a pre offset has to be built accordingly
	int pre_offset = sizeof(unsigned int) * 1;

	//start should contain the distance in (double), thus we adjust it to be a proper offset
	MPI_Offset offset = sizeof(double) * displs[rank] + pre_offset;

	//change the file view
	MPI_File_set_view(file, offset, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

	if (rank == 0) {
		std::cout << "Writing function evaluations on file...";
	}

	//previous writing method
	MPI_File_write(file, diffs_high, counts[rank], MPI_DOUBLE, MPI_STATUS_IGNORE);

	MPI_File_close(&file);

	if (rank == 0) {
		std::cout << " Done!\n";
	}

	//////////////////////////////////////////////////////////////////

	delete[] diffs_high;

	delete[] counts;
	delete[] displs;

	// test(fake_instance);

	procs_instance.SFQED_finalize_PHTN_emission();
	procs_instance.SFQED_finalize_PAIR_creation();

	////////////////////////////////////////////////

	MPI_Finalize();
}
