#include "SFQED_Processes.h"
// #include "SFQED_Functions.h"

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

//****************************************************************************************
//MAIN FOR SPECTRA SAMPLING

void sample_photon_spectrum(SFQED_Processes &procs_instance,
                                int &loop_start,
                                int &loop_end,
                                double *e_nrgs,
                                double& chi,
                                double& gamma,
								double dr){


        //variables needed in the main loop                                                                                                                                         
        double rnd, tempo;                                                                                                                                                                    

        //simulation loop (remember that nested parallelism is turned off by default and you have a parallel region inside SFQED_LCFA_emitted_photon_energy)                        
        #pragma omp parallel for private(rnd, tempo)
        for(int i = loop_start; i < loop_end; i++){
                //generate random number                                                                                                                                            
                // rnd = dis(gen);

				rnd = i * dr;

				tempo = procs_instance.SFQED_LCFA_emitted_photon_energy(gamma, chi, rnd);

				e_nrgs[i - loop_start] = tempo;
				
        }
}

void sample_pair_spectrum(SFQED_Processes &procs_instance,
                                int &loop_start,
                                int &loop_end,
                                double *e_nrgs,
                                double& chi,
                                double& gamma,
								double dr){


        //variables needed in the main loop                                                                                                                                         
        double rnd, tempo;                                                                                                                                                                    

        //simulation loop (remember that nested parallelism is turned off by default and you have a parallel region inside SFQED_LCFA_emitted_photon_energy)                        
        #pragma omp parallel for private(rnd, tempo)
        for(int i = loop_start; i < loop_end; i++){
                //generate random number                                                                                                                                            
                // rnd = dis(gen);

				rnd = i * dr;

				tempo = procs_instance.SFQED_PAIR_created_electron_energy(gamma, chi, rnd);

				e_nrgs[i - loop_start] = tempo;
				
        }
}


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

	/////////////////////////////////
	//SFQED processess initialization
	/////////////////////////////////
	
	SFQED_Processes procs_instance;

    if (rank == 0) {
		std::cout << "Initializing simulation... " << std::flush;
	}

	string SFQEDtoolkit_location_var_name = "SFQED_TOOLKIT";
	string SFQEDtoolkit_location = GetEnv(SFQEDtoolkit_location_var_name);

	//BEFORE WE WERE USING
	// "/home/smonte/extend/spectrum_tester/new_coefficients/"
	procs_instance.SFQED_init_PHTN_emission(SFQEDtoolkit_location + "/coefficients/");
	procs_instance.SFQED_init_PAIR_creation(SFQEDtoolkit_location + "/coefficients/");

    if (rank == 0) {
		std::cout << " Done!\n" << std::flush;
	}
	
    ofstream out_electrons("./electrons.txt");

    out_electrons << std::setprecision(16) << std::scientific;

	double chi_min = 0.0, chi_max = 2.0;

	// double chi_min = 0.01, chi_max = 0.3;
	// double chi_min = 0.3, chi_max = 2.0;

	// double chi_min = 2.0, chi_max = 20.0;
	// double chi_min = 20.0, chi_max = 80.0;
	// double chi_min = 80.0, chi_max = 600.0;
	// double chi_min = 600.0, chi_max = 2000.0;

	double r_min = 0.0, r_max = 1.0;
	// double r_min = 0.5, r_max = 1.0;

	//this should represent the quantum parameter of your
	//particles, change this accordingly
	double chi = 0.5 * (chi_max + chi_min);
	
	//this is the energy of the particle
	double nrg_ref = 1;
    
	//number of points
	unsigned int M = 1000000;
    unsigned int N_electrons = 100000; //10 * M;	

	////////////////////////////////////////////////
	//simulation simulation region

	/*workload splitting section (rank and size have already been determined above)
	* every 'count' or 'displacement' array encountered below
	* in the MPI collective functions will refer to some
	* variable defined within this section
	*/

	//vectors used to store the counts and displacements informations
	//used during process' communication
	int* counts = new int[size];
	int* displs = new int[size];

	//set the first displacement to 0 and manually set the first entry of the function counts
	displs[0] = 0;
	counts[0] = (N_electrons) / size + (0 < ((N_electrons) % size));

	//the following for section cannot be parallelized,
	//as each iteration employs the result of the one
	//supposedly ended right before
	for (int i = 1; i < size; i++) {
		counts[i] = N_electrons / size + (i < (N_electrons % size));
		displs[i] = counts[i - 1] + displs[i - 1];
	}

	//range on function evals loop
	int loop_start = displs[rank];
	int loop_end = loop_start + counts[rank];

	//this pointer will store the energy values
	double *e_nrgs = new double[counts[rank]];
    double *p_nrgs = new double[counts[rank]];

	if (rank == 0) {
		std::cout << "Entering simulation Region\n" << std::flush;
	}

	// std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	//spectra

	///////////////////////////////////////////////////////////////////////////
	// UNCOMMENT THIS TO SAMPLE THE ENTIRETY OF THE PAIR'S ELECTRON SPECTRUM //
	///////////////////////////////////////////////////////////////////////////
	sample_pair_spectrum(procs_instance, loop_start, loop_end, e_nrgs, chi, nrg_ref, 1. / N_electrons);

	///////////////////////////////////////////////////////////////////////////
	// UNCOMMENT THIS TO SAMPLE THE ENTIRETY OF THE SYNCHROPHOTON SPECTRUM //
	///////////////////////////////////////////////////////////////////////////
	// sample_photon_spectrum(procs_instance, loop_start, loop_end, e_nrgs, chi, nrg_ref, 1. / N_electrons);


    // std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	// if(rank == 0){
	// 	std::cout << "pair loops executed in t = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl << std::flush;
	// }

	//////////////////////////////////////////////////
	
	double *all_e_nrgs;

    if(rank == 0){
		std::cout << "Transmitting nrgs..." << std::flush;
		//only the root process initializes the pointers
        all_e_nrgs = new double[N_electrons];
    }

    MPI_Gather(e_nrgs, counts[rank], MPI_DOUBLE, all_e_nrgs, counts[rank], MPI_DOUBLE, 0, newcomm);

	if (rank == 0) {
		
        std::cout << " Done!\nWriting photon energies on file... " << std::flush;

        //write on file section
        out_electrons << chi << '\n' 
                << gamma << '\n';
         

        for(int i = 0; i < N_electrons; i++){
            out_electrons << all_e_nrgs[i] << '\n';
        }

		std::cout << " Done!\n" << std::flush;

        delete[] all_e_nrgs;
	}

	delete[] e_nrgs;

	out_electrons.close();

	procs_instance.SFQED_finalize_PHTN_emission();
	procs_instance.SFQED_finalize_PAIR_creation();

	////////////////////////////////////////////////

	MPI_Finalize();
}
