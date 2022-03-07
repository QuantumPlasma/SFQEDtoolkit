#include "SFQED_Processes.h"
#include "RelativisticSolver.h"
#include "SFQED_Admin_Constants.h"

#include <iostream>
#include <iomanip>
#include <omp.h>
#include <chrono>
#include <random>
#include <vector>

using namespace std;

// The space for the static variable is allocated only one time and
// this is used for the entirety of the program.
// Once this variable is declared, it exists till the program executes.
// So, the lifetime of a static variable is the lifetime of the program.
double rnd_generator(){
	static std::random_device rd;  // Will be used to obtain a seed for the random number engine
    static std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    static std::uniform_real_distribution<> dis(0.0, 1.0);

	return dis(gen);
}

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

	////////////////////////////////////////////////
	//simulation parameters
    ////////////////////////////////////////////////

	double wave_length = 0.8e-6; //[m]

	//extract the reference length
	double reference_length = wave_length / (2. * PI);//c / w_r; //c / 3.30768e+15; //

    //reference angular frequency
    double w_r = c / reference_length;// [1/s] //3.30768e+15;//

	double normalized_sim_time = 1.257e2;

	double time_step = 0.1;//0.1;

    //number of pusher iterations
    int iterations = normalized_sim_time / time_step;

    ///////////////////////////////
	//plane wave or fields creation
	///////////////////////////////
	
	PlaneWave3D wave;

    double a0 = 2 * 10.; // 10.1433; //
	double fwhm = 5e-15 * w_r; //[normalized]
	double peak_position = 30. * 2. * PI;
	//direction propagation of the wave
	Vector3D<double> k(0., 0., 1.);
	Vector3D<double> A(1., 0., 0.);

    //usual plane wave constructor
	// wave = PlaneWave3D(k, A);
	// wave.seta0(100.0);
	// wave.setDelta(1.0);

    //the wave created through this constructor is a gaussian plane wave
	wave = PlaneWave3D(k, A, a0, fwhm, peak_position);

    /////////////////////////////////
    // workload splitting section
    /////////////////////////////////

    //number of electrons
    unsigned int M = 1000000;
    unsigned int N_electrons = 100 * M;

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

    ///////////////////
	//particles creation
	///////////////////

    //position of the particle (we place the
    //particle where the fields are maximum)
    //[REMEMBER THE LENGTHS ARE NORMALIZED TO THE PLANE WAVE \lambda / (2 \pi),
    // and thus the normalized wave length of the wave is trivially 2 \pi]
	Vector3D<double> pos(0., 0., PI * 100, "Position");

	// Vector3D<double> pos(0., 0., 0., "Position");

	//we want the particle to be an electron
    double mass = 1.,
            charge = -1.;

	//particle reduced momentum (we want a particle 
    //propagating against the plane wave with an initial energy of 10 GeV)
    double particle_gamma = 2 * 19570.513;//9784.756875; // 19569.5; // 
    Vector3D<double> u(0., 0., -particle_gamma);

	// Vector3D<double> u = wave.extractEightParticleMomentumLinear(pos, 0.0, mass, charge);
    
	//this pointer will store the particles that every
    //process needs
	Particle3D *particle;	
	
	/////////////////////////
    //file creation
    /////////////////////////

    string procs_suffix = to_string(rank);

	//KEEP IN MIND THAT THE FILE WILL NOW CONTAIN
	//THE ENERGIES TO THE EMITTED PHOTONS OR TO THE
	//CREATED ELECTRONS (DEPENDING ON WHICH METHOD
	// YOU ARE TESTING) 
    string name = "./nrgs" + procs_suffix + ".txt";

    ofstream out_file(name);

    out_file << std::setprecision(16) << std::scientific;

    out_file << particle_gamma << "\n";

    /////////////////////////////////
	//SFQED processess initialization
	/////////////////////////////////
	
	SFQED_Processes procs_instance;
    if (rank == 0) {
		std::cout << "Initializing simulation object... " << std::flush;
	}

	string SFQEDtoolkit_location_var_name = "SFQED_TOOLKIT";
	string SFQEDtoolkit_location = GetEnv(SFQEDtoolkit_location_var_name);


	procs_instance.SFQED_set_reference_length(reference_length);

	//BEFORE WE WERE USING
	// "/home/smonte/extend/spectrum_tester/new_coefficients/"
	procs_instance.SFQED_init_PHTN_emission(SFQEDtoolkit_location + "/coefficients/");
	procs_instance.SFQED_init_PAIR_creation(SFQEDtoolkit_location + "/coefficients/");

    procs_instance.SFQED_set_time_step(time_step);
    
    ////////////////////////////////////////////////
	//actual simulation
    ////////////////////////////////////////////////

    if (rank == 0) {
		std::cout << "Done!\nEntering simulation Region. Recap:\n"
					<< "reference length = " << reference_length
                    << "\nSimulation time = " << normalized_sim_time
                    << "\ntime_step = " << time_step << "\n" 
                    << std::flush;
	}

	//process parallelized loop
    for(int index = 0; index < counts[rank]; index++){
		// out_file << "#" << index << '\n';
		particle = new Particle3D(pos, u, mass, charge);
        particle->init_optical_depth(-log(1. - rnd_generator()));

		///////////////////////////////////////////////////////////
		// UNCOMMENT THIS TO TEST SYNCHROPHOTON EMISSION IN LCFA //
		///////////////////////////////////////////////////////////
        // RelativisticRungeKuttaSolver_LCFA(out_file, *particle, wave, time_step, iterations, procs_instance, &rnd_generator);
		///////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////
		// UNCOMMENT THIS TO TEST SYNCHROPHOTON EMISSION IN BLCFA //
		////////////////////////////////////////////////////////////
        // RelativisticRungeKuttaSolver_BLCFA(out_file, *particle, wave, time_step, iterations, procs_instance, &rnd_generator);
		///////////////////////////////////////////////////////////

		///////////////////////////////////////////////////////////////////
		// UNCOMMENT THIS TO TEST BREIT-WHEELER PAIR PRODUCTION IN BLCFA //
		///////////////////////////////////////////////////////////////////
		RelativisticRungeKuttaSolverPhotons(out_file, *particle, wave, time_step, iterations, procs_instance, &rnd_generator);
		///////////////////////////////////////////////////////////////////

		delete particle;
    }

	out_file.close();

	MPI_Barrier(newcomm);
	
    if (rank == 0) {
		std::cout << "Simulation Done!\n" << std::flush;
		std::cout << "files closed and array of particles deleted!\n" << std::flush;
	}

	procs_instance.SFQED_finalize_PHTN_emission();
	procs_instance.SFQED_finalize_PAIR_creation();

	MPI_Barrier(newcomm);
	
    if (rank == 0) {
		std::cout << "object killed successfully!\n" << std::flush;
	}

	////////////////////////////////////////////////

	MPI_Finalize();

	// std::cout << "exiting program!\n" << std::flush;
	
}