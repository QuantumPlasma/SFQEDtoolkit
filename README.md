Hi!

If you want to see whether SFQEDtoolkit works or not you simply have to use
the makefile we provided (and that you downloaded together with the repository).
The instructions are quite simple and summarized in the following points.

If you instead want to simply start using then the only things you need are the 
src and the coefficients folders.

If something goes wrong or it is not clear please do not hesitate to reach out to us!

0. Before starting, make sure you possess the openMPI wrapper of the c++ compiler
mpic++. If you don't, then please install openMPI and add to the PATH variable its
folder location. The fact is that the test applications we are about to compile
are meant to work in parallel. Moreover, for the compilation of the test files it
is required that you export the location of the toolkit to the enironment variable
SFQED_TOOLKIT. Therefore open the .bash_profile file in your home directory and add
this line

export SFQED_TOOLKIT=location/of/SFQEDtoolkit

where location/of/SFQEDtoolkit corresponds to the folder where the SFQEDtoolkit's src
directory is placed

After that, source your .bash_profile with

source .bash_profile 

1. **IMPORTANT:**
**this step is necessary only if you want to build the "accuracytester" target of point 2,**
**which is the only one needing the GSL library. The other two targets do not need it**
**and neither does the actual SFQEDtoolkit (obviously). Skip this step if you want**
**to reproduce the synchrophoton emission or pair production spectra, or if you want**
**to use test the MonteCarlo pusher.**

You need to install the GNU Scientific Library (GSL) and export its folder
location to the environment variable GSL_INSTALLATION_DIR. The library can be found
at https://www.gnu.org/software/gsl/. After the installation go to your home folder and
open the file

.bash_profile

Add the following line

export GSL_INSTALLATION_DIR=/path/to/the/GSL

then save the file and source it.

2. Go to the SFQEDtoolkit repository and then you can run either one of the following
three targets

a) make pushertester: to compile the Monte-Carlo pusher for the comparison of the LCFA and BLCFA models. 
    Please, look inside the main located
    at ./Test_Modules/Mains/Pusher_Main/SFQEDtoolkit_tester.cpp to see what you are simulating.
    The modified pusher is found at ./Test_Modules/Pusher_Modules/RelativisticSolver.cpp.
    Read it and change the probability approach if you like.
    After it is compiled you can run it with

    mpirun -n ${PROCS} -x OMP_NUM_THREADS=${THREADS} -x I_MPI_PIN_DOMAIN=omp pusher_tester

    where PROCS and THREADS are environment variables representing the number of cores and the
    number of threads respectively.
    You can plot the results of this execution by employing the python plotting script inside 
    ./python_LCFA-BLCFA_comparison_plotter.
    

b) make accuracytester: compile the executable that compares the actual energies with those
    returned by the toolkit. Its main is located at ./Test_Modules/Mains/Accuracy_Main/SFQEDtoolkit_accuracy_tester.cpp.
    After it is compiled you can run it with

    mpirun -n ${PROCS} -x OMP_NUM_THREADS=${THREADS} -x I_MPI_PIN_DOMAIN=omp accuracy_tester

    where PROCS and THREADS are environment variables representing the number of cores and the
    number of threads respectively.
    You can plot the results of this execution by employing the python plotting script inside 
    ./python_accuracy_plotter_for_pairelectron_and_synchrophoton.

c) make spectrumtester: compile the executable that uses the energies returned by the toolkit
    to create the emission spectrum (main located at ./Test_Modules/Mains/Accuracy_Main/SFQEDtoolkit_spectrum_sampler.cpp.)
    After it is compiled you can run it with

    mpirun -n ${PROCS} -x OMP_NUM_THREADS=${THREADS} -x I_MPI_PIN_DOMAIN=omp spectrum_tester

    where PROCS and THREADS are environment variables representing the number of cores and the
    number of threads respectively.
    You can plot the results of this execution by employing the python plotting script inside 
    ./python_photon_and_electron_nrgs_sampler.



