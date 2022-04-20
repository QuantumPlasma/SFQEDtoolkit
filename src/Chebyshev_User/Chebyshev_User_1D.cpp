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

#include "Chebyshev_User_1D.h"

#include <iostream>
#include <sstream>
#include <string>

Chebyshev_User_1D::Chebyshev_User_1D(unsigned int order, double a, double b){
    this->evaluation_order = order;
    this->a = a;
    this->b = b;
    this->domain_half_range = (b-a)*0.5;
    this->domain_middle_point = (b+a)*0.5;
}

// Chebyshev_User_1D Chebyshev_User_1D::init_from_txt_file(std::ifstream& input_file){
//     unsigned int order;
//     double a, b;

//     input_file >> order >> a >> b;

//     Chebyshev_User_1D cheby_user_1d(order, a, b);

// 	cheby_user_1d.last_coeffs = new double[order];

// 	for (unsigned int i = 0; i < order; i++) {
// 		input_file >> cheby_user_1d.last_coeffs[i];
// 	}

//     //by returning the obj in this way, you are requesting
//     //to COPY the cheby_user_1d object from the stack of the
//     //current function to the heap of the main program
//     //(or wherever this method is called).
//     return cheby_user_1d;

//     //PLEASE USE RVO OR NRVO
//     //The common idea of these two optimizations is to allow the compiler
//     //to use the memory space of this object T t, which is outside the function,
//     //to directly construct the object being initialized inside the function
//     //and that is returned from it. This effectively removes the need for copying intermediary objects.

//     //But for the RVO to be applied, the returned object has to be constructed on a return statement.
//     //Therefore this object does not have a name.
//     //The NRVO (Named-RVO) goes one step further: it can remove the intermediary objects
//     //even if the returned object has a name and is therefore not constructed on the return statement.
// }

Chebyshev_User_1D Chebyshev_User_1D::init_from_txt_file(std::ifstream& input_file){
    unsigned int order[2];

	std::string buffer;
	std::getline(input_file, buffer);
	std::istringstream converted_buffer(buffer);

	while (converted_buffer.peek() == ':' || converted_buffer.peek() == ' ' || converted_buffer.peek() == ','){
        converted_buffer.ignore();
    }
	for (int index = 0; index < 2 && !converted_buffer.eof(); index++) {

        converted_buffer >> order[index];
		
		order[1] = order[index];

        while (converted_buffer.peek() == ':' || converted_buffer.peek() == ' ' || converted_buffer.peek() == ','){
            converted_buffer.ignore();
        }  
        
    }

    double a, b;

    input_file >> a >> b;

    Chebyshev_User_1D cheby_user_1d(order[1], a, b);

	cheby_user_1d.last_coeffs = new double[order[1]];

	for (unsigned int i = 0; i < order[1]; i++) {
		input_file >> cheby_user_1d.last_coeffs[i];
	}

	// std::cout << order[1] << '\n' << a << '\n' << b << '\n';
	// for(unsigned int i = 0; i < order[1]; i++) {
	// 	std::cout << cheby_user_1d.last_coeffs[i] << '\n';
	// }

    //by returning the obj in this way, you are requesting
    //to COPY the cheby_user_1d object from the stack of the
    //current function to the heap of the main program
    //(or wherever this method is called).
    return cheby_user_1d;

    //PLEASE USE RVO OR NRVO
    //The common idea of these two optimizations is to allow the compiler
    //to use the memory space of this object T t, which is outside the function,
    //to directly construct the object being initialized inside the function
    //and that is returned from it. This effectively removes the need for copying intermediary objects.

    //But for the RVO to be applied, the returned object has to be constructed on a return statement.
    //Therefore this object does not have a name.
    //The NRVO (Named-RVO) goes one step further: it can remove the intermediary objects
    //even if the returned object has a name and is therefore not constructed on the return statement.
}

Chebyshev_User_1D Chebyshev_User_1D::init_from_bin_file(MPI_File& file, MPI_Comm &communicator) {

    //we will employ all the procs available in communicator to read this binary file
	int rank;
	MPI_Comm_rank(communicator, &rank);

	int nProcs;
	MPI_Comm_size(communicator, &nProcs);

	if (rank == 0) {
		std::cout << "Reading order and domain from file...";
	}

    unsigned int order;
    double a, b;

	//let's read the infos obj related informations we have available on the file.
	//All the procs will do so. THIS WILL BE VERY UNELEGANT (MAYBE)...
	MPI_File_read(file, &(order), 1, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_read(file, &(a), 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_read(file, &(b), 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

    Chebyshev_User_1D cheby_user_1d(order, a, b);

	if (rank == 0) {
		std::cout << " Done!\nSplitting workload...";
	}

	int* counts = new int[nProcs];
	int* displs = new int[nProcs];

	counts[0] = order / nProcs + (0 < (order % nProcs));
	displs[0] = 0;

	for (int i = 1; i < nProcs; i++) {
		counts[i] = (order / nProcs + (i < order% nProcs));
		displs[i] = (displs[i - 1] + counts[i - 1]);
	}

	int start = displs[rank];

	//we already read something so a pre offset has to be built accordingly
	int pre_offset = sizeof(unsigned int) * 1 + sizeof(double) * 2;

	//start should contain the distance in (double), thus we adjust it to be a proper offset
	MPI_Offset offset = sizeof(double) * start + pre_offset;

	MPI_File_set_view(file, offset, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

	//allocate the memory for our coefficients
	double* coefficients = new double[order];

	//in case you need to reset the file pointer
	//MPI_File_seek(file, offset, MPI_SEEK_SET);

	if (rank == 0) {
		std::cout << " Done!\nReading coefficients from file...";
	}

	//read
	MPI_File_read(file, &(coefficients[start]), counts[rank], MPI_DOUBLE, MPI_STATUS_IGNORE);

	//gather what a proc read to all other procs
	MPI_Allgatherv(MPI_IN_PLACE, counts[rank], MPI_DOUBLE, coefficients, counts, displs, MPI_DOUBLE, communicator);

	/********************* TESTING PART *******************
	//DO NOT DELETE THIS PART (it may come in hand)
	//now we test what we wrote to the file
	//by rewriting the stuff we just read to a txt file, only first proc
	if (rank == 0) {

		std::ofstream out_file("infos_read_test.txt");

		out_file << std::setprecision(16) << std::scientific;

		out_file <<
			infos.last_computed_order_N << ' ' <<
			infos.a << ' ' << infos.b << '\n';

		out_file << "********************\n****  infos 1D  ****\n********************\n";

		for (int i = 0; i < number_of_coeffs; i++) {
			out_file << coefficients[i] << '\n';
		}

		out_file.close();
	}*/
	/********************* END OF TESTING ********************/

	if (rank == 0) {
		std::cout << " Done!\n";
	}

	delete[] counts;
	delete[] displs;

    cheby_user_1d.last_coeffs = coefficients;

	return cheby_user_1d;

}