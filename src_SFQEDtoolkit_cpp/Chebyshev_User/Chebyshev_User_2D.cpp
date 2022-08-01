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

// #include <omp.h>
#include <iostream>
#include <sstream>
#include <string>

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

        //second Clenshaw algorithm application (on x)
        //ATTENTION: not parallelizable
        for (int i = evaluation_order_N - 1; i >= 0; i--) {
            tmp = clensh_j;
            clensh_j = d_i[i] - clensh_j_plus_1 + x_alg * clensh_j;
            clensh_j_plus_1 = tmp;
        }

        delete[] d_i;

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

	return d_i;
}

//comments omitted (simplified version of evaluate)
double* Chebyshev_User_2D::evaluate_x(double x) const{

	double X = (x - domain_middle_point_N) / domain_half_range_N;

	double x_alg = X * 2.;

	double tmp, clensh_j_plus_1 = 0., clensh_j = 0.;

	//Clenshaw coefficients
	double* d_i = new double[evaluation_order_M];

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

	return d_i;
}


Chebyshev_User_2D Chebyshev_User_2D::init_from_txt_file(std::ifstream& in_file){
    unsigned int order_N[2], order_M[2];
	
	std::string buffer;
	std::getline(in_file, buffer);
	std::istringstream converted_buffer(buffer);

	while (converted_buffer.peek() == ':' || converted_buffer.peek() == ' ' || converted_buffer.peek() == ','){
        converted_buffer.ignore();
    }
	for (int index = 0; index < 2 && !converted_buffer.eof(); index++) {
        
		converted_buffer >> order_N[index];
		
		order_N[1] = order_N[index];

        while (converted_buffer.peek() == ':' || converted_buffer.peek() == ' ' || converted_buffer.peek() == ','){
            converted_buffer.ignore();
        }  
        
    }

    double a, b, c, d, tmp;

    in_file >> a >> b;

	//remove the endline character following b
	std::getline(in_file, buffer);

	std::string buffer1;
	std::getline(in_file, buffer1);
	std::istringstream converted_buffer1(buffer1);

	while (converted_buffer1.peek() == ':' || converted_buffer1.peek() == ' ' || converted_buffer1.peek() == ','){
        converted_buffer1.ignore();
    }
	for (int index = 0; index < 2 && !converted_buffer1.eof(); index++) {
        
		converted_buffer1 >> order_M[index];

		order_M[1] = order_M[index];

        while (converted_buffer1.peek() == ':' || converted_buffer1.peek() == ' ' || converted_buffer1.peek() == ','){
            converted_buffer1.ignore();
        }  
        
    }
			
	in_file >> c >> d;

    Chebyshev_User_2D cheby_user_2d(order_N[1], a, b, order_M[1], c, d);

	unsigned int order = order_N[1]*order_M[1];

	cheby_user_2d.last_coeffs = new double[order];

	for(unsigned int i = 0; i < order_N[1]; i++){
		for(unsigned int j = 0; j < order_M[0]; j++){
			in_file >> tmp;

			if(j < order_M[1]){
				cheby_user_2d.last_coeffs[i * order_M[1] + j] = tmp;
			}
		}
	}

    //for debug
	// std::cout << buffer << '\n';
	// std::cout << order_N[1] << '\n' << a << '\n' << b << '\n';
	// std::cout << order_M[1] << '\n' << c << '\n' << d << '\n';
	// for(unsigned int i = 0; i < order_N[1]*order_M[1]; i++) {
	// 	std::cout << cheby_user_2d.last_coeffs[i] << '\n';
	// }

    return cheby_user_2d;
    
}
