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

#include "SFQED_Processes.h"
#include "SFQED_User_Constants.h"

#include <iostream>



void SFQED_Processes::SFQED_set_reference_length(const double& ref_length){
    Lambda = ref_length;

    //reference angular frequency = c / Lambda
    omega_r = light_speed / Lambda;

    //normalized compton time
    norm_Compton_time = compton_time * omega_r;

    //2 \pi \lambda_C / \lambda
    twopiComptonDivLambda = twopiCompton / Lambda;
    //\lambda / 2 \pi \lambda_C
    LambdaDivtwopiCompton = Lambda / twopiCompton;

    //\lambda_C / \lambda
    ComptonDivLambda = compton_length / Lambda;
    // \lambda / \lambda_C 
    LambdaDivCompton = Lambda / compton_length;

    //(\alpha / \sqrt{3} \pi) * (\lambda / 2 \pi \lambda_C)
    //coef_rate_Wrad = coef_W_rad*LambdaDivtwopiCompton;
    //(\alpha / \sqrt{3} \pi) * (\lambda / \lambda_C)
    coef_rate_Wrad = coef_W_rad*LambdaDivCompton;
    
    //(\alpha / 3 \sqrt{3} \pi) * (\lambda / 2 \pi \lambda_C)
    // coef_rate_Wpair = coef_W_pair*LambdaDivtwopiCompton;
    //(\alpha / 3 \sqrt{3} \pi) * (\lambda / \lambda_C)
    coef_rate_Wpair = coef_W_pair*LambdaDivCompton;

}

void SFQED_Processes::SFQED_set_reference_angular_frequency(const double& ref_freq){

    omega_r = ref_freq;
    
    //reference angular frequency = c / Lambda
    Lambda = light_speed / omega_r;

    //normalized compton time
    norm_Compton_time = compton_time * omega_r;

    //2 \pi \lambda_C / \lambda
    twopiComptonDivLambda = twopiCompton / Lambda;
    //\lambda / 2 \pi \lambda_C
    LambdaDivtwopiCompton = Lambda / twopiCompton;

    //\lambda_C / \lambda
    ComptonDivLambda = compton_length / Lambda;
    // \lambda / \lambda_C 
    LambdaDivCompton = Lambda / compton_length;

    //(\alpha / \sqrt{3} \pi) * (\lambda / 2 \pi \lambda_C)
    //coef_rate_Wrad = coef_W_rad*LambdaDivtwopiCompton;
    //(\alpha / \sqrt{3} \pi) * (\lambda / \lambda_C)
    coef_rate_Wrad = coef_W_rad*LambdaDivCompton;
    
    //(\alpha / 3 \sqrt{3} \pi) * (\lambda / 2 \pi \lambda_C)
    // coef_rate_Wpair = coef_W_pair*LambdaDivtwopiCompton;
    //(\alpha / 3 \sqrt{3} \pi) * (\lambda / \lambda_C)
    coef_rate_Wpair = coef_W_pair*LambdaDivCompton;

}

void SFQED_Processes::SFQED_set_all_to_one(){

    omega_r = 1;
    
    //reference angular frequency = c / Lambda
    Lambda = 1;

    //normalized compton time
    norm_Compton_time = 1;

    //2 \pi \lambda_C / \lambda
    twopiComptonDivLambda = 1;
    //\lambda / 2 \pi \lambda_C
    LambdaDivtwopiCompton = 1;

    //\lambda_C / \lambda
    ComptonDivLambda = 1;
    // \lambda / \lambda_C 
    LambdaDivCompton = 1;

    //(\alpha / \sqrt{3} \pi) * (\lambda / 2 \pi \lambda_C)
    //coef_rate_Wrad = coef_W_rad*LambdaDivtwopiCompton;
    //(\alpha / \sqrt{3} \pi) * (\lambda / \lambda_C)
    coef_rate_Wrad = 1;
    
    //(\alpha / 3 \sqrt{3} \pi) * (\lambda / 2 \pi \lambda_C)
    // coef_rate_Wpair = coef_W_pair*LambdaDivtwopiCompton;
    //(\alpha / 3 \sqrt{3} \pi) * (\lambda / \lambda_C)
    coef_rate_Wpair = 1;

}