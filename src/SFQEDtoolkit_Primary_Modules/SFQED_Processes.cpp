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

    norm_Compton_time_2 = norm_Compton_time * norm_Compton_time;

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

    norm_Compton_time_2 = norm_Compton_time * norm_Compton_time;

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

void SFQED_Processes::SFQED_set_time_step(const double& dt){
    one_over_dt = 1. / dt;
    one_over_dt2 = one_over_dt * one_over_dt;
}

// double dot_prod(const double v0[3], const double v1[3]){
//     return v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2];
// }

//MIND THE GAUSS UNITS
//remember that, if the momentum of the fermion p_in is given in terms of "m_e c", than
//gamma = sqrt(1 + |p_in|^2)
// this function should be used with p_in corresponding to the e- or e+ momentum at half-step
// calculate the \chi parameter. Momentum and fields MUST BE in units of
// "m_e c" and of "m_e c \omega / |e|" for EE and "m_e \omega / |e|" for BB, respectively
// double SFQED_Processes::SFQED_compute_quantum_param(double gamma, const double p_in[3], const double EE[3], const double BB[3]) {
    
//     double cross[3], ff[3];
    
//     // the calculation of \chi or \kappa in this form provides
//     // the most accurate result with also very good performance
//     cross[0] = p_in[1]*BB[2] - p_in[2]*BB[1];
//     cross[1] = p_in[2]*BB[0] - p_in[0]*BB[2];
//     cross[2] = p_in[0]*BB[1] - p_in[1]*BB[0];
    
//     ff[0] = gamma * EE[0] + cross[0];
//     ff[1] = gamma * EE[1] + cross[1];
//     ff[2] = gamma * EE[2] + cross[2];
    
//     double arg = dot_prod(p_in, EE);
//     arg = dot_prod(ff, ff) - (arg*arg);
    
//     //if ff is parallel to p_in and large, then round-off errors
//     //(or an invalid input gam factor) may render arg negative
//     //=> NaN is returned instead of zero
//     //we avoid this by checking the sign of arg
//     arg = arg < 0. ? 0. : arg;
    
//     //pay attention to the UNITS
//     //return twopiComptonDivLambda * sqrt(arg);
//     return ComptonDivLambda * sqrt(arg);
// }

// double SFQED_Processes::SFQED_compute_quantum_param(const double &gamma,
//                                                 const double &p_in_x,
//                                                 const double &p_in_y,
//                                                 const double &p_in_z,
//                                                 const double &EE_x,
//                                                 const double &EE_y,
//                                                 const double &EE_z,
//                                                 const double &BB_x,
//                                                 const double &BB_y,
//                                                 const double &BB_z) {
    
//     double cross[3], ff[3];
    
//     // the calculation of \chi or \kappa in this form provides
//     // the most accurate result with also very good performance
//     cross[0] = p_in_y*BB_z - p_in_z*BB_y;
//     cross[1] = p_in_z*BB_x - p_in_x*BB_z;
//     cross[2] = p_in_x*BB_y - p_in_y*BB_x;
    
//     ff[0] = gamma * EE_x + cross[0];
//     ff[1] = gamma * EE_y + cross[1];
//     ff[2] = gamma * EE_z + cross[2];
    
//     double arg = p_in_x * EE_x + p_in_y * EE_y + p_in_z * EE_z;
//     arg = dot_prod(ff, ff) - (arg*arg);
    
//     //if ff is parallel to p_in and large, then round-off errors
//     //(or an invalid input gam factor) may render arg negative
//     //=> NaN is returned instead of zero
//     //we avoid this by checking the sign of arg
//     arg = arg < 0. ? 0. : arg;
    
//     //pay attention to the UNITS
//     //return twopiComptonDivLambda * sqrt(arg);
//     return ComptonDivLambda * sqrt(arg);
// }

/*****************/
/*PHOTON EMISSION*/
/*****************/

int stupid_power_of_two(int exponent){
    int ret = 1;

    for(int i = 0; i < exponent; i++){
        ret *= 2;
    }

    return ret;
}

void SFQED_Processes::SFQED_init_PHTN_emission(std::string path_to_coeffs){

    init_phtn_mx_tables();

    int tmp_i;

    //photon emission rate
    std::string emission_name = path_to_coeffs + "full_rate/phtn_emission_rate_0-2.txt";
    std::ifstream in_file(emission_name.c_str());
    phtn_mx_rt_0_2 = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "full_rate/phtn_emission_rate_2-20.txt";
    in_file.open(emission_name.c_str());
    phtn_mx_rt_2_20 = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "full_rate/phtn_emission_rate_20-80.txt";
    in_file.open(emission_name.c_str());
    phtn_mx_rt_20_80 = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "full_rate/phtn_emission_rate_80-600.txt";
    in_file.open(emission_name.c_str());
    phtn_mx_rt_80_600 = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "full_rate/phtn_emission_rate_600-2000.txt";
    in_file.open(emission_name.c_str());
    phtn_mx_rt_600_2000 = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

     //prepare photon emission rate lookup table
    Chebyshev_User_1D *phtn_mx_rates[4] = {&phtn_mx_rt_80_600,
                                            &phtn_mx_rt_20_80,
                                            &phtn_mx_rt_2_20,
                                            &phtn_mx_rt_0_2};

    look_up_table_phtn_mx_rt[0] = &phtn_mx_rt_600_2000;                                            

    for(int i = 0; i < 4; i++){
        tmp_i = 1 << i;
        for(int j = 0; j < tmp_i; j++){
            look_up_table_phtn_mx_rt[tmp_i + j] = phtn_mx_rates[i];
        }
    }

    // look_up_table_phtn_mx_rt[0] = &phtn_mx_rt_600_2000;
    // look_up_table_phtn_mx_rt[1] = &phtn_mx_rt_80_600;
    // look_up_table_phtn_mx_rt[2] = &phtn_mx_rt_20_80;
    // look_up_table_phtn_mx_rt[3] = &phtn_mx_rt_2_20;
    // look_up_table_phtn_mx_rt[4] = &phtn_mx_rt_0_2;


    //photon partial emission rate g(w,\chi) = \int_0^w{Wrad(\chi,w_{integrated})dw_{integrated}}
    emission_name = path_to_coeffs + "prtl_rate/phtn_partial_rate_0-2.txt";
    in_file.open(emission_name.c_str());
    phtn_prtl_rt_0_2 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "prtl_rate/phtn_partial_rate_2-20.txt";
    in_file.open(emission_name.c_str());
    phtn_prtl_rt_2_20 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "prtl_rate/phtn_partial_rate_20-80.txt";
    in_file.open(emission_name.c_str());
    phtn_prtl_rt_20_80 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "prtl_rate/phtn_partial_rate_80-600.txt";
    in_file.open(emission_name.c_str());
    phtn_prtl_rt_80_600 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "prtl_rate/phtn_partial_rate_600-2000.txt";
    in_file.open(emission_name.c_str());
    phtn_prtl_rt_600_2000 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    //prepare photon emission rate lookup table
    Chebyshev_User_2D *phtn_prtl_rates[4] = {&phtn_prtl_rt_80_600,
                                            &phtn_prtl_rt_20_80,
                                            &phtn_prtl_rt_2_20,
                                            &phtn_prtl_rt_0_2};

    look_up_table_phtn_prtl_rt[0] = &phtn_prtl_rt_600_2000;                                            

    for(int i = 0; i < 4; i++){
        tmp_i = 1 << i;
        for(int j = 0; j < tmp_i; j++){
            look_up_table_phtn_prtl_rt[tmp_i + j] = phtn_prtl_rates[i];
        }
    }

    // look_up_table_phtn_prtl_rt[0] = &phtn_prtl_rt_600_2000;
    // look_up_table_phtn_prtl_rt[1] = &phtn_prtl_rt_80_600;
    // look_up_table_phtn_prtl_rt[2] = &phtn_prtl_rt_20_80;
    // look_up_table_phtn_prtl_rt[3] = &phtn_prtl_rt_2_20;
    // look_up_table_phtn_prtl_rt[4] = &phtn_prtl_rt_0_2;


    //photon energy in w
    emission_name = path_to_coeffs + "w/phtn_w_nrg_0-2.txt";
    in_file.open(emission_name.c_str());
    phtn_momw_0_2 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "w/phtn_w_nrg_2-20.txt";
    in_file.open(emission_name.c_str());
    phtn_momw_2_20 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "w/phtn_w_nrg_20-80.txt";
    in_file.open(emission_name.c_str());
    phtn_momw_20_80 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "w/phtn_w_nrg_80-600.txt";
    in_file.open(emission_name.c_str());
    phtn_momw_80_600 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "w/phtn_w_nrg_600-2000.txt";
    in_file.open(emission_name.c_str());
    phtn_momw_600_2000 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    //prepare photon emission nrg lookup table
    Chebyshev_User_2D *phtn_momw_rates[4] = {&phtn_momw_80_600,
                                            &phtn_momw_20_80,
                                            &phtn_momw_2_20,
                                            &phtn_momw_0_2};

    look_up_table_phtn_momw[0] = &phtn_momw_600_2000;                                            

    for(int i = 0; i < 4; i++){
        tmp_i = 1 << i;
        for(int j = 0; j < tmp_i; j++){
            look_up_table_phtn_momw[tmp_i + j] = phtn_momw_rates[i];
        }
    }

    // look_up_table_phtn_momw[0] = &phtn_momw_600_2000;
    // look_up_table_phtn_momw[1] = &phtn_momw_80_600,
    // look_up_table_phtn_momw[2] = &phtn_momw_20_80,
    // look_up_table_phtn_momw[3] = &phtn_momw_2_20,
    // look_up_table_phtn_momw[4] = &phtn_momw_0_2;


    //photon energy in 1/w
    emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_0-2.txt";
    in_file.open(emission_name.c_str());
    phtn_1_over_w_0_2 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_2-20.txt";
    in_file.open(emission_name.c_str());
    phtn_1_over_w_2_20 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_20-80.txt";
    in_file.open(emission_name.c_str());
    phtn_1_over_w_20_80 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_80-600.txt";
    in_file.open(emission_name.c_str());
    phtn_1_over_w_80_600 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_600-2000.txt";
    in_file.open(emission_name.c_str());
    phtn_1_over_w_600_2000 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    //prepare photon emission nrg lookup table
    Chebyshev_User_2D *phtn_1_over_w[4] = {&phtn_1_over_w_80_600,
                                            &phtn_1_over_w_20_80,
                                            &phtn_1_over_w_2_20,
                                            &phtn_1_over_w_0_2};

    look_up_table_1_over_w[0] = &phtn_1_over_w_600_2000;                                            

    for(int i = 0; i < 4; i++){
        tmp_i = 1 << i;
        for(int j = 0; j < tmp_i; j++){
            look_up_table_1_over_w[tmp_i + j] = phtn_1_over_w[i];
        }
    }

    // look_up_table_1_over_w[0] = &phtn_1_over_w_600_2000;
    // look_up_table_1_over_w[1] = &phtn_1_over_w_80_600;
    // look_up_table_1_over_w[2] = &phtn_1_over_w_20_80;
    // look_up_table_1_over_w[3] = &phtn_1_over_w_2_20;
    // look_up_table_1_over_w[4] = &phtn_1_over_w_0_2;



    phtn_1_over_w_proj_0_2.a = phtn_1_over_w_0_2.a;
    phtn_1_over_w_proj_0_2.b = phtn_1_over_w_0_2.b;
    phtn_1_over_w_proj_0_2.domain_half_range = phtn_1_over_w_0_2.domain_half_range_N;
    phtn_1_over_w_proj_0_2.domain_middle_point = phtn_1_over_w_0_2.domain_middle_point_N;
    phtn_1_over_w_proj_0_2.evaluation_order = phtn_1_over_w_0_2.evaluation_order_N;
    phtn_1_over_w_proj_0_2.last_coeffs = phtn_1_over_w_0_2.evaluate_y(heval_r_0_2);

    phtn_1_over_w_proj_2_20.a = phtn_1_over_w_2_20.a;
    phtn_1_over_w_proj_2_20.b = phtn_1_over_w_2_20.b;
    phtn_1_over_w_proj_2_20.domain_half_range = phtn_1_over_w_2_20.domain_half_range_N;
    phtn_1_over_w_proj_2_20.domain_middle_point = phtn_1_over_w_2_20.domain_middle_point_N;
    phtn_1_over_w_proj_2_20.evaluation_order = phtn_1_over_w_2_20.evaluation_order_N;
    phtn_1_over_w_proj_2_20.last_coeffs = phtn_1_over_w_2_20.evaluate_y(heval_r_2_20);

    phtn_1_over_w_proj_20_80.a = phtn_1_over_w_20_80.a;
    phtn_1_over_w_proj_20_80.b = phtn_1_over_w_20_80.b;
    phtn_1_over_w_proj_20_80.domain_half_range = phtn_1_over_w_20_80.domain_half_range_N;
    phtn_1_over_w_proj_20_80.domain_middle_point = phtn_1_over_w_20_80.domain_middle_point_N;
    phtn_1_over_w_proj_20_80.evaluation_order = phtn_1_over_w_20_80.evaluation_order_N;
    phtn_1_over_w_proj_20_80.last_coeffs = phtn_1_over_w_20_80.evaluate_y(heval_r_20_80);

    phtn_1_over_w_proj_80_600.a = phtn_1_over_w_80_600.a;
    phtn_1_over_w_proj_80_600.b = phtn_1_over_w_80_600.b;
    phtn_1_over_w_proj_80_600.domain_half_range = phtn_1_over_w_80_600.domain_half_range_N;
    phtn_1_over_w_proj_80_600.domain_middle_point = phtn_1_over_w_80_600.domain_middle_point_N;
    phtn_1_over_w_proj_80_600.evaluation_order = phtn_1_over_w_80_600.evaluation_order_N;
    phtn_1_over_w_proj_80_600.last_coeffs = phtn_1_over_w_80_600.evaluate_y(heval_r_80_600);

    phtn_1_over_w_proj_600_2000.a = phtn_1_over_w_600_2000.a;
    phtn_1_over_w_proj_600_2000.b = phtn_1_over_w_600_2000.b;
    phtn_1_over_w_proj_600_2000.domain_half_range = phtn_1_over_w_600_2000.domain_half_range_N;
    phtn_1_over_w_proj_600_2000.domain_middle_point = phtn_1_over_w_600_2000.domain_middle_point_N;
    phtn_1_over_w_proj_600_2000.evaluation_order = phtn_1_over_w_600_2000.evaluation_order_N;
    phtn_1_over_w_proj_600_2000.last_coeffs = phtn_1_over_w_600_2000.evaluate_y(heval_r_600_2000);

    Chebyshev_User_1D *phtn_1_o_w_projs[4] = {&phtn_1_over_w_proj_80_600,
                                            &phtn_1_over_w_proj_20_80,
                                            &phtn_1_over_w_proj_2_20,
                                            &phtn_1_over_w_proj_0_2};

    look_up_table_phtn_1_o_w_proj[0] = &phtn_1_over_w_proj_600_2000;                                            

    for(int i = 0; i < 4; i++){
        tmp_i = 1 << i;
        for(int j = 0; j < tmp_i; j++){
            look_up_table_phtn_1_o_w_proj[tmp_i + j] = phtn_1_o_w_projs[i];
        }
    }

    look_up_table_phtn_1_o_w_proj[0] = &phtn_1_over_w_proj_600_2000;  
    look_up_table_phtn_1_o_w_proj[1] = &phtn_1_over_w_proj_80_600;
    look_up_table_phtn_1_o_w_proj[2] = &phtn_1_over_w_proj_20_80;
    look_up_table_phtn_1_o_w_proj[3] = &phtn_1_over_w_proj_2_20;
    look_up_table_phtn_1_o_w_proj[4] = &phtn_1_over_w_proj_0_2;

    //differential cross section
    emission_name = path_to_coeffs + "dP_over_dt_dw/phtn_diff_crss_sctn_0-2.txt";
    // emission_name = path_to_coeffs + "prtl_rate_complete/phtn_partial_rate_0-2.txt";
    in_file.open(emission_name.c_str());
    phtn_diff_crss_sctn_0_2 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "dP_over_dt_dw/phtn_diff_crss_sctn_2-20.txt";
    // emission_name = path_to_coeffs + "prtl_rate_complete/phtn_partial_rate_2-20.txt";
    in_file.open(emission_name.c_str());
    phtn_diff_crss_sctn_2_20 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "dP_over_dt_dw/phtn_diff_crss_sctn_20-80.txt";
    // emission_name = path_to_coeffs + "prtl_rate_complete/phtn_partial_rate_20-80.txt";
    in_file.open(emission_name.c_str());
    phtn_diff_crss_sctn_20_80 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "dP_over_dt_dw/phtn_diff_crss_sctn_80-600.txt";
    // emission_name = path_to_coeffs + "prtl_rate_complete/phtn_partial_rate_80-600.txt";
    in_file.open(emission_name.c_str());
    phtn_diff_crss_sctn_80_600 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "dP_over_dt_dw/phtn_diff_crss_sctn_600-2000.txt";
    // emission_name = path_to_coeffs + "prtl_rate_complete/phtn_partial_rate_600-2000.txt";
    in_file.open(emission_name.c_str());
    phtn_diff_crss_sctn_600_2000 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    //prepare photon emission nrg lookup table
    Chebyshev_User_2D *phtn_diff_crss_sctn[4] = {&phtn_diff_crss_sctn_80_600,
                                            &phtn_diff_crss_sctn_20_80,
                                            &phtn_diff_crss_sctn_2_20,
                                            &phtn_diff_crss_sctn_0_2};

    look_up_table_phtn_diff_crss_sctn[0] = &phtn_diff_crss_sctn_600_2000;                                            

    for(int i = 0; i < 4; i++){
        tmp_i = 1 << i;
        for(int j = 0; j < tmp_i; j++){
            look_up_table_phtn_diff_crss_sctn[tmp_i + j] = phtn_diff_crss_sctn[i];
        }
    }

    // look_up_table_phtn_diff_crss_sctn[0] = &phtn_diff_crss_sctn_600_2000;   
    // look_up_table_phtn_diff_crss_sctn[1] = &phtn_diff_crss_sctn_80_600;
    // look_up_table_phtn_diff_crss_sctn[2] = &phtn_diff_crss_sctn_20_80;
    // look_up_table_phtn_diff_crss_sctn[3] = &phtn_diff_crss_sctn_2_20;
    // look_up_table_phtn_diff_crss_sctn[4] = &phtn_diff_crss_sctn_0_2;
}

//the function below returns the rate of photon emission for an electron
//having  a certain chi and gamma (gamma is the relativistic factor of 
//the electron, but it can also be interpreted as its energy normalized in 
//units of m_e c^2)
double SFQED_Processes::SFQED_PHTN_emission_rate(const double &gamma, const double &chi) const{
    
    //coefficient for the rate of emission
    //multiply by \tilde{W_rad} to get the rate of emission
    double coefW_rad = coef_rate_Wrad * chi / gamma;
    
    //low impact Branch to calculate \tilde{W_rad} and then the rate
    //---------------------------------------------------
    /**/
    //lookup tables boolean selectors
    bool chi_0_2 = chi <= bound_chi_1st,
        chi_2_20 = chi <= bound_chi_2nd,
        chi_20_80 = chi <= bound_chi_3rd,
        chi_80_600 = chi <= bound_chi_4th;
        //chi_600_2000 = chi <= 2000.;

    int lookup_index = chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600;
    //lookup_index = (1 - std::signbit(lookup_index-1))*(32 - __builtin_clz(lookup_index));

    //IMPORTANT: through this constructs we means
    //to default every chi > 600 case into the
    //600 < chi <= 2000 one
    const Chebyshev_User_1D * const chebyshev_phtn_rate = (look_up_table_phtn_mx_rt[lookup_index]);
    

    /*********************/
    /* NO LOOK UP TABLES */
    /*********************/
    //if you decide to use this version please comment every
    //statement relying on lookup_index

    // static Chebyshev_User_1D chebyshev_phtn_rate;

    // if(chi <= bound_chi_1st){
    //     chebyshev_phtn_rate = phtn_mx_rt_0_2;
    // }
    // else if(chi <= bound_chi_2nd){
    //     chebyshev_phtn_rate = phtn_mx_rt_2_20;
    // }
    // else if(chi <= bound_chi_3rd){
    //     chebyshev_phtn_rate = phtn_mx_rt_20_80;
    // }
    // else if(chi <= bound_chi_4th){
    //     chebyshev_phtn_rate = phtn_mx_rt_80_600;
    // }
    // else{
    //     chebyshev_phtn_rate = phtn_mx_rt_600_2000;
    // }
    //////////////////////////////////////

    return coefW_rad * chebyshev_phtn_rate->evaluate(chi);//

    //standard if branching
    // if (chi_0_2){
    //     return coefW_rad*phtn_mx_rt_0_2.evaluate(chi);
    // }
    
    // if (chi_2_20){
    //     return coefW_rad*phtn_mx_rt_2_20.evaluate(chi);
    // }
    
    // if (chi_20_80){
    //     return coefW_rad*phtn_mx_rt_20_80.evaluate(chi);
    // }
       
    // if (chi_80_600) {
    //     return coefW_rad*phtn_mx_rt_80_600.evaluate(chi);
    // }
    
    // //default case
    // return coefW_rad*phtn_mx_rt_600_2000.evaluate(chi);
}

double SFQED_Processes::SFQED_LCFA_phtn_nrg_aux(const double &chi, const double &rnd, const int& lookup_index) const{

    double v;

    //std::cout << rnd << " " << lookup_index << " 0\n" << std::flush;

    /*********************/
    /* NO LOOK UP TABLES */
    /*********************/
    /*
    //if you decide to use this version please comment every
    //statement relying on lookup_index

    static double low_r_limit,
            inverse_r_limit,
            high_r_limit;
    static Chebyshev_User_1D chebyshev_phtn_rate, chebyshev_phtn_1_o_w_proj;
    static Chebyshev_User_2D chebyshev_phtn_momw, chebyshev_phtn_1_o_w, chebyshev_phtn_prtl_rate;

    if(chi <= bound_chi_1st){
        low_r_limit = lower_r_0_2;
        inverse_r_limit = inverse_zone_r_0_2;
        high_r_limit = upper_r_0_2;
        chebyshev_phtn_rate = phtn_mx_rt_0_2;
        chebyshev_phtn_1_o_w_proj = phtn_1_over_w_proj_0_2;
        chebyshev_phtn_momw = phtn_momw_0_2;
        chebyshev_phtn_1_o_w = phtn_1_over_w_0_2;
        chebyshev_phtn_prtl_rate = phtn_prtl_rt_0_2;
    }
    else if(chi <= bound_chi_2nd){
        low_r_limit = lower_r_2_20;
        inverse_r_limit = inverse_zone_r_2_20;
        high_r_limit = upper_r_2_20;
        chebyshev_phtn_rate = phtn_mx_rt_2_20;
        chebyshev_phtn_1_o_w_proj = phtn_1_over_w_proj_2_20;
        chebyshev_phtn_momw = phtn_momw_2_20;
        chebyshev_phtn_1_o_w = phtn_1_over_w_2_20;
        chebyshev_phtn_prtl_rate = phtn_prtl_rt_2_20;
    }
    else if(chi <= bound_chi_3rd){
        low_r_limit = lower_r_20_80;
        inverse_r_limit = inverse_zone_r_20_80;
        high_r_limit = upper_r_20_80;
        chebyshev_phtn_rate = phtn_mx_rt_20_80;
        chebyshev_phtn_1_o_w_proj = phtn_1_over_w_proj_20_80;
        chebyshev_phtn_momw = phtn_momw_20_80;
        chebyshev_phtn_1_o_w = phtn_1_over_w_20_80;
        chebyshev_phtn_prtl_rate = phtn_prtl_rt_20_80;
    }
    else if(chi <= bound_chi_4th){
        low_r_limit = lower_r_80_600;
        inverse_r_limit = inverse_zone_r_80_600;
        high_r_limit = upper_r_80_600;
        chebyshev_phtn_rate = phtn_mx_rt_80_600;
        chebyshev_phtn_1_o_w_proj = phtn_1_over_w_proj_80_600;
        chebyshev_phtn_momw = phtn_momw_80_600;
        chebyshev_phtn_1_o_w = phtn_1_over_w_80_600;
        chebyshev_phtn_prtl_rate = phtn_prtl_rt_80_600;
    }
    else{
        low_r_limit = lower_r_600_2000;
        inverse_r_limit = inverse_zone_r_600_2000;
        high_r_limit = upper_r_600_2000;
        chebyshev_phtn_rate = phtn_mx_rt_600_2000;
        chebyshev_phtn_1_o_w_proj = phtn_1_over_w_proj_600_2000;
        chebyshev_phtn_momw = phtn_momw_600_2000;
        chebyshev_phtn_1_o_w = phtn_1_over_w_600_2000;
        chebyshev_phtn_prtl_rate = phtn_prtl_rt_600_2000;
    }
    */
    //////////////////////////////////////

    const double low_r_limit = *(lu_table_lower_r_bounds[lookup_index]),
            inverse_r_limit = *(lu_table_inverse_r_bounds[lookup_index]),
            high_r_limit = *(lu_table_upper_r_bounds[lookup_index]);

    
    //std::cout << " 1\n" << std::flush;

    //low energy tail: if the random number is smaller than 
    //the corresponding low energy limit for r, then we can employ the
    //soft photon approximation
    if(rnd < low_r_limit){

        //std::cout << "2\n" << std::flush;

        //////////////////////////////////////
        //this part is savagely repeated in the first
        //two 'if' scenarios. Even though this appears
        //quite unelegant, we expect the branch
        //prediction to always end up in the last case,
        //where this section is not needed
        const Chebyshev_User_1D * const chebyshev_phtn_rate =
                (look_up_table_phtn_mx_rt[lookup_index]);

        //std::cout << "3\n" << std::flush;

        //photon emission rate (just the integral without the coefficient \chi / gamma)
        //double tildeWrad = chebyshev_phtn_rate.evaluate(chi);
        double rnd_tildeWrad = rnd * chebyshev_phtn_rate->evaluate(chi);
        //////////////////////////////////////

        //std::cout << "4\n" << std::flush;

        //this is w
        v = soft_ph_coef * rnd_tildeWrad;
        //this is v = w^3
        // v *= v*v;
    }
    //otherwise we use the inverse equation Chebyshev approximation
    else if (rnd < inverse_r_limit){

        //std::cout << "11\n" << std::flush;

        //this evaluation, inverse of g(\chi, w) - r*W_{rad}(\chi) = 0,
        //where \chi and r are fixed at this stage, will give a w value
        const Chebyshev_User_2D * const chebyshev_phtn_momw =
                (look_up_table_phtn_momw[lookup_index]);

        //std::cout << "12\n" << std::flush;

        //this is w
        v = chebyshev_phtn_momw->evaluate(chi, rnd);

        //std::cout << "13\n" << std::flush;

        //this is v
        // v *= v*v;

        //debug
        // std::cout << "[MEDIUM] for chi = " << chi << ", r = " << rnd 
        //         << ", we evaluated a v = " << v << "\n" << std::flush;
    }
    //in case the random number corresponds to some high energy value we resort to the inverse Chebyshev approx
    else if(rnd < high_r_limit){
        const Chebyshev_User_2D * const chebyshev_phtn_1_o_w =
                (look_up_table_1_over_w[lookup_index]);

        //this is w
        v = 1. / chebyshev_phtn_1_o_w->evaluate(chi, rnd);

        //this is v
        // v *= v*v;
    }
    //in the off chance that we happen to be beyond the
    //high energy tail we will employ the exponential 
    //approximation (see notes)
    else{

        //std::cout << "5\n" << std::flush;

        //////////////////////////////////////
        //this part is savagely repeated in the first
        //two 'if' scenarios. Even though this appears
        //quite unelegant, we expect the branch
        //prediction to always end up in the last case,
        //where this section is not needed
        const Chebyshev_User_1D * const chebyshev_phtn_rate_high =
                (look_up_table_phtn_mx_rt[lookup_index]);

        //std::cout << "6\n" << std::flush;

        //photon emission rate (just the integral without the coefficient \chi / gamma)
        double tildeWrad = chebyshev_phtn_rate_high->evaluate(chi);
        double rnd_tildeWrad = rnd * tildeWrad;
        //////////////////////////////////////

        //std::cout << "7\n" << std::flush;

        //get the partial photon emission rate, by integrating
        //the differential probability from 0 to a certain w
        //(keep chi fixed). Please, notice we are evaluating
        //this on the upper bound of validity of the approximation
        const Chebyshev_User_2D * const chebyshev_phtn_prtl_rate =
                (look_up_table_phtn_prtl_rt[lookup_index]);

        //std::cout << "8\n" << std::flush;

        //fixed pivots******************************
        // double w0 = *(lu_high_tail_w_bounds[lookup_index]),
        //         v0 = *(lu_high_tail_v_bounds[lookup_index]);
        //t*****************************************
        
        //newly computed pivots/////////////////////
        //this new attempt requires us to compute w through the inverse (maybe it is too much time consuming?)
        const Chebyshev_User_1D * const chebyshev_phtn_1_o_w_proj =
                (look_up_table_phtn_1_o_w_proj[lookup_index]);
        double w0 = 1. / chebyshev_phtn_1_o_w_proj->evaluate(chi);
        double v0 = w0 * w0 * w0;

        // Chebyshev_User_2D chebyshev_phtn_momw =
        //         *(look_up_table_phtn_momw[lookup_index]);
        // double w0 =  chebyshev_phtn_momw.evaluate(chi, high_r_limit);
        // double v0 = w0 * w0 * w0;
        ////////////////////////////////////////////

        double g_of_w_chi = chebyshev_phtn_prtl_rate->evaluate(chi, w0);

        //std::cout << "9\n" << std::flush;

        //be careful to when v is null
        v = std::abs((tildeWrad - rnd_tildeWrad) / (tildeWrad - g_of_w_chi));
                
        //we use the cubic root to retrieve the w value
        v = cbrt(v0 - log(v));

        //debug
        // std::cout << "[HIGH] for chi = " << chi << ", r = " << rnd << " and w0 = " << w0
        //         << ", we had a tildeWrad = " << tildeWrad << ", a g_of_w_chi = " << g_of_w_chi <<
        //         ", and thus a v = " << v << "\n" << std::flush;
    }

    //return the photon energy disguised as w value
    return v;
}  

void SFQED_Processes::SFQED_collinear_momentum(const double &gamma_out, const double p_in[3], double p_out[3]){
    
    double v = sqrt( (gamma_out*gamma_out - 1.) / (p_out[0]*p_out[0] + p_out[1]*p_out[1] + p_out[2]*p_out[2]) );

    //set the emitted photon momentum by considering that
    //p_\gamma \approx \varepsilon_\gamma (p / |p|) \approx
    //(\varepsilon_\gamma / \varepsilon) p = [ 3 \chi v / (2 + 3 \chi v) ] * p
    //co-linear emission
    p_out[0] = v*p_in[0];
    p_out[1] = v*p_in[1];
    p_out[2] = v*p_in[2];

    return;
}

//return the emitted photon energy (normalized in units of m_e c^2)
double SFQED_Processes::SFQED_LCFA_emitted_photon_energy(const double &gamma, const double &chi, const double &rnd) const{

    //BRANCH ACCORDING TO THE VALUE OF \chi
    //-------------------------------------
    //lookup tables boolean selectors
    bool chi_0_2 = chi <= bound_chi_1st,
        chi_2_20 = chi <= bound_chi_2nd,
        chi_20_80 = chi <= bound_chi_3rd,
        chi_80_600 = chi <= bound_chi_4th;
        //chi_600_2000 = chi <= 2000.;

    int lookup_index = chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600;
    //lookup_index = (1 - std::signbit(lookup_index-1))*(32 - __builtin_clz(lookup_index));

    //alternative way to compute the index and use the lutables
    //int lookup_index = static_cast<int>(log2(chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600 + 1));

    //IMPORTANT: through the above construct we mean
    //to default every chi > 600 case into the
    //600 < chi <= 2000 one
    
    double v = SFQED_LCFA_phtn_nrg_aux(chi, rnd, lookup_index);

    //CALCULATE THE PHOTON MOMENTUM AND ENERGY
    //in normalized units (i.e., in units of m_e c^2) by using:
    //( \varepsilon_\gamma / \varepsilon ) = 3 \chi v / (2 + 3 \chi v)
    v = 3.0*chi*v*v*v;
    v = gamma * v / (2.0 + v);

    return v;
}

double determinant(const double& w, const double& chi){
    double aux = w / (1. + 1.5 * chi * w * w * w);

    return aux * aux;
}

double determinant1(const double& phtn_nrg, const double& chi, const double& gamma){
    return gamma * cbrt((2.)/(81. * phtn_nrg * phtn_nrg * chi * (gamma - phtn_nrg) * (gamma - phtn_nrg) * (gamma - phtn_nrg) * (gamma - phtn_nrg)));
}

//\gamma here represents the normalized energy of the particle \epsilon,
//that in case of an electron \epsilon = \gamma. In general it is not
//the sheer gamma factor
//if an increasing tail appears at low energies (log scale) please check the conditions
double SFQED_Processes::SFQED_BLCFA_emitted_photon_energy(const double& LCFA_limit,
                                                            const double& gamma, 
                                                            const double& chi, 
                                                            const double& rnd, 
                                                            const double& rnd2) const{

    // static double increment = 0.01;

    //key generating section
    bool chi_0_2 = chi <= bound_chi_1st,
        chi_2_20 = chi <= bound_chi_2nd,
        chi_20_80 = chi <= bound_chi_3rd,
        chi_80_600 = chi <= bound_chi_4th;
        //chi_600_2000 = chi <= 2000.;

    int lookup_index = chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600;
    //lookup_index = (1 - std::signbit(lookup_index-1))*(32 - __builtin_clz(lookup_index));
    
    //transform the threshold energy into a w value
    double LCFA_limit_w  = cbrt(LCFA_limit / (1.5 * chi * (gamma - LCFA_limit)));
    

    //in case the phtn nrg threshold is very close to the particle nrg
    //or it is above the upper bound of the approximated domain
    //we don't want the emission to occur 
    if(LCFA_limit > 0.75 * gamma || LCFA_limit_w > *(lu_table_upper_w_bounds[lookup_index])){
        return 0.;
    }

    //compute the w-value energy of the LCFA emitted photon
    double emitted_phtn_LCFA_nrg_w = SFQED_LCFA_phtn_nrg_aux(chi, rnd, lookup_index),
            diff_probability,
            diff_prob_LCFA,
            rs_dw_over_depsgamma;

    /////////////////////////////////////////
    //  USUAL METHOD

    // double w_validity = *(lu_table_upper_w_bounds[lookup_index]);
    // //if the w value of the energy threshold is above the limit
    // //of validity the emission is suppressed by default
    // if(LCFA_limit_w > w_validity){
    //     return 0.;
    // }

    //BLCFA section (it starts with a comparison)
    if(emitted_phtn_LCFA_nrg_w < LCFA_limit_w){

        const Chebyshev_User_2D * const phtn_diff_crss_sctn =
                (look_up_table_phtn_diff_crss_sctn[lookup_index]);

        //compute the threshold value
        rs_dw_over_depsgamma = (1. + 1.5 * chi * LCFA_limit_w * LCFA_limit_w * LCFA_limit_w) / LCFA_limit_w;
        diff_prob_LCFA = phtn_diff_crss_sctn->evaluate(chi, LCFA_limit_w) * rs_dw_over_depsgamma * rs_dw_over_depsgamma;// / determinant(LCFA_limit_w, chi);

        //compute the differential probability associated to the LCFA supposed phtn nrg
        //(pay attention not to overwrite the w value energy)
        rs_dw_over_depsgamma = (1. + 1.5 * chi * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w) / emitted_phtn_LCFA_nrg_w;
        diff_probability = phtn_diff_crss_sctn->evaluate(chi, emitted_phtn_LCFA_nrg_w) * rs_dw_over_depsgamma * rs_dw_over_depsgamma;// / determinant(emitted_phtn_LCFA_nrg_w, chi);

        //now I use the rejection sampling method (employing rnd2) to determine if I have
        //to keep or reject the just occurred emission
        // if(diff_probability * rnd2 > diff_prob_LCFA ){
        if(diff_probability * rnd2 > diff_prob_LCFA){
            return 0.;
        }

    }

    //convert to a proper photon energy
    emitted_phtn_LCFA_nrg_w = 3.0 * chi * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w;
    emitted_phtn_LCFA_nrg_w = gamma * emitted_phtn_LCFA_nrg_w / (2.0 + emitted_phtn_LCFA_nrg_w);

    return emitted_phtn_LCFA_nrg_w;
    ///////////////////////////////////////////////

    /////////////////////////////////////////
    //  2nd METHOD

    // double emitted_phtn_LCFA_nrg = 3.0 * chi * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w;
    // emitted_phtn_LCFA_nrg = gamma * emitted_phtn_LCFA_nrg / (2.0 + emitted_phtn_LCFA_nrg);

    // //BLCFA section (it starts with a comparison)
    // if(emitted_phtn_LCFA_nrg < LCFA_limit){

    //     Chebyshev_User_2D phtn_diff_crss_sctn =
    //             *(look_up_table_phtn_diff_crss_sctn[lookup_index]);

    //     //copmute the threshold value
    //     diff_prob_LCFA = phtn_diff_crss_sctn.evaluate(chi, LCFA_limit_w) * determinant1(LCFA_limit, chi, gamma);

    //     //compute the differential probability associated to the LCFA supposed phtn nrg
    //     //(pay attention not to overwrite the w value energy)
    //     diff_probability = phtn_diff_crss_sctn.evaluate(chi, emitted_phtn_LCFA_nrg_w) * determinant1(emitted_phtn_LCFA_nrg, chi, gamma);

    //     //now I use the rejection sampling method (employing rnd2) to determine if I have
    //     //to keep or reject the just occurred emission
    //     if(diff_probability * rnd2 > diff_prob_LCFA ){
    //         return 0.;
    //     }

    // }

    // return emitted_phtn_LCFA_nrg;
    ///////////////////////////////////////////////

    /////////////////////////////////////////
    //  3rd METHOD (this worked)

    // double emitted_phtn_LCFA_nrg = 3.0 * chi * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w;
    // emitted_phtn_LCFA_nrg = gamma * emitted_phtn_LCFA_nrg / (2.0 + emitted_phtn_LCFA_nrg);

    // //BLCFA section (it starts with a comparison)
    // if(emitted_phtn_LCFA_nrg < LCFA_limit){

    //     Chebyshev_User_2D phtn_diff_crss_sctn =
    //             *(look_up_table_phtn_diff_crss_sctn[lookup_index]);

    //     double en_up, v_up, w_up, en, v, w;
    //     double f_w_up, f_w;
    //     double dW_dw, dw_deg, dWrad_en_cut, dWrad_ph_en;

    //     //**********work with the threshold*****************
    //     en_up = LCFA_limit + increment;
    //     v_up = (2.0 * en_up) / (3.0 * chi * (gamma - en_up));
    //     w_up = cbrt(v_up);
        
    //     en = LCFA_limit;
    //     v = (2.0 * en) / (3.0 * chi * (gamma - en));
    //     w = cbrt(v);

    //     f_w_up = phtn_diff_crss_sctn.evaluate(chi, w_up);
    //     f_w = phtn_diff_crss_sctn.evaluate(chi, w);
  
    //     dW_dw = (f_w_up - f_w) / (w_up - w);
    //     dw_deg = 2.0 / (81.0 * chi * en * en * (gamma - en) * (gamma - en) * (gamma - en) * (gamma - en));
    //     dw_deg = cbrt(dw_deg) * gamma;
    //     dWrad_en_cut = dW_dw * dw_deg;

    //     //***********work with the actual LCFA nrg******************
    //     en_up = emitted_phtn_LCFA_nrg + increment;
    //     v_up = (2.0 * en_up) / (3.0 * chi * (gamma - en_up));
    //     w_up = cbrt(v_up);
          
    //     en = emitted_phtn_LCFA_nrg;
    //     v = (2.0 * en) / (3.0 * chi * (gamma - en));
    //     w = cbrt(v);
          
    //     f_w_up = phtn_diff_crss_sctn.evaluate(chi, w_up);;
    //     f_w = phtn_diff_crss_sctn.evaluate(chi, w);
          
    //     dW_dw = (f_w_up - f_w) / (w_up - w);
    //     dw_deg = 2.0 / (81.0 * chi * en * en * (gamma - en) * (gamma - en) * (gamma - en) * (gamma - en));
    //     dw_deg = cbrt(dw_deg) * gamma;
    //     dWrad_ph_en = dW_dw * dw_deg;

    //     //now I use the rejection sampling method (employing rnd2) to determine if I have
    //     //to keep or reject the just occurred emission
    //     if((rnd2 * dWrad_ph_en) > dWrad_en_cut){
    //         return 0.;
    //     }

    // }

    // return emitted_phtn_LCFA_nrg;
    ///////////////////////////////////////////////

    /////////////////////////////////////////
    //  4th METHOD

    // double emitted_phtn_LCFA_nrg = 3.0 * chi * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w;
    // emitted_phtn_LCFA_nrg = gamma * emitted_phtn_LCFA_nrg / (2.0 + emitted_phtn_LCFA_nrg);

    // //BLCFA section (it starts with a comparison)
    // if(emitted_phtn_LCFA_nrg < LCFA_limit){


    //     double en, v, w;
    //     double dW_dw, dw_deg, dWrad_en_cut, dWrad_ph_en;

    //     //**********work with the threshold*****************

    //     en = LCFA_limit;
    //     v = (2.0 * en) / (3.0 * chi * (gamma - en));
    //     //change w (use the one above directly)
    //     w = cbrt(v);

    //     //change function (use your approximation)
    //     dW_dw = differential_cross_section_phtn_emission(chi, w, NULL);
    //     //change determinant (use the function)
    //     dw_deg = 2.0 / (81.0 * chi * en * en * (gamma - en) * (gamma - en) * (gamma - en) * (gamma - en));
    //     dw_deg = cbrt(dw_deg) * gamma;
    //     dWrad_en_cut = dW_dw * dw_deg;

    //     //***********work with the actual LCFA nrg******************
          
    //     en = emitted_phtn_LCFA_nrg;
    //     v = (2.0 * en) / (3.0 * chi * (gamma - en));
    //     w = cbrt(v);
          
    //     dW_dw = differential_cross_section_phtn_emission(chi, w, NULL);
    //     dw_deg = 2.0 / (81.0 * chi * en * en * (gamma - en) * (gamma - en) * (gamma - en) * (gamma - en));
    //     dw_deg = cbrt(dw_deg) * gamma;
    //     dWrad_ph_en = dW_dw * dw_deg;

    //     //now I use the rejection sampling method (employing rnd2) to determine if I have
    //     //to keep or reject the just occurred emission
    //     if((rnd2 * dWrad_ph_en) > dWrad_en_cut){
    //         return 0.;
    //     }

    // }

    // return emitted_phtn_LCFA_nrg;
    ///////////////////////////////////////////////

    /////////////////////////////////////////
    //  5th METHOD

    // double emitted_phtn_LCFA_nrg = 3.0 * chi * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w;
    // emitted_phtn_LCFA_nrg = gamma * emitted_phtn_LCFA_nrg / (2.0 + emitted_phtn_LCFA_nrg);

    // //BLCFA section (it starts with a comparison)
    // if(emitted_phtn_LCFA_nrg < LCFA_limit){


    //     Chebyshev_User_2D phtn_diff_crss_sctn = phtn_diff_crss_sctn_0_2;
    //             // *(look_up_table_phtn_diff_crss_sctn[lookup_index]);

    //     double en, v, w;
    //     double dW_dw, dw_deg, dWrad_en_cut, dWrad_ph_en;

    //     //**********work with the threshold*****************

    //     en = LCFA_limit;
    //     v = (2.0 * en) / (3.0 * chi * (gamma - en));
    //     //change w (use the one above directly)
    //     w = cbrt(v);

    //     //change function (use your approximation)
    //     dW_dw = phtn_diff_crss_sctn.evaluate(chi, w);
    //     //change determinant (use the function)
    //     dw_deg = 2.0 / (81.0 * chi * en * en * (gamma - en) * (gamma - en) * (gamma - en) * (gamma - en));
    //     dw_deg = cbrt(dw_deg) * gamma;
    //     dWrad_en_cut = dW_dw * dw_deg;

    //     //***********work with the actual LCFA nrg******************
          
    //     en = emitted_phtn_LCFA_nrg;
    //     v = (2.0 * en) / (3.0 * chi * (gamma - en));
    //     w = cbrt(v);
          
    //     dW_dw = phtn_diff_crss_sctn.evaluate(chi, w);
    //     dw_deg = 2.0 / (81.0 * chi * en * en * (gamma - en) * (gamma - en) * (gamma - en) * (gamma - en));
    //     dw_deg = cbrt(dw_deg) * gamma;
    //     dWrad_ph_en = dW_dw * dw_deg;

    //     //now I use the rejection sampling method (employing rnd2) to determine if I have
    //     //to keep or reject the just occurred emission
    //     if((rnd2 * dWrad_ph_en) > dWrad_en_cut){
    //         return 0.;
    //     }

    // }

    // return emitted_phtn_LCFA_nrg;
    ///////////////////////////////////////////////
}


bool BLCFA_Object::SFQED_BLCFA_update_entities_quantities(const SFQED_Processes& proc,
                                                const double* const pushed_momentum,
                                                const double* const momentum,
                                                double* Lorentz_F_t_der,
                                                double* Lorentz_F_tt_der,
                                                double& part_gamma,
                                                double& part_chi){

    // ///////////////////////////////////////////
    // // FIRST METHOD ///////////////////////////
    // ///////////////////////////////////////////

    //even though this external pointer variable is meant to store the
    //force derivative value, we will use it to temporarily deal with
    //the actual Lorentz force vector
    Lorentz_F_t_der[0] = proc.one_over_dt * (pushed_momentum[0] - momentum[0]),
    Lorentz_F_t_der[1] = proc.one_over_dt * (pushed_momentum[1] - momentum[1]),
    Lorentz_F_t_der[2] = proc.one_over_dt * (pushed_momentum[2] - momentum[2]);

    double momentum_sum_x = (pushed_momentum[0] + momentum[0]),
            momentum_sum_y = (pushed_momentum[1] + momentum[1]),
            momentum_sum_z = (pushed_momentum[2] + momentum[2]);

    //compute the proper gamma factor. we use this to approximate the momentum modulus
    //and avoid possible numerical issue
    part_gamma = 4. + (momentum_sum_x*momentum_sum_x + momentum_sum_y*momentum_sum_y + momentum_sum_z*momentum_sum_z);

    //this variable is called chi, but right now it stores the
    //projection factor of the Lorentz force along the 
    //momentum_sum vector direction
    part_chi = (Lorentz_F_t_der[0]*momentum_sum_x + Lorentz_F_t_der[1]*momentum_sum_y + Lorentz_F_t_der[2]*momentum_sum_z) / part_gamma;

    //compute the proper gamma [we also update a variable passed from outside]
    part_gamma = std::sqrt(part_gamma*0.25);

    //now we replace the actual Lorentz force comps with those
    //perpendicular to the momentum
    Lorentz_F_t_der[0] = Lorentz_F_t_der[0] - part_chi * momentum_sum_x,
    Lorentz_F_t_der[1] = Lorentz_F_t_der[1] - part_chi * momentum_sum_y,
    Lorentz_F_t_der[2] = Lorentz_F_t_der[2] - part_chi * momentum_sum_z;

    //compute the squared modulus of the orthogonal Lorentz force
    //multiplied by the squared gamma factor
    double L_F_Modulus = Lorentz_F_t_der[0]*Lorentz_F_t_der[0] +
                                Lorentz_F_t_der[1]*Lorentz_F_t_der[1] + 
                                Lorentz_F_t_der[2]*Lorentz_F_t_der[2];

    //we now use the previously defined momentum sum components to store
    //the new orthogonal lorentz force difference
    momentum_sum_x = Lorentz_F_t_der[0] - this->Lorentz_F_Old[0],
    momentum_sum_y = Lorentz_F_t_der[1] - this->Lorentz_F_Old[1],
    momentum_sum_z = Lorentz_F_t_der[2] - this->Lorentz_F_Old[2];

    //update the old lorentz force quantities
    this->Lorentz_F_Old[0] = Lorentz_F_t_der[0];
    this->Lorentz_F_Old[1] = Lorentz_F_t_der[1];
    this->Lorentz_F_Old[2] = Lorentz_F_t_der[2];

    //compute second derivatives
    Lorentz_F_tt_der[0] = proc.one_over_dt2 * (momentum_sum_x - this->Delta_Lorentz_F_Old[0]),
    Lorentz_F_tt_der[1] = proc.one_over_dt2 * (momentum_sum_y - this->Delta_Lorentz_F_Old[1]),
    Lorentz_F_tt_der[2] = proc.one_over_dt2 * (momentum_sum_z - this->Delta_Lorentz_F_Old[2]);

    //update old lorentz force difference
    this->Delta_Lorentz_F_Old[0] = momentum_sum_x;
    this->Delta_Lorentz_F_Old[1] = momentum_sum_y;
    this->Delta_Lorentz_F_Old[2] = momentum_sum_z;

    //compute the actual first derivatives
    Lorentz_F_t_der[0] = proc.one_over_dt * momentum_sum_x;
    Lorentz_F_t_der[1] = proc.one_over_dt * momentum_sum_y;
    Lorentz_F_t_der[2] = proc.one_over_dt * momentum_sum_z;

    //compute the ultimate chi [do we update also a variable passed from outside?]
    part_chi = part_gamma * std::sqrt(proc.norm_Compton_time_2 * L_F_Modulus);//* (mass / charge); we decided to provide these methods for positrons and electrons only

    //if the particle has just been created we exit 
    //immediately after we updated the FORCES, GAMMA and CHI! [IMPORTANT!!!]

    //this return value is used to infer whether a particle has just been created
    if(this->just_created){
       this->just_created = false;
       return false; 
    }
    
    return true;

    ////////////////////////////////////////////
    // SECOND METHOD ///////////////////////////
    ////////////////////////////////////////////
    /*
    //even though this external pointer variable is meant to store the
    //force derivative value, we will use it to temporarily deal with
    //the actual Lorentz force vector
    Lorentz_F_t_der[0] = proc.one_over_dt * (pushed_momentum[0] - momentum[0]),
    Lorentz_F_t_der[1] = proc.one_over_dt * (pushed_momentum[1] - momentum[1]),
    Lorentz_F_t_der[2] = proc.one_over_dt * (pushed_momentum[2] - momentum[2]);

    //MID MOMENTUM
    double momentum_sum_x = (pushed_momentum[0] + momentum[0]) * 0.5,
            momentum_sum_y = (pushed_momentum[1] + momentum[1]) * 0.5,
            momentum_sum_z = (pushed_momentum[2] + momentum[2]) * 0.5;

    //square of the gamma
    part_gamma = 1. + (momentum_sum_x*momentum_sum_x + momentum_sum_y*momentum_sum_y + momentum_sum_z*momentum_sum_z);

    //this variable is called chi, but right now it stores the
    //projection factor of the Lorentz force along the 
    //momentum_sum vector direction
    part_chi = (Lorentz_F_t_der[0]*momentum_sum_x + Lorentz_F_t_der[1]*momentum_sum_y + Lorentz_F_t_der[2]*momentum_sum_z) / part_gamma;

    //compute the proper gamma [we also update a variable passed from outside]
    part_gamma = std::sqrt(part_gamma);

    //now we replace the actual Lorentz force comps with those
    //perpendicular to the momentum
    Lorentz_F_t_der[0] = Lorentz_F_t_der[0] - part_chi * momentum_sum_x,
    Lorentz_F_t_der[1] = Lorentz_F_t_der[1] - part_chi * momentum_sum_y,
    Lorentz_F_t_der[2] = Lorentz_F_t_der[2] - part_chi * momentum_sum_z;

    //compute the squared modulus of the orthogonal Lorentz force
    //multiplied by the squared gamma factor
    double L_F_Modulus = Lorentz_F_t_der[0]*Lorentz_F_t_der[0] +
                                Lorentz_F_t_der[1]*Lorentz_F_t_der[1] + 
                                Lorentz_F_t_der[2]*Lorentz_F_t_der[2];

    //compute second derivatives
    Lorentz_F_tt_der[0] = proc.one_over_dt2 * (Lorentz_F_t_der[0] - 2. * this->Lorentz_F_Old[0] + this->Delta_Lorentz_F_Old[0]);
    Lorentz_F_tt_der[1] = proc.one_over_dt2 * (Lorentz_F_t_der[1] - 2. * this->Lorentz_F_Old[1] + this->Delta_Lorentz_F_Old[1]);
    Lorentz_F_tt_der[2] = proc.one_over_dt2 * (Lorentz_F_t_der[2] - 2. * this->Lorentz_F_Old[2] + this->Delta_Lorentz_F_Old[2]);

    //we now use the previously defined momentum sum components to store
    //the first derivative
    momentum_sum_x = proc.one_over_dt * (Lorentz_F_t_der[0] - this->Lorentz_F_Old[0]);
    momentum_sum_y = proc.one_over_dt * (Lorentz_F_t_der[1] - this->Lorentz_F_Old[1]);
    momentum_sum_z = proc.one_over_dt * (Lorentz_F_t_der[2] - this->Lorentz_F_Old[2]);

    //update old lorentz force difference
    this->Delta_Lorentz_F_Old[0] = this->Lorentz_F_Old[0];
    this->Delta_Lorentz_F_Old[1] = this->Lorentz_F_Old[1];
    this->Delta_Lorentz_F_Old[2] = this->Lorentz_F_Old[2];

    //update the old lorentz force quantities
    this->Lorentz_F_Old[0] = Lorentz_F_t_der[0];
    this->Lorentz_F_Old[1] = Lorentz_F_t_der[1];
    this->Lorentz_F_Old[2] = Lorentz_F_t_der[2];

    //register derivative
    Lorentz_F_t_der[0] = momentum_sum_x;
    Lorentz_F_t_der[1] = momentum_sum_y;
    Lorentz_F_t_der[2] = momentum_sum_z;

    //compute the ultimate chi [do we update also a variable passed from outside?]
    part_chi = part_gamma * std::sqrt(proc.norm_Compton_time_2 * L_F_Modulus);//* (mass / charge); we decided to provide these methods for positrons and electrons only

    //if the particle has just been created we exit 
    //immediately after we updated the FORCES, GAMMA and CHI! [IMPORTANT!!!]

    //this return value is used to infer whether a particle has just been created
    if(this->just_created){
       this->just_created = false;
       return false; 
    }
    
    return true;
    */
}

double BLCFA_Object::SFQED_BLCFA_find_energy_threshold(const SFQED_Processes& proc,
                                                        const double* const Lorentz_F_t_der,
                                                        const double* const Lorentz_F_tt_der,
                                                        const double& part_gamma,
                                                        const double& part_chi) const{
    //first of all: compute the delta (see article)
    double delta = proc.norm_Compton_time_2 * 
                        (Lorentz_F_t_der[0]*Lorentz_F_t_der[0] + Lorentz_F_t_der[1]*Lorentz_F_t_der[1] + Lorentz_F_t_der[2]*Lorentz_F_t_der[2] +
                         std::abs(this->Lorentz_F_Old[0]*Lorentz_F_tt_der[0] + this->Lorentz_F_Old[1]*Lorentz_F_tt_der[1] + this->Lorentz_F_Old[2]*Lorentz_F_tt_der[2]));

    //remember that this represents the squared modulus
    //of the Lorentz force perpendicular to the momentum
    double L_F_Modulus_Sq = this->Lorentz_F_Old[0]*this->Lorentz_F_Old[0]
                            + this->Lorentz_F_Old[1]*this->Lorentz_F_Old[1]
                            + this->Lorentz_F_Old[2]*this->Lorentz_F_Old[2];

    //this variable will be used to store the quantity (\tau)/(\tau_c) * (\chi)/(\gamma)
    //but for the moment it simply stores \tau_C^2 |F_{\perp}|^4 (needed in the control below)
    double formation_time_ratio_mult_chi_over_gamma = proc.norm_Compton_time_2 * L_F_Modulus_Sq * L_F_Modulus_Sq;//multiply this by (m/q)^2 if you want to deal with a general particle

    //the first condition is needed to check if chi is big enough
    //while the second keeps us from encountering numerical issues
    //(like division by 0 or similar) when the fields are constant
    //and the force derivatives are zero.
    if(part_chi <= zeta || delta / zeta_2 <= formation_time_ratio_mult_chi_over_gamma){
        //in this case we return 0 and the LCFA will be used
        return 0.;
    }

    //update the variable formation_time_ratio_mult_chi_over_gamma
    //with that value it is meant to store
    formation_time_ratio_mult_chi_over_gamma = twopi * std::sqrt(formation_time_ratio_mult_chi_over_gamma / delta);

    //debug
    // std::cout << delta << " " << part_gamma << " " << part_chi << " " << formation_time_ratio_mult_chi_over_gamma << " ";

    //this is the normalized LCFA energy threshold
    //[if we want this method to hold for a generic
    // particle of mass m, we should multiply the
    // following expression by m (remember that an electron has m = 1)]
    double LCFA_limit = thresh_factor * part_gamma * part_chi / (part_chi + four_over_threepi * sinh(3. * asinh(formation_time_ratio_mult_chi_over_gamma / 8.)));

    return LCFA_limit;

    // //maybe do this here
    // double LCFA_limit_w  = cbrt(LCFA_limit / (1.5 * part_chi * (part_gamma - LCFA_limit)));
    // //key generating section
    // // bool chi_0_2 = part_chi <= bound_chi_1st,
    // //     chi_2_20 = part_chi <= bound_chi_2nd,
    // //     chi_20_80 = part_chi <= bound_chi_3rd,
    // //     chi_80_600 = part_chi <= bound_chi_4th;
    //     //chi_600_2000 = chi <= 2000.;

    // //in case the phtn nrg threshold is very close to the particle nrg
    // //we don't want the emission to occur 
    // // if(gamma - phtn_nrg < 1e-3 * gamma){
    // if(LCFA_limit > 0.75 * part_gamma){
    // //     LCFA_limit = -1.;
    // // if(LCFA_limit > 0.75 * part_gamma || LCFA_limit_w > *(lu_table_upper_w_bounds[chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600])){
    //     LCFA_limit_w = -1.;
    // }

    // return LCFA_limit_w; //LCFA_limit;
}

void BLCFA_Object::init_optical_depth(double tau_0){
    this->optical_depth_ref = tau_0;
    this->optical_depth = 0.0;
}

void BLCFA_Object::reset_optical_depth(){
    this->optical_depth = 0.0;
}

double BLCFA_Object::increase_optical_depth(double increment){
    return this->optical_depth += increment;
}

bool BLCFA_Object::check_optical_depth_for_emission(){
    return this->optical_depth >= this->optical_depth_ref;
}


void SFQED_Processes::SFQED_finalize_PHTN_emission(){
    //purge emission rates
    delete[] phtn_mx_rt_0_2.last_coeffs;
    delete[] phtn_mx_rt_2_20.last_coeffs;
    delete[] phtn_mx_rt_20_80.last_coeffs;
    delete[] phtn_mx_rt_80_600.last_coeffs;
    delete[] phtn_mx_rt_600_2000.last_coeffs;

    //purge partial (from 0 to w) emission rates
    delete[] phtn_prtl_rt_0_2.last_coeffs;
    delete[] phtn_prtl_rt_2_20.last_coeffs;
    delete[] phtn_prtl_rt_20_80.last_coeffs;
    delete[] phtn_prtl_rt_80_600.last_coeffs;
    delete[] phtn_prtl_rt_600_2000.last_coeffs;

    //purge photon emission energy
    delete[] phtn_momw_0_2.last_coeffs;
    delete[] phtn_momw_2_20.last_coeffs;
    delete[] phtn_momw_20_80.last_coeffs;
    delete[] phtn_momw_80_600.last_coeffs;
    delete[] phtn_momw_600_2000.last_coeffs;

    //purge 1/w
    delete[] phtn_1_over_w_0_2.last_coeffs;
    delete[] phtn_1_over_w_2_20.last_coeffs;
    delete[] phtn_1_over_w_20_80.last_coeffs;
    delete[] phtn_1_over_w_80_600.last_coeffs;
    delete[] phtn_1_over_w_600_2000.last_coeffs;

    //purge 1D projection of 1/w
    delete[] phtn_1_over_w_proj_0_2.last_coeffs;
    delete[] phtn_1_over_w_proj_2_20.last_coeffs;
    delete[] phtn_1_over_w_proj_20_80.last_coeffs;
    delete[] phtn_1_over_w_proj_80_600.last_coeffs;
    delete[] phtn_1_over_w_proj_600_2000.last_coeffs;

    //purge diff crss sctn
    delete[] phtn_diff_crss_sctn_0_2.last_coeffs;
    delete[] phtn_diff_crss_sctn_2_20.last_coeffs;
    delete[] phtn_diff_crss_sctn_20_80.last_coeffs;
    delete[] phtn_diff_crss_sctn_80_600.last_coeffs;
    delete[] phtn_diff_crss_sctn_600_2000.last_coeffs;

    //comment this shit when not needed
    // delete[] intdWrad_0chi2;
}

/*****************/
/*PAIR PRODUCTION*/
/*****************/

double asympt_pair_emission_rate_low_k(double k){
    return (27. * pi * k)/(16. * sqrt(2.)) * exp(-(8/(3 * k))) * (1. - 11./64.*k + 7585./73728.*k*k);
}

double zero_1D(double x){
    return 0.;
}

double zero_2D(double x, double y){
    return 0.;
}

void SFQED_Processes::SFQED_init_PAIR_creation(std::string path_to_coeffs){

    init_pair_crtn_tables();

    int tmp_i;

    //pair creation rate

    pair_prd_rt_001_03.func_to_execute = &asympt_pair_emission_rate_low_k;
    pair_prd_rt_001_03_fast.func_to_execute = &zero_1D;
    
    std::string emission_name = path_to_coeffs + "pair_rate/pair_production_rate_0-2.txt";
    std::ifstream in_file(emission_name.c_str());
    pair_prd_rt_0_2 = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_rate/pair_production_rate_2-20.txt";
    in_file.open(emission_name.c_str());
    pair_prd_rt_2_20 = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_rate/pair_production_rate_20-80.txt";
    in_file.open(emission_name.c_str());
    pair_prd_rt_20_80 = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_rate/pair_production_rate_80-600.txt";
    in_file.open(emission_name.c_str());
    pair_prd_rt_80_600 = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_rate/pair_production_rate_600-2000.txt";
    in_file.open(emission_name.c_str());
    pair_prd_rt_600_2000 = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    //prepare photon emission rate lookup table
    Chebyshev_User_1D *pair_rates[5] = {&pair_prd_rt_80_600,
                                            &pair_prd_rt_20_80,
                                            &pair_prd_rt_2_20,
                                            &pair_prd_rt_0_2,
                                            &pair_prd_rt_001_03};

    look_up_table_pair_prd_rt[0] = &pair_prd_rt_600_2000;                                            

    for(int i = 0; i < 5; i++){
        tmp_i = 1 << i;
        for(int j = 0; j < tmp_i; j++){
            look_up_table_pair_prd_rt[tmp_i + j] = pair_rates[i];
        }
    }

    // look_up_table_pair_prd_rt[0] = &pair_prd_rt_600_2000;
    // look_up_table_pair_prd_rt[1] = &pair_prd_rt_80_600;
    // look_up_table_pair_prd_rt[2] = &pair_prd_rt_20_80;
    // look_up_table_pair_prd_rt[3] = &pair_prd_rt_2_20;
    // look_up_table_pair_prd_rt[4] = &pair_prd_rt_0_2;
    // look_up_table_pair_prd_rt[5] = &pair_prd_rt_001_03;

    //prepare fast lookup table
    pair_rates[4] = &pair_prd_rt_001_03_fast;

    look_up_table_pair_prd_rt_fast[0] = &pair_prd_rt_600_2000;                                            

    for(int i = 0; i < 5; i++){
        tmp_i = 1 << i;
        for(int j = 0; j < tmp_i; j++){
            look_up_table_pair_prd_rt_fast[tmp_i + j] = pair_rates[i];
        }
    }

    // look_up_table_pair_prd_rt_fast[0] = &pair_prd_rt_600_2000;
    // look_up_table_pair_prd_rt_fast[1] = &pair_prd_rt_80_600;
    // look_up_table_pair_prd_rt_fast[2] = &pair_prd_rt_20_80;
    // look_up_table_pair_prd_rt_fast[3] = &pair_prd_rt_2_20;
    // look_up_table_pair_prd_rt_fast[4] = &pair_prd_rt_0_2;
    // look_up_table_pair_prd_rt_fast[5] = &pair_prd_rt_001_03_fast;

//
    pair_prtl_rt_001_03.func_to_execute = &zero_2D;
//

    emission_name = path_to_coeffs + "pair_prtl_rate/pair_partial_rate_0-2.txt";
    in_file.open(emission_name.c_str());
    pair_prtl_rt_0_2 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_prtl_rate/pair_partial_rate_2-20.txt";
    in_file.open(emission_name.c_str());
    pair_prtl_rt_2_20 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_prtl_rate/pair_partial_rate_20-80.txt";
    in_file.open(emission_name.c_str());
    pair_prtl_rt_20_80 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_prtl_rate/pair_partial_rate_80-600.txt";
    in_file.open(emission_name.c_str());
    pair_prtl_rt_80_600 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_prtl_rate/pair_partial_rate_600-2000.txt";
    in_file.open(emission_name.c_str());
    pair_prtl_rt_600_2000 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    //prepare photon emission rate lookup table
    Chebyshev_User_2D *pair_prtl_rates[5] = {&pair_prtl_rt_80_600,
                                            &pair_prtl_rt_20_80,
                                            &pair_prtl_rt_2_20,
                                            &pair_prtl_rt_0_2,
                                            &pair_prtl_rt_001_03};

    look_up_table_pair_prtl_rt[0] = &pair_prtl_rt_600_2000;                                            

    for(int i = 0; i < 5; i++){
        tmp_i = 1 << i;
        for(int j = 0; j < tmp_i; j++){
            look_up_table_pair_prtl_rt[tmp_i + j] = pair_prtl_rates[i];
        }
    }

    // look_up_table_pair_prtl_rt[0] = &pair_prtl_rt_600_2000;  
    // look_up_table_pair_prtl_rt[1] = &pair_prtl_rt_80_600,
    // look_up_table_pair_prtl_rt[2] = &pair_prtl_rt_20_80,
    // look_up_table_pair_prtl_rt[3] = &pair_prtl_rt_2_20,
    // look_up_table_pair_prtl_rt[4] = &pair_prtl_rt_0_2,
    // look_up_table_pair_prtl_rt[5] = &pair_prtl_rt_001_03;


    //photon energy in v

//
    emission_name = path_to_coeffs + "v/pair_nrgs_v_001-03.txt";
    in_file.open(emission_name.c_str());
    pair_v_nrgs_001_03 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();
//

    emission_name = path_to_coeffs + "v/pair_nrgs_v_0-2.txt";
    in_file.open(emission_name.c_str());
    pair_v_nrgs_0_2 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "v/pair_nrgs_v_2-20.txt";
    in_file.open(emission_name.c_str());
    pair_v_nrgs_2_20 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "v/pair_nrgs_v_20-80.txt";
    in_file.open(emission_name.c_str());
    pair_v_nrgs_20_80 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "v/pair_nrgs_v_80-600.txt";
    in_file.open(emission_name.c_str());
    pair_v_nrgs_80_600 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "v/pair_nrgs_v_600-2000.txt";
    in_file.open(emission_name.c_str());
    pair_v_nrgs_600_2000 = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    //prepare photon emission nrg lookup table
    Chebyshev_User_2D *pair_v_nrgs_rates[5] = {&pair_v_nrgs_80_600,
                                            &pair_v_nrgs_20_80,
                                            &pair_v_nrgs_2_20,
                                            &pair_v_nrgs_0_2,
                                            &pair_v_nrgs_001_03};

    look_up_table_pair_v_nrgs[0] = &pair_v_nrgs_600_2000;                                            

    for(int i = 0; i < 5; i++){
        tmp_i = 1 << i;
        for(int j = 0; j < tmp_i; j++){
            look_up_table_pair_v_nrgs[tmp_i + j] = pair_v_nrgs_rates[i];
        }
    }

    // look_up_table_pair_v_nrgs[0] = &pair_v_nrgs_600_2000;  
    // look_up_table_pair_v_nrgs[1] = &pair_v_nrgs_80_600,
    // look_up_table_pair_v_nrgs[2] = &pair_v_nrgs_20_80,
    // look_up_table_pair_v_nrgs[3] = &pair_v_nrgs_2_20,
    // look_up_table_pair_v_nrgs[4] = &pair_v_nrgs_0_2,
    // look_up_table_pair_v_nrgs[5] = &pair_v_nrgs_001_03;

    //photon energy in v high part

//
    emission_name = path_to_coeffs + "v/pair_nrgs_v_high_001-03.txt";
    in_file.open(emission_name.c_str());
    pair_v_nrgs_001_03_high = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();
//

    emission_name = path_to_coeffs + "v/pair_nrgs_v_high_0-2.txt";
    in_file.open(emission_name.c_str());
    pair_v_nrgs_0_2_high = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "v/pair_nrgs_v_high_2-20.txt";
    in_file.open(emission_name.c_str());
    pair_v_nrgs_2_20_high = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "v/pair_nrgs_v_high_20-80.txt";
    in_file.open(emission_name.c_str());
    pair_v_nrgs_20_80_high = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "v/pair_nrgs_v_high_80-600.txt";
    in_file.open(emission_name.c_str());
    pair_v_nrgs_80_600_high = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "v/pair_nrgs_v_high_600-2000.txt";
    in_file.open(emission_name.c_str());
    pair_v_nrgs_600_2000_high = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    //prepare photon emission nrg lookup table
    Chebyshev_User_2D *pair_v_nrgs_rates_high[5] = {&pair_v_nrgs_80_600_high,
                                            &pair_v_nrgs_20_80_high,
                                            &pair_v_nrgs_2_20_high,
                                            &pair_v_nrgs_0_2_high,
                                            &pair_v_nrgs_001_03_high};

    look_up_table_pair_v_nrgs_high[0] = &pair_v_nrgs_600_2000_high;                                            

    for(int i = 0; i < 5; i++){
        tmp_i = 1 << i;
        for(int j = 0; j < tmp_i; j++){
            look_up_table_pair_v_nrgs_high[tmp_i + j] = pair_v_nrgs_rates_high[i];
        }
    }

    // look_up_table_pair_v_nrgs_high[0] = &pair_v_nrgs_600_2000_high;
    // look_up_table_pair_v_nrgs_high[1] = &pair_v_nrgs_80_600_high;
    // look_up_table_pair_v_nrgs_high[2] = &pair_v_nrgs_20_80_high;
    // look_up_table_pair_v_nrgs_high[3] = &pair_v_nrgs_2_20_high;
    // look_up_table_pair_v_nrgs_high[4] = &pair_v_nrgs_0_2_high;
    // look_up_table_pair_v_nrgs_high[5] = &pair_v_nrgs_001_03_high;

    //inverse 1D projection (to be fixed)

//
    pair_v_nrg_proj_001_03.a = pair_v_nrgs_001_03_high.a;
    pair_v_nrg_proj_001_03.b = pair_v_nrgs_001_03_high.b;
    pair_v_nrg_proj_001_03.domain_half_range = pair_v_nrgs_001_03_high.domain_half_range_N;
    pair_v_nrg_proj_001_03.domain_middle_point = pair_v_nrgs_001_03_high.domain_middle_point_N;
    pair_v_nrg_proj_001_03.evaluation_order = pair_v_nrgs_001_03_high.evaluation_order_N;
    pair_v_nrg_proj_001_03.last_coeffs = pair_v_nrgs_001_03_high.evaluate_y(pair_r_exp_limit_001_03);
//

    pair_v_nrg_proj_0_2.a = pair_v_nrgs_0_2_high.a;
    pair_v_nrg_proj_0_2.b = pair_v_nrgs_0_2_high.b;
    pair_v_nrg_proj_0_2.domain_half_range = pair_v_nrgs_0_2_high.domain_half_range_N;
    pair_v_nrg_proj_0_2.domain_middle_point = pair_v_nrgs_0_2_high.domain_middle_point_N;
    pair_v_nrg_proj_0_2.evaluation_order = pair_v_nrgs_0_2_high.evaluation_order_N;
    pair_v_nrg_proj_0_2.last_coeffs = pair_v_nrgs_0_2_high.evaluate_y(pair_r_exp_limit_0_2);

    pair_v_nrg_proj_2_20.a = pair_v_nrgs_2_20_high.a;
    pair_v_nrg_proj_2_20.b = pair_v_nrgs_2_20_high.b;
    pair_v_nrg_proj_2_20.domain_half_range = pair_v_nrgs_2_20_high.domain_half_range_N;
    pair_v_nrg_proj_2_20.domain_middle_point = pair_v_nrgs_2_20_high.domain_middle_point_N;
    pair_v_nrg_proj_2_20.evaluation_order = pair_v_nrgs_2_20_high.evaluation_order_N;
    pair_v_nrg_proj_2_20.last_coeffs = pair_v_nrgs_2_20_high.evaluate_y(pair_r_exp_limit_2_20);

    pair_v_nrg_proj_20_80.a = pair_v_nrgs_20_80_high.a;
    pair_v_nrg_proj_20_80.b = pair_v_nrgs_20_80_high.b;
    pair_v_nrg_proj_20_80.domain_half_range = pair_v_nrgs_20_80_high.domain_half_range_N;
    pair_v_nrg_proj_20_80.domain_middle_point = pair_v_nrgs_20_80_high.domain_middle_point_N;
    pair_v_nrg_proj_20_80.evaluation_order = pair_v_nrgs_20_80_high.evaluation_order_N;
    pair_v_nrg_proj_20_80.last_coeffs = pair_v_nrgs_20_80_high.evaluate_y(pair_r_exp_limit_20_80);

    pair_v_nrg_proj_80_600.a = pair_v_nrgs_80_600_high.a;
    pair_v_nrg_proj_80_600.b = pair_v_nrgs_80_600_high.b;
    pair_v_nrg_proj_80_600.domain_half_range = pair_v_nrgs_80_600_high.domain_half_range_N;
    pair_v_nrg_proj_80_600.domain_middle_point = pair_v_nrgs_80_600_high.domain_middle_point_N;
    pair_v_nrg_proj_80_600.evaluation_order = pair_v_nrgs_80_600_high.evaluation_order_N;
    pair_v_nrg_proj_80_600.last_coeffs = pair_v_nrgs_80_600_high.evaluate_y(pair_r_exp_limit_80_600);

    pair_v_nrg_proj_600_2000.a = pair_v_nrgs_600_2000_high.a;
    pair_v_nrg_proj_600_2000.b = pair_v_nrgs_600_2000_high.b;
    pair_v_nrg_proj_600_2000.domain_half_range = pair_v_nrgs_600_2000_high.domain_half_range_N;
    pair_v_nrg_proj_600_2000.domain_middle_point = pair_v_nrgs_600_2000_high.domain_middle_point_N;
    pair_v_nrg_proj_600_2000.evaluation_order = pair_v_nrgs_600_2000_high.evaluation_order_N;
    pair_v_nrg_proj_600_2000.last_coeffs = pair_v_nrgs_600_2000_high.evaluate_y(pair_r_exp_limit_600_2000);

    Chebyshev_User_1D *pair_projs[5] = {&pair_v_nrg_proj_80_600,
                                            &pair_v_nrg_proj_20_80,
                                            &pair_v_nrg_proj_2_20,
                                            &pair_v_nrg_proj_0_2,
                                            &pair_v_nrg_proj_001_03};

    look_up_table_pair_v_nrgs_proj[0] = &pair_v_nrg_proj_600_2000;                                            

    for(int i = 0; i < 5; i++){
        tmp_i = 1 << i;
        for(int j = 0; j < tmp_i; j++){
            look_up_table_pair_v_nrgs_proj[tmp_i + j] = pair_projs[i];
        }
    }

    // look_up_table_pair_v_nrgs_proj[0] = &pair_v_nrg_proj_600_2000; 
    // look_up_table_pair_v_nrgs_proj[1] = &pair_v_nrg_proj_80_600;
    // look_up_table_pair_v_nrgs_proj[2] = &pair_v_nrg_proj_20_80;
    // look_up_table_pair_v_nrgs_proj[3] = &pair_v_nrg_proj_2_20;
    // look_up_table_pair_v_nrgs_proj[4] = &pair_v_nrg_proj_0_2;
    // look_up_table_pair_v_nrgs_proj[5] = &pair_v_nrg_proj_001_03;
}

double SFQED_Processes::SFQED_PAIR_creation_rate(const double &gamma, const double &chi) const{
    
    //coefficient for the rate of emission
    //multiply by \tilde{W_rad} to get the rate of emission
    double coefW_pair = coef_rate_Wpair / gamma;
    
    //low impact Branch to calculate \tilde{W_rad} and then the rate
    //---------------------------------------------------
    /**/
    //lookup tables boolean selectors
    bool chi_001_03 = chi <= bound_kappa_0th, 
        chi_0_2 = chi <= bound_kappa_1st,
        chi_2_20 = chi <= bound_kappa_2nd,
        chi_20_80 = chi <= bound_kappa_3rd,
        chi_80_600 = chi <= bound_kappa_4th;
        //chi_600_2000 = chi <= 2000.;

    //IMPORTANT: through this constructs we mean
    //to default every chi > 600 case into the
    //600 < chi <= 2000 one
    const Chebyshev_User_1D * const chebyshev_pair_rate = 
                        (look_up_table_pair_prd_rt[chi_001_03*16 + chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600]);
    
    return coefW_pair * chebyshev_pair_rate->evaluate(chi);//
}

double SFQED_Processes::SFQED_PAIR_creation_rate_fast(const double &gamma, const double &chi) const{
    
    //coefficient for the rate of emission
    //multiply by \tilde{W_rad} to get the rate of emission
    double coefW_pair = coef_rate_Wpair / gamma;
    
    //low impact Branch to calculate \tilde{W_rad} and then the rate
    //---------------------------------------------------
    /**/
    //lookup tables boolean selectors
    bool chi_001_03 = chi <= bound_kappa_0th, 
        chi_0_2 = chi <= bound_kappa_1st,
        chi_2_20 = chi <= bound_kappa_2nd,
        chi_20_80 = chi <= bound_kappa_3rd,
        chi_80_600 = chi <= bound_kappa_4th;
        //chi_600_2000 = chi <= 2000.;

    //IMPORTANT: through this constructs we mean
    //to default every chi > 600 case into the
    //600 < chi <= 2000 one
    const Chebyshev_User_1D * const chebyshev_pair_rate = 
                        (look_up_table_pair_prd_rt_fast[chi_001_03*16 + chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600]);
    
    return coefW_pair * chebyshev_pair_rate->evaluate(chi);//
}

double SFQED_Processes::SFQED_PAIR_emitted_electron_energy_aux(const int &lookup_index, const double &chi, const double &rescaled_rnd) const{

    const double low_r_limit = *(lu_table_pair_lower_r_bounds[lookup_index]),
            high_r_limit = *(lu_table_pair_upper_r_bounds[lookup_index]);

    // std::cout << "low limit = " << low_r_limit << " and high limit = " << high_r_limit << '\n';

    // std::cout << "1\n";

    double v;

    if(rescaled_rnd < low_r_limit){

        // std::cout << "3\n";

        const Chebyshev_User_2D * const chebyshev_pair_v_nrgs_inv =
                (look_up_table_pair_v_nrgs[lookup_index]);

        // std::cout << "4\n";

        v = chebyshev_pair_v_nrgs_inv->evaluate(chi, rescaled_rnd);

        // std::cout << "6\n";
    }
    else if(rescaled_rnd < high_r_limit){

        // std::cout << "3\n";

        const Chebyshev_User_2D * const chebyshev_pair_v_nrgs_inv_high =
                (look_up_table_pair_v_nrgs_high[lookup_index]);

        // std::cout << "order in N = " << chebyshev_pair_v_nrgs_inv_high->evaluation_order_N <<
        //             "order in M = " << chebyshev_pair_v_nrgs_inv_high->evaluation_order_M << '\n';

        v = chebyshev_pair_v_nrgs_inv_high->evaluate(chi, rescaled_rnd);

        // std::cout << "6\n";
    }
    //in the off chance case that we happen to be beyond the
    //high r threshold, we will assume to be dealing directly
    //with v = 0. Since v = (\epsilon_{\gamma} - 2 \epsilon_e) / \epsilon_{\gamma}
    else{

        // std::cout << "7\n";

        //retrieve the appropriate pair rate set of coefficients
        const Chebyshev_User_1D * const chebyshev_pair_rate =
                (look_up_table_pair_prd_rt[lookup_index]);

        // std::cout << "8\n";

        //compute the rate of pair emission
        double tildeWrad = chebyshev_pair_rate->evaluate(chi) / 3.;
                            
        double rnd_tildeWrad = rescaled_rnd * tildeWrad;

        // std::cout << "9\n";
        
        //compute the pivot at which the exponential approximation will be applied
        const Chebyshev_User_1D * const chebyshev_pair_v_nrgs_proj =
                (look_up_table_pair_v_nrgs_proj[lookup_index]);

        // std::cout << "10\n";

        double v0 = chebyshev_pair_v_nrgs_proj->evaluate(chi);
        v = 8.0 / (3.0 * chi * (1.0 - v0*v0));

        // std::cout << "11\n";

        //integrate the pair creation differential probability
        //from the just computed pivot up to 1
        const Chebyshev_User_2D * const chebyshev_pair_prtl_rate =
                (look_up_table_pair_prtl_rt[lookup_index]);

        // std::cout << "12\n";

        double tail_integral = chebyshev_pair_prtl_rate->evaluate(chi, v0);

        // std::cout << "13\n";

        //now invert the exponential approximation and extract the v
        v = v - log((tildeWrad - rnd_tildeWrad) / (tildeWrad - tail_integral));
        v = sqrt(1.0 - 8.0 / (3.0 * chi * v));

        // std::cout << "14\n";
    }

    return v;
}

//this method returns the nrg of the emitted electron only!
//please consider that the relation phtn_nrg = e_nrg + p_nrg
//holds if you are interested in the positron energy
double SFQED_Processes::SFQED_PAIR_created_electron_energy(const double &nrg, const double &chi, const double &rnd) const{

    //lookup tables boolean selectors
    bool chi_001_03 = chi <= bound_kappa_0th, 
        chi_0_2 = chi <= bound_kappa_1st,
        chi_2_20 = chi <= bound_kappa_2nd,
        chi_20_80 = chi <= bound_kappa_3rd,
        chi_80_600 = chi <= bound_kappa_4th;
        //chi_600_2000 = chi <= 2000.;

    int lookup_index = chi_001_03*16 + chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600;
    //lookup_index = (1 - std::signbit(lookup_index-1))*(32 - __builtin_clz(lookup_index));

    double rescaled_rnd = 2. * rnd - 1.;
    double sgn = -2. * std::signbit(rescaled_rnd) + 1.;
    //extract the absolute value
    rescaled_rnd *= sgn;

    //v to return
    double v = SFQED_PAIR_emitted_electron_energy_aux(lookup_index, chi, rescaled_rnd);

    //return the electron energy
    return nrg * 0.5 * (1. + sgn * v);
}

double SFQED_Processes::SFQED_PAIR_created_electron_energy_fast(const double &nrg, const double &chi, const double &rnd) const{

    //lookup tables boolean selectors
    bool chi_0_2 = chi <= bound_kappa_1st,
        chi_2_20 = chi <= bound_kappa_2nd,
        chi_20_80 = chi <= bound_kappa_3rd,
        chi_80_600 = chi <= bound_kappa_4th;
        //chi_600_2000 = chi <= 2000.;

    int lookup_index = chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600;
    //lookup_index = (1 - std::signbit(lookup_index-1))*(32 - __builtin_clz(lookup_index));

    double rescaled_rnd = 2. * rnd - 1.;
    double sgn = -2. * std::signbit(rescaled_rnd) + 1.;
    //extract the absolute value
    rescaled_rnd *= sgn;

    //v to return
    double v = SFQED_PAIR_emitted_electron_energy_aux(lookup_index, chi, rescaled_rnd);

    //return the electron energy
    return nrg * 0.5 * (1. + sgn * v);
} 


void SFQED_Processes::SFQED_finalize_PAIR_creation(){
    //purge emission rates
    delete[] pair_prd_rt_0_2.last_coeffs;
    delete[] pair_prd_rt_2_20.last_coeffs;
    delete[] pair_prd_rt_20_80.last_coeffs;
    delete[] pair_prd_rt_80_600.last_coeffs;
    delete[] pair_prd_rt_600_2000.last_coeffs;

    delete[] pair_prtl_rt_0_2.last_coeffs;
    delete[] pair_prtl_rt_2_20.last_coeffs;
    delete[] pair_prtl_rt_20_80.last_coeffs;
    delete[] pair_prtl_rt_80_600.last_coeffs;
    delete[] pair_prtl_rt_600_2000.last_coeffs;

    delete[] pair_v_nrgs_001_03.last_coeffs;
    delete[] pair_v_nrgs_0_2.last_coeffs;
    delete[] pair_v_nrgs_2_20.last_coeffs;
    delete[] pair_v_nrgs_20_80.last_coeffs;
    delete[] pair_v_nrgs_80_600.last_coeffs;
    delete[] pair_v_nrgs_600_2000.last_coeffs;

    delete[] pair_v_nrgs_001_03_high.last_coeffs;
    delete[] pair_v_nrgs_0_2_high.last_coeffs;
    delete[] pair_v_nrgs_2_20_high.last_coeffs;
    delete[] pair_v_nrgs_20_80_high.last_coeffs;
    delete[] pair_v_nrgs_80_600_high.last_coeffs;
    delete[] pair_v_nrgs_600_2000_high.last_coeffs;

    delete[] pair_v_nrg_proj_001_03.last_coeffs;
    delete[] pair_v_nrg_proj_0_2.last_coeffs;
    delete[] pair_v_nrg_proj_2_20.last_coeffs;
    delete[] pair_v_nrg_proj_20_80.last_coeffs;
    delete[] pair_v_nrg_proj_80_600.last_coeffs;
    delete[] pair_v_nrg_proj_600_2000.last_coeffs;
}


/////////////////////
// DEBUG FUNCTIONS //
/////////////////////

// double SFQED_Processes::SFQED_BLCFA_emitted_photon_energy_derivative(const double& LCFA_limit,
//                                                                         const double& gamma, 
//                                                                         const double& chi, 
//                                                                         const double& rnd, 
//                                                                         const double& rnd2){

//     static double increment = 0.001;

//     //key generating section
//     bool chi_0_2 = chi <= bound_chi_1st,
//         chi_2_20 = chi <= bound_chi_2nd,
//         chi_20_80 = chi <= bound_chi_3rd,
//         chi_80_600 = chi <= bound_chi_4th;
//         //chi_600_2000 = chi <= 2000.;
//     int lookup_index = chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600;

//     //transform the threshold energy into a w value
//     double LCFA_limit_w  = cbrt(LCFA_limit / (1.5 * chi * (gamma - LCFA_limit)));

//     //in case the phtn nrg threshold is very close to the particle nrg
//     //or it is above the upper bound of the approximated domain
//     //we don't want the emission to occur 
//     if(LCFA_limit > 0.75 * gamma || LCFA_limit_w > *(lu_table_upper_w_bounds[lookup_index])){
//         return 0.;
//     }

//     //compute the w-value energy of the LCFA emitted photon
//     double emitted_phtn_LCFA_nrg_w = SFQED_LCFA_phtn_nrg_aux(chi, rnd, lookup_index),
//             rs_dw_over_depsgamma;

//     //convert the phtn LCFA energy limit (LCFA_limit) to a w value
//     //maybe it would be smart doing this before entering this function
//     // double LCFA_limit_w  = cbrt(LCFA_limit / (1.5 * chi * (gamma - LCFA_limit)));

//     /////////////////////////////////////////
//     //  CUMULATIVE WITH DERIVATIVE METHOD 

//     // double w_validity = *(lu_table_upper_w_bounds[lookup_index]);
//     // //if the w value of the energy threshold is above the limit
//     // //of validity the emission is suppressed by default
//     // if(LCFA_limit_w > w_validity){
//     //     return 0.;
//     // }

//     //BLCFA section (it starts with a comparison)
//     if(emitted_phtn_LCFA_nrg_w < LCFA_limit_w){

//         //actully it is the integral
//         Chebyshev_User_2D phtn_diff_crss_sctn =
//                 *(look_up_table_phtn_diff_crss_sctn[lookup_index]);

//         Chebyshev_User_1D collapsed_x(phtn_diff_crss_sctn.evaluation_order_M, phtn_diff_crss_sctn.c, phtn_diff_crss_sctn.d);
//         collapsed_x.last_coeffs = phtn_diff_crss_sctn.evaluate_x(chi);

//         double f_w_up, f_w;
//         double dWrad_en_cut, dWrad_ph_en;

//         //**********work with the threshold*****************

//         f_w_up = collapsed_x.evaluate(LCFA_limit_w + increment);
//         f_w = collapsed_x.evaluate(LCFA_limit_w);
//         rs_dw_over_depsgamma = (1. + 1.5 * chi * LCFA_limit_w * LCFA_limit_w * LCFA_limit_w) / LCFA_limit_w;
        
//         dWrad_en_cut = (f_w_up - f_w) / (increment) * rs_dw_over_depsgamma * rs_dw_over_depsgamma; // / determinant(LCFA_limit_w, chi);

//         //***********work with the actual LCFA nrg******************
          
//         f_w_up = collapsed_x.evaluate(emitted_phtn_LCFA_nrg_w + increment);
//         f_w = collapsed_x.evaluate(emitted_phtn_LCFA_nrg_w);
//         rs_dw_over_depsgamma = (1. + 1.5 * chi * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w) / emitted_phtn_LCFA_nrg_w;
          
//         dWrad_ph_en = (f_w_up - f_w) / (increment) * rs_dw_over_depsgamma * rs_dw_over_depsgamma; // / determinant(emitted_phtn_LCFA_nrg_w, chi);

//         //now I use the rejection sampling method (employing rnd2) to determine if I have
//         //to keep or reject the just occurred emission
//         if((rnd2 * dWrad_ph_en) > dWrad_en_cut){
//             return 0.;
//         }

//     }

//     emitted_phtn_LCFA_nrg_w = 3.0 * chi * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w;
//     emitted_phtn_LCFA_nrg_w = gamma * emitted_phtn_LCFA_nrg_w / (2.0 + emitted_phtn_LCFA_nrg_w);


//     return emitted_phtn_LCFA_nrg_w;
//     ///////////////////////////////////////////////
// }

// double SFQED_Processes::SFQED_BLCFA_emitted_photon_energy_no_approx(const double& LCFA_limit,
//                                                                     const double& gamma, 
//                                                                     const double& chi, 
//                                                                     const double& rnd, 
//                                                                     const double& rnd2){

//     //key generating section
//     bool chi_0_2 = chi <= bound_chi_1st,
//         chi_2_20 = chi <= bound_chi_2nd,
//         chi_20_80 = chi <= bound_chi_3rd,
//         chi_80_600 = chi <= bound_chi_4th;
//         //chi_600_2000 = chi <= 2000.;
//     int lookup_index = chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600;
    
//     //in case the phtn nrg threshold is very close to the particle nrg
//     //or it is above the upper bound of the approximated domain
//     //we don't want the emission to occur 
//     if(LCFA_limit > 0.75 * gamma){
//         return 0.;
//     }


//     //compute the w-value energy of the LCFA threshold and emitted photon
//     double LCFA_limit_w  = cbrt(LCFA_limit / (1.5 * chi * (gamma - LCFA_limit))),
//             emitted_phtn_LCFA_nrg_w = SFQED_LCFA_phtn_nrg_aux(chi, rnd, lookup_index),
//             rs_dw_over_depsgamma;

//     /////////////////////////////////////////
//     //  NO APPROX METHOD

//     // double LCFA_limit_w = cbrt((LCFA_limit) / (1.5 * chi * (gamma - LCFA_limit)));

//     //BLCFA section (it starts with a comparison)
//     if(emitted_phtn_LCFA_nrg_w < LCFA_limit_w){

//         double dWrad_en_cut, dWrad_ph_en;

//         //**********work with the threshold*****************
        
//         dWrad_en_cut = differential_cross_section_phtn_emission(chi, LCFA_limit_w, NULL);
//         rs_dw_over_depsgamma = (1. + 1.5 * chi * LCFA_limit_w * LCFA_limit_w * LCFA_limit_w) / LCFA_limit_w;
//         dWrad_en_cut *= rs_dw_over_depsgamma * rs_dw_over_depsgamma;// /= determinant(LCFA_limit_w, chi);

//         //***********work with the actual LCFA nrg******************
          
//         dWrad_ph_en = differential_cross_section_phtn_emission(chi, emitted_phtn_LCFA_nrg_w, NULL);
//         rs_dw_over_depsgamma = (1. + 1.5 * chi * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w) / emitted_phtn_LCFA_nrg_w;
//         dWrad_ph_en *= rs_dw_over_depsgamma * rs_dw_over_depsgamma;// /= determinant(emitted_phtn_LCFA_nrg_w, chi);

//         //now I use the rejection sampling method (employing rnd2) to determine if I have
//         //to keep or reject the just occurred emission
//         if((rnd2 * dWrad_ph_en) > dWrad_en_cut){
//             return 0.;
//         }
//     }

//     emitted_phtn_LCFA_nrg_w = 3.0 * chi * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w;
//     emitted_phtn_LCFA_nrg_w = gamma * emitted_phtn_LCFA_nrg_w / (2.0 + emitted_phtn_LCFA_nrg_w);

//     return emitted_phtn_LCFA_nrg_w;
//     ///////////////////////////////////////////////
// }

// double SFQED_Processes::SFQED_BLCFA_matteo(const double& threshold,
//                                             const double& gam,
//                                             const double& chi, 
//                                             const double& rnd, 
//                                             const double& rnd2){
    
//     static double b_1st = 2.;
//     // static double b_2nd = 20.;
//     // static double b_3rd = 80.;
    
//     static double b_w_1st = 1.6;
//     // static double b_w_2nd = 1.33;
//     // static double b_w_3rd = 1.12;

//     static double delta = 1.0e-2;

//     //method 1
//     // double LCFA_nrg = this->SFQED_LCFA_emitted_photon_energy(gam, chi, rnd);
//     /////////////////////////////////

//     //method 2
//     //key generating section
//     bool chi_0_2 = chi <= bound_chi_1st,
//         chi_2_20 = chi <= bound_chi_2nd,
//         chi_20_80 = chi <= bound_chi_3rd,
//         chi_80_600 = chi <= bound_chi_4th;
//         //chi_600_2000 = chi <= 2000.;
//     int lookup_index = chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600;

//     //compute the w-value energy of the LCFA emitted photon
//     double LCFA_nrg = SFQED_LCFA_phtn_nrg_aux(chi, rnd, lookup_index);
//     LCFA_nrg = 3.0 * chi * LCFA_nrg * LCFA_nrg * LCFA_nrg;
//     LCFA_nrg = gam * LCFA_nrg / (2.0 + LCFA_nrg);
//     //////////////////////////////////

//     double en_up, v_up, w_up, en, v, w;
//     double normal_chi, normal_w_up, normal_w, f_w_up, f_w;
//     double dW_dw, dw_deg, dWrad_en_cut, dWrad_ph_en, constant;

//     double* table;

//     if(LCFA_nrg < threshold){

//         //project the 2d cheby coeffs onto a 1d array (using the chi value)
//         normal_chi =  scale_interval_to_normal(0., b_1st, chi);
//         //remember to purge this
//         table = clenshaw_2d_step1(intdWrad_0chi2, dim_N_0chi2, dim_M_0chi2, normal_chi);

//         //**********work with the threshold*****************
//         en_up = threshold + delta;
//         v_up = (2.0 * en_up) / (3.0 * chi * (gam - en_up));
//         w_up = cbrt(v_up);
        
//         en = threshold;
//         v = (2.0 * en) / (3.0 * chi * (gam - en));
//         w = cbrt(v);
        
//         normal_w_up = scale_interval_to_normal(0., b_w_1st, w_up);
//         normal_w = scale_interval_to_normal(0., b_w_1st, w);

//         f_w_up = clenshaw_1d(table, 0., dim_N_0chi2, normal_w_up);
//         f_w = clenshaw_1d(table, 0., dim_N_0chi2, normal_w);
  
//         dW_dw = (f_w_up - f_w) / (w_up - w);
//         dw_deg = 2.0 / (81.0 * chi * en * en * (gam - en) * (gam - en) * (gam - en) * (gam - en));
//         dw_deg = cbrt(dw_deg) * gam;
//         dWrad_en_cut = dW_dw * dw_deg;

//         //***********work with the actual LCFA nrg******************
//         en_up = LCFA_nrg + delta;
//         v_up = (2.0 * en_up) / (3.0 * chi * (gam - en_up));
//         w_up = cbrt(v_up);
          
//         en = LCFA_nrg;
//         v = (2.0 * en) / (3.0 * chi * (gam - en));
//         w = cbrt(v);
          
//         normal_w_up = scale_interval_to_normal(0.0, b_w_1st, w_up);
//         normal_w = scale_interval_to_normal(0.0, b_w_1st, w);
          
//         f_w_up = clenshaw_1d(table, 0, dim_N_0chi2, normal_w_up);
//         f_w = clenshaw_1d(table, 0, dim_N_0chi2, normal_w);
          
//         dW_dw = (f_w_up - f_w) / (w_up - w);
//         dw_deg = 2.0 / (81.0 * chi * en * en * (gam - en) * (gam - en) * (gam - en) * (gam - en));
//         dw_deg = cbrt(dw_deg) * gam;
//         dWrad_ph_en = dW_dw * dw_deg;

//         //purge the table
//         delete[] table;

//         if((rnd2 * dWrad_ph_en) > dWrad_en_cut){
//             return 0.;
//         }
//     }

//     return LCFA_nrg;
// }

// double SFQED_Processes::SFQED_BLCFA_DEBUG(const double& LCFA_limit,
//                             const double& gamma,
//                             const double& chi, 
//                             const double& rnd, 
//                             const double& rnd2,
//                             std::ofstream& out_file){

//     static double increment = 0.01;

//     //key generating section
//     bool chi_0_2 = chi <= bound_chi_1st,
//         chi_2_20 = chi <= bound_chi_2nd,
//         chi_20_80 = chi <= bound_chi_3rd,
//         chi_80_600 = chi <= bound_chi_4th;
//         //chi_600_2000 = chi <= 2000.;
//     int lookup_index = chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600;

//     //compute the w-value energy of the LCFA emitted photon
//     double emitted_phtn_LCFA_nrg_w = SFQED_LCFA_phtn_nrg_aux(chi, rnd, lookup_index),
//             diff_probability;

//     //convert the phtn LCFA energy limit (LCFA_limit) to a w value
//     double LCFA_limit_w  = cbrt(LCFA_limit / (1.5 * chi * (gamma - LCFA_limit))),
//             diff_prob_LCFA;

//     double emitted_phtn_LCFA_nrg = 3.0 * chi * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w;
//     emitted_phtn_LCFA_nrg = gamma * emitted_phtn_LCFA_nrg / (2.0 + emitted_phtn_LCFA_nrg);

//     double* table;

//     //BLCFA section (it starts with a comparison)
//     if(emitted_phtn_LCFA_nrg < LCFA_limit){


//         Chebyshev_User_2D phtn_diff_crss_sctn =
//                 *(look_up_table_phtn_diff_crss_sctn[lookup_index]);

//         double en_up, v_up, w_up, en, v, w;
//         double normal_chi, normal_w_up, normal_w, f_w_up, f_w;
//         double dW_dw, dw_deg, dWrad_en_cut, dWrad_en_cut_1, dWrad_en_cut_2, dWrad_ph_en, dWrad_ph_en_1, dWrad_ph_en_2;

//         //project the 2d cheby coeffs onto a 1d array (using the chi value)
//         normal_chi =  scale_interval_to_normal(0., 2., chi);
//         //remember to purge this
//         table = clenshaw_2d_step1(intdWrad_0chi2, dim_N_0chi2, dim_M_0chi2, normal_chi);

//         //**********work with the threshold*****************

//         //w value from gamma
//         en = LCFA_limit;
//         v = (2.0 * en) / (3.0 * chi * (gamma - en));
//         //change w (use the one above directly)
//         w = cbrt(v);
        
//         //determinant
//         //change determinant (use the function)
//         dw_deg = 2.0 / (81.0 * chi * en * en * (gamma - en) * (gamma - en) * (gamma - en) * (gamma - en));
//         dw_deg = cbrt(dw_deg) * gamma;

//         // out_file << chi << ' ' << w << ' ';
//         out_file << w << ' ';

//         //MY//////////////

//         //change function (use your approximation)
//         dW_dw = phtn_diff_crss_sctn.evaluate(chi, w);
//         dWrad_en_cut = dW_dw * dw_deg;

//         out_file << dWrad_en_cut << ' ';

//         //MATTEO//////////
//         en_up = LCFA_limit + increment;
//         v_up = (2.0 * en_up) / (3.0 * chi * (gamma - en_up));
//         w_up = cbrt(v_up);
        
//         normal_w_up = scale_interval_to_normal(0., 1.6, w_up);
//         normal_w = scale_interval_to_normal(0., 1.6, w);

//         f_w_up = clenshaw_1d(table, 0., dim_N_0chi2, normal_w_up);
//         f_w = clenshaw_1d(table, 0., dim_N_0chi2, normal_w);
  
//         dW_dw = (f_w_up - f_w) / (w_up - w);
//         dWrad_en_cut_1 = dW_dw * dw_deg;

//         out_file << dWrad_en_cut_1 << ' ';

//         //ORIGINAL
//         dW_dw = differential_cross_section_phtn_emission(chi, w, NULL);
//         dWrad_en_cut_2 = dW_dw * dw_deg;
        
//         out_file << dWrad_en_cut_2 << ' ';


//         //***********work with the actual LCFA nrg******************

//         //w value from gamma
//         en = emitted_phtn_LCFA_nrg;
//         v = (2.0 * en) / (3.0 * chi * (gamma - en));
//         w = cbrt(v);

//         //determinant
//         dw_deg = 2.0 / (81.0 * chi * en * en * (gamma - en) * (gamma - en) * (gamma - en) * (gamma - en));
//         dw_deg = cbrt(dw_deg) * gamma;

//         out_file << w << ' ';
          
//         //MY  
//         dW_dw = phtn_diff_crss_sctn.evaluate(chi, w);
//         dWrad_ph_en = dW_dw * dw_deg;

//         out_file << dWrad_ph_en << ' ';

//         if((rnd2 * dWrad_ph_en) > dWrad_en_cut){
//             out_file << 0. << ' ';
//         }
//         else{
//             out_file << emitted_phtn_LCFA_nrg << ' ';
//         }

//         //MATTEO
//         en_up = emitted_phtn_LCFA_nrg + increment;
//         v_up = (2.0 * en_up) / (3.0 * chi * (gamma - en_up));
//         w_up = cbrt(v_up);
          
//         normal_w_up = scale_interval_to_normal(0.0, 1.6, w_up);
//         normal_w = scale_interval_to_normal(0.0, 1.6, w);
          
//         f_w_up = clenshaw_1d(table, 0, dim_N_0chi2, normal_w_up);
//         f_w = clenshaw_1d(table, 0, dim_N_0chi2, normal_w);
          
//         dW_dw = (f_w_up - f_w) / (w_up - w);
//         dWrad_ph_en_1 = dW_dw * dw_deg;

//         out_file << dWrad_ph_en_1 << ' ';

//         //purge the table
//         delete[] table;

//         if((rnd2 * dWrad_ph_en_1) > dWrad_en_cut_1){
//             out_file << 0. << ' ';
//         }
//         else{
//             out_file << emitted_phtn_LCFA_nrg << ' ';
//         }

//         //ORIGINAL
//         dW_dw = differential_cross_section_phtn_emission(chi, w, NULL);
//         dWrad_ph_en_2 = dW_dw * dw_deg;
        
//         out_file << dWrad_ph_en_2 << ' ';

//         if((rnd2 * dWrad_ph_en_2) > dWrad_en_cut_2){
//             out_file << 0. << '\n';
//             return 0.;
//         }
//         else{
//             out_file << emitted_phtn_LCFA_nrg << '\n';
//         }

//     }

//     return emitted_phtn_LCFA_nrg;

// };

// //read chebyshev coefficients from simple txt file
// double* load_cheb_coefficients(std::ifstream& file, int ubound1, int ubound2){

//     short n_elements = 0;
//     double tmp_coefficient;

//     Coeff_Cheby* list_head = NULL;
//     Coeff_Cheby* list_tail = NULL;
//     Coeff_Cheby* list_element = NULL;

//     double* data_to_return = new double[ubound1*ubound2];

//     //read elements from file
//     while(true){
//         file >> tmp_coefficient;

//         if(file.eof()){
//             break;
//         }

//         n_elements++;

//         list_element = new Coeff_Cheby(tmp_coefficient);

//         if(!list_head){
//             list_head = list_element;
//             list_tail = list_element;
//         }
//         else{
//             list_tail->pointer_to_next = list_element;
//             list_tail = list_element;
//         }
//     }
    
//     for(int i = 0; i < n_elements; i++){
//         list_element = list_head;
//         list_head = list_element->pointer_to_next;
//         data_to_return[i] = list_element->coeff;
//         delete list_element;
//     }

//     return data_to_return;

// }

// double* clenshaw_2d_step1(const double* const coefficients, const int& N, const int& M, const double& evaluation_point){
//     double tmp_up, tmp_down;

//     double* data_to_return = new double[N];

//     for(int j = 0; j < N; j++){
//         tmp_up = coefficients[j * (M)];
//         tmp_down = coefficients[j * (M) + 1] + 2. * evaluation_point * tmp_up;

//         for(int i = 2; i < M; i+=2){
//             tmp_up = coefficients[i + j*M] - tmp_up + 2. * evaluation_point * tmp_down;
//             tmp_down = coefficients[i + 1 + j*M] - tmp_down + 2. * evaluation_point * tmp_up;
//         }

//         data_to_return[N - 1 - j] = tmp_down - evaluation_point * tmp_up;
//     }

//     return data_to_return;
    
// }

// double clenshaw_1d(const double* const coeffs, const int& start, const int& size, const double& evaluation_point){

//     double result, up, down;

//     up = coeffs[start];
//     down = coeffs[start + 1] + 2. * evaluation_point * up;
//     for(int j = start + 2; j < size; j+=2){
//         up = coeffs[j] - up + 2. * evaluation_point * down;
//         down = coeffs[j + 1] - down + 2. * evaluation_point * up;
//     }
    
//     result = down - evaluation_point * up;

//     return result;

// }

// void SFQED_Processes::initialize_intdWrad(std::string path_to_coeffs){
//     std::string emission_name = path_to_coeffs + "intdWrad_0chi2.dat";
//     std::ifstream in_file(emission_name.c_str());
//     intdWrad_0chi2 = load_cheb_coefficients(in_file, dim_N_0chi2, dim_M_0chi2);
//     in_file.close();
//     in_file.clear();
// }

// double scale_interval_to_normal(const double& a, const double& b, const double& x){
//     return (2. * x - a - b) / (b - a);
// }