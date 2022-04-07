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

void SFQED_Processes::SFQED_set_all_to_one(){

    omega_r = 1;
    
    //reference angular frequency = c / Lambda
    Lambda = 1;

    //normalized compton time
    norm_Compton_time = 1;

    norm_Compton_time_2 = 1;

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

void SFQED_Processes::SFQED_set_time_step(const double& dt){
    one_over_dt = 1. / dt;
    one_over_dt2 = one_over_dt * one_over_dt;
}

/*****************/
/*PHOTON EMISSION*/
/*****************/


void SFQED_Processes::SFQED_init_PHTN_emission(std::string path_to_coeffs){

    init_phtn_mx_tables();

    int tmp_i;

    //photon emission rate
    std::string emission_name = path_to_coeffs + "full_rate/phtn_emission_rate_0-2.txt";
    std::ifstream in_file(emission_name.c_str());
    look_up_table_phtn_mx_rt[4] = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "full_rate/phtn_emission_rate_2-20.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_mx_rt[3] = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "full_rate/phtn_emission_rate_20-80.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_mx_rt[2] = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "full_rate/phtn_emission_rate_80-600.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_mx_rt[1] = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "full_rate/phtn_emission_rate_600-2000.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_mx_rt[0] = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();


    //photon partial emission rate g(w,\chi) = \int_0^w{Wrad(\chi,w_{integrated})dw_{integrated}}
    emission_name = path_to_coeffs + "prtl_rate/phtn_partial_rate_0-2.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_prtl_rt[4] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "prtl_rate/phtn_partial_rate_2-20.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_prtl_rt[3] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "prtl_rate/phtn_partial_rate_20-80.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_prtl_rt[2] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "prtl_rate/phtn_partial_rate_80-600.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_prtl_rt[1] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "prtl_rate/phtn_partial_rate_600-2000.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_prtl_rt[0] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();


    //photon energy in w
    emission_name = path_to_coeffs + "w/phtn_w_nrg_0-2.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_momw[4] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "w/phtn_w_nrg_2-20.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_momw[3] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "w/phtn_w_nrg_20-80.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_momw[2] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "w/phtn_w_nrg_80-600.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_momw[1] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "w/phtn_w_nrg_600-2000.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_momw[0] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    
    //photon energy in 1/w
    emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_0-2.txt";
    in_file.open(emission_name.c_str());
    look_up_table_1_over_w[4] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_2-20.txt";
    in_file.open(emission_name.c_str());
    look_up_table_1_over_w[3] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_20-80.txt";
    in_file.open(emission_name.c_str());
    look_up_table_1_over_w[2] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_80-600.txt";
    in_file.open(emission_name.c_str());
    look_up_table_1_over_w[1] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_600-2000.txt";
    in_file.open(emission_name.c_str());
    look_up_table_1_over_w[0] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();


    look_up_table_phtn_1_o_w_proj[4].a = look_up_table_1_over_w[4].a;
    look_up_table_phtn_1_o_w_proj[4].b = look_up_table_1_over_w[4].b;
    look_up_table_phtn_1_o_w_proj[4].domain_half_range = look_up_table_1_over_w[4].domain_half_range_N;
    look_up_table_phtn_1_o_w_proj[4].domain_middle_point = look_up_table_1_over_w[4].domain_middle_point_N;
    look_up_table_phtn_1_o_w_proj[4].evaluation_order = look_up_table_1_over_w[4].evaluation_order_N;
    look_up_table_phtn_1_o_w_proj[4].last_coeffs = look_up_table_1_over_w[4].evaluate_y(heval_r_0_2);

    look_up_table_phtn_1_o_w_proj[3].a = look_up_table_1_over_w[3].a;
    look_up_table_phtn_1_o_w_proj[3].b = look_up_table_1_over_w[3].b;
    look_up_table_phtn_1_o_w_proj[3].domain_half_range = look_up_table_1_over_w[3].domain_half_range_N;
    look_up_table_phtn_1_o_w_proj[3].domain_middle_point = look_up_table_1_over_w[3].domain_middle_point_N;
    look_up_table_phtn_1_o_w_proj[3].evaluation_order = look_up_table_1_over_w[3].evaluation_order_N;
    look_up_table_phtn_1_o_w_proj[3].last_coeffs = look_up_table_1_over_w[3].evaluate_y(heval_r_2_20);

    look_up_table_phtn_1_o_w_proj[2].a = look_up_table_1_over_w[2].a;
    look_up_table_phtn_1_o_w_proj[2].b = look_up_table_1_over_w[2].b;
    look_up_table_phtn_1_o_w_proj[2].domain_half_range = look_up_table_1_over_w[2].domain_half_range_N;
    look_up_table_phtn_1_o_w_proj[2].domain_middle_point = look_up_table_1_over_w[2].domain_middle_point_N;
    look_up_table_phtn_1_o_w_proj[2].evaluation_order = look_up_table_1_over_w[2].evaluation_order_N;
    look_up_table_phtn_1_o_w_proj[2].last_coeffs = look_up_table_1_over_w[2].evaluate_y(heval_r_20_80);

    look_up_table_phtn_1_o_w_proj[1].a = look_up_table_1_over_w[1].a;
    look_up_table_phtn_1_o_w_proj[1].b = look_up_table_1_over_w[1].b;
    look_up_table_phtn_1_o_w_proj[1].domain_half_range = look_up_table_1_over_w[1].domain_half_range_N;
    look_up_table_phtn_1_o_w_proj[1].domain_middle_point = look_up_table_1_over_w[1].domain_middle_point_N;
    look_up_table_phtn_1_o_w_proj[1].evaluation_order = look_up_table_1_over_w[1].evaluation_order_N;
    look_up_table_phtn_1_o_w_proj[1].last_coeffs = look_up_table_1_over_w[1].evaluate_y(heval_r_80_600);

    look_up_table_phtn_1_o_w_proj[0].a = look_up_table_1_over_w[0].a;
    look_up_table_phtn_1_o_w_proj[0].b = look_up_table_1_over_w[0].b;
    look_up_table_phtn_1_o_w_proj[0].domain_half_range = look_up_table_1_over_w[0].domain_half_range_N;
    look_up_table_phtn_1_o_w_proj[0].domain_middle_point = look_up_table_1_over_w[0].domain_middle_point_N;
    look_up_table_phtn_1_o_w_proj[0].evaluation_order = look_up_table_1_over_w[0].evaluation_order_N;
    look_up_table_phtn_1_o_w_proj[0].last_coeffs = look_up_table_1_over_w[0].evaluate_y(heval_r_600_2000);


    //differential cross section
    emission_name = path_to_coeffs + "dP_over_dt_dw/phtn_diff_crss_sctn_0-2.txt";
    // emission_name = path_to_coeffs + "prtl_rate_complete/phtn_partial_rate_0-2.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_diff_crss_sctn[4] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "dP_over_dt_dw/phtn_diff_crss_sctn_2-20.txt";
    // emission_name = path_to_coeffs + "prtl_rate_complete/phtn_partial_rate_2-20.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_diff_crss_sctn[3] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "dP_over_dt_dw/phtn_diff_crss_sctn_20-80.txt";
    // emission_name = path_to_coeffs + "prtl_rate_complete/phtn_partial_rate_20-80.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_diff_crss_sctn[2] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "dP_over_dt_dw/phtn_diff_crss_sctn_80-600.txt";
    // emission_name = path_to_coeffs + "prtl_rate_complete/phtn_partial_rate_80-600.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_diff_crss_sctn[1] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "dP_over_dt_dw/phtn_diff_crss_sctn_600-2000.txt";
    // emission_name = path_to_coeffs + "prtl_rate_complete/phtn_partial_rate_600-2000.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_diff_crss_sctn[0] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();
}

//the function below returns the rate of photon emission for an electron
//having  a certain chi and gamma (gamma is the relativistic factor of 
//the electron, but it can also be interpreted as its energy normalized in 
//units of m_e c^2)
// double SFQED_Processes::SFQED_PHTN_emission_rate(const double &gamma, const double &chi) const {
        
//     //coefficient for the rate of emission
//     //multiply by \tilde{W_rad} to get the rate of emission
//     double coefW_rad = coef_rate_Wrad * chi / gamma;
                
//     //low impact Branch to calculate \tilde{W_rad} and then the rate
//     //---------------------------------------------------
//     /**/
//     //lookup tables boolean selectors
//     // bool chi_0_2 = chi <= bound_chi_1st,
//     //     chi_2_20 = chi <= bound_chi_2nd,
//     //     chi_20_80 = chi <= bound_chi_3rd,
//     //     chi_80_600 = chi <= bound_chi_4th;
//     //     //chi_600_2000 = chi <= 2000.;

//     // int lookup_index = chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600;
//     // lookup_index = (1 - std::signbit(lookup_index-1))*(32 - __builtin_clz(lookup_index));

//     int lookup_index;

//     if(chi <= bound_chi_1st){
//         lookup_index = 4;
//     }else if(chi <= bound_chi_2nd){
//         lookup_index = 3;
//     }else if(chi <= bound_chi_3rd){
//         lookup_index = 2;
//     }else if(chi <= bound_chi_4th){
//         lookup_index = 1;
//     }else{
//         lookup_index = 0;
//     }

//     //IMPORTANT: through this constructs we means
//     //to default every chi > 600 case into the
//     //600 < chi <= 2000 one
//     const Chebyshev_User_1D chebyshev_phtn_rate = (look_up_table_phtn_mx_rt[lookup_index]);
                
//     return coefW_rad * chebyshev_phtn_rate.evaluate(chi);
// }

double SFQED_Processes::SFQED_LCFA_phtn_nrg_aux(const double &chi, const double &rnd, const int& lookup_index) const{

    double v;

    const double low_r_limit = (lu_table_lower_r_bounds[lookup_index]),
            inverse_r_limit = (lu_table_inverse_r_bounds[lookup_index]),
            high_r_limit = (lu_table_upper_r_bounds[lookup_index]);

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
        const Chebyshev_User_1D chebyshev_phtn_rate =
                (look_up_table_phtn_mx_rt[lookup_index]);

        //std::cout << "3\n" << std::flush;

        //photon emission rate (just the integral without the coefficient \chi / gamma)
        //double tildeWrad = chebyshev_phtn_rate.evaluate(chi);
        double rnd_tildeWrad = rnd * chebyshev_phtn_rate.evaluate(chi);
        //////////////////////////////////////

        //std::cout << "4\n" << std::flush;

        //this is w
        v = soft_ph_coef * rnd_tildeWrad;
    }
    //otherwise we use the inverse equation Chebyshev approximation
    else if (rnd < inverse_r_limit){

        //std::cout << "11\n" << std::flush;

        //this evaluation, inverse of g(\chi, w) - r*W_{rad}(\chi) = 0,
        //where \chi and r are fixed at this stage, will give a w value
        const Chebyshev_User_2D chebyshev_phtn_momw =
                (look_up_table_phtn_momw[lookup_index]);

        //std::cout << "12\n" << std::flush;

        //this is w
        v = chebyshev_phtn_momw.evaluate(chi, rnd);

        //std::cout << "13\n" << std::flush;
    }
    //in case the random number corresponds to some high energy value we resort to the inverse Chebyshev approx
    else if(rnd < high_r_limit){
        const Chebyshev_User_2D chebyshev_phtn_1_o_w =
                (look_up_table_1_over_w[lookup_index]);

        //this is w
        v = 1. / chebyshev_phtn_1_o_w.evaluate(chi, rnd);
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
        const Chebyshev_User_1D chebyshev_phtn_rate_high =
                (look_up_table_phtn_mx_rt[lookup_index]);

        //std::cout << "6\n" << std::flush;

        //photon emission rate (just the integral without the coefficient \chi / gamma)
        double tildeWrad = chebyshev_phtn_rate_high.evaluate(chi);
        double rnd_tildeWrad = rnd * tildeWrad;
        //////////////////////////////////////

        //std::cout << "7\n" << std::flush;

        //get the partial photon emission rate, by integrating
        //the differential probability from 0 to a certain w
        //(keep chi fixed). Please, notice we are evaluating
        //this on the upper bound of validity of the approximation
        const Chebyshev_User_2D chebyshev_phtn_prtl_rate =
                (look_up_table_phtn_prtl_rt[lookup_index]);

        //std::cout << "8\n" << std::flush;

        //fixed pivots******************************
        // double w0 = *(lu_high_tail_w_bounds[lookup_index]),
        //         v0 = *(lu_high_tail_v_bounds[lookup_index]);
        //t*****************************************
        
        //newly computed pivots/////////////////////
        //this new attempt requires us to compute w through the inverse (maybe it is too much time consuming?)
        const Chebyshev_User_1D chebyshev_phtn_1_o_w_proj =
                (look_up_table_phtn_1_o_w_proj[lookup_index]);

        double w0 = 1. / chebyshev_phtn_1_o_w_proj.evaluate(chi);
        double v0 = w0 * w0 * w0;

        // Chebyshev_User_2D chebyshev_phtn_momw =
        //         *(look_up_table_phtn_momw[lookup_index]);
        // double w0 =  chebyshev_phtn_momw.evaluate(chi, high_r_limit);
        // double v0 = w0 * w0 * w0;
        ////////////////////////////////////////////

        // double g_of_w_chi = chebyshev_phtn_prtl_rate->evaluate(chi, w0);
        double g_of_w_chi = chebyshev_phtn_prtl_rate.evaluate(chi, w0);

        //std::cout << "9\n" << std::flush;

        //be careful to when v is null
        v = std::abs((tildeWrad - rnd_tildeWrad) / (tildeWrad - g_of_w_chi));
                
        //we use the cubic root to retrieve the w value
        v = cbrt(v0 - log(v));

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
    // bool chi_0_2 = chi <= bound_chi_1st,
    //     chi_2_20 = chi <= bound_chi_2nd,
    //     chi_20_80 = chi <= bound_chi_3rd,
    //     chi_80_600 = chi <= bound_chi_4th;
    //     //chi_600_2000 = chi <= 2000.;

    // int lookup_index = chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600;
    // lookup_index = (1 - std::signbit(lookup_index-1))*(32 - __builtin_clz(lookup_index));

    int lookup_index;

    if(chi <= bound_chi_1st){
        lookup_index = 4;
    }else if(chi <= bound_chi_2nd){
        lookup_index = 3;
    }else if(chi <= bound_chi_3rd){
        lookup_index = 2;
    }else if(chi <= bound_chi_4th){
        lookup_index = 1;
    }else{
        lookup_index = 0;
    }

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
    // bool chi_0_2 = chi <= bound_chi_1st,
    //     chi_2_20 = chi <= bound_chi_2nd,
    //     chi_20_80 = chi <= bound_chi_3rd,
    //     chi_80_600 = chi <= bound_chi_4th;
    //     //chi_600_2000 = chi <= 2000.;

    // int lookup_index = chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600;
    // lookup_index = (1 - std::signbit(lookup_index-1))*(32 - __builtin_clz(lookup_index));

    int lookup_index;

    if(chi <= bound_chi_1st){
        lookup_index = 4;
    }else if(chi <= bound_chi_2nd){
        lookup_index = 3;
    }else if(chi <= bound_chi_3rd){
        lookup_index = 2;
    }else if(chi <= bound_chi_4th){
        lookup_index = 1;
    }else{
        lookup_index = 0;
    }
    
    //transform the threshold energy into a w value
    double LCFA_limit_w  = cbrt(LCFA_limit / (1.5 * chi * (gamma - LCFA_limit)));
    

    //in case the phtn nrg threshold is very close to the particle nrg
    //or it is above the upper bound of the approximated domain
    //we don't want the emission to occur 
    if(LCFA_limit > 0.75 * gamma || LCFA_limit_w > (lu_table_upper_w_bounds[lookup_index])){
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

        // const Chebyshev_User_2D * const phtn_diff_crss_sctn =
        //         (look_up_table_phtn_diff_crss_sctn[lookup_index]);
        const Chebyshev_User_2D phtn_diff_crss_sctn =
                (look_up_table_phtn_diff_crss_sctn[lookup_index]);

        //compute the threshold value
        rs_dw_over_depsgamma = (1. + 1.5 * chi * LCFA_limit_w * LCFA_limit_w * LCFA_limit_w) / LCFA_limit_w;
        diff_prob_LCFA = phtn_diff_crss_sctn.evaluate(chi, LCFA_limit_w) * rs_dw_over_depsgamma * rs_dw_over_depsgamma;

        //compute the differential probability associated to the LCFA supposed phtn nrg
        //(pay attention not to overwrite the w value energy)
        rs_dw_over_depsgamma = (1. + 1.5 * chi * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w) / emitted_phtn_LCFA_nrg_w;
        diff_probability = phtn_diff_crss_sctn.evaluate(chi, emitted_phtn_LCFA_nrg_w) * rs_dw_over_depsgamma * rs_dw_over_depsgamma;

        //now I use the rejection sampling method (employing rnd2) to determine if I have
        //to keep or reject the just occurred emission
        if(diff_probability * rnd2 > diff_prob_LCFA){
            return 0.;
        }

    }

    //convert to a proper photon energy
    emitted_phtn_LCFA_nrg_w = 3.0 * chi * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w;
    emitted_phtn_LCFA_nrg_w = gamma * emitted_phtn_LCFA_nrg_w / (2.0 + emitted_phtn_LCFA_nrg_w);

    return emitted_phtn_LCFA_nrg_w;
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
    delete[] look_up_table_phtn_mx_rt[4].last_coeffs;
    delete[] look_up_table_phtn_mx_rt[3].last_coeffs;
    delete[] look_up_table_phtn_mx_rt[2].last_coeffs;
    delete[] look_up_table_phtn_mx_rt[1].last_coeffs;
    delete[] look_up_table_phtn_mx_rt[0].last_coeffs;

    //purge partial (from 0 to w) emission rates
    delete[] look_up_table_phtn_prtl_rt[4].last_coeffs;
    delete[] look_up_table_phtn_prtl_rt[3].last_coeffs;
    delete[] look_up_table_phtn_prtl_rt[2].last_coeffs;
    delete[] look_up_table_phtn_prtl_rt[1].last_coeffs;
    delete[] look_up_table_phtn_prtl_rt[0].last_coeffs;

    //purge photon emission energy
    delete[] look_up_table_phtn_momw[4].last_coeffs;
    delete[] look_up_table_phtn_momw[3].last_coeffs;
    delete[] look_up_table_phtn_momw[2].last_coeffs;
    delete[] look_up_table_phtn_momw[1].last_coeffs;
    delete[] look_up_table_phtn_momw[0].last_coeffs;

    //purge 1/w
    delete[] look_up_table_1_over_w[4].last_coeffs;
    delete[] look_up_table_1_over_w[3].last_coeffs;
    delete[] look_up_table_1_over_w[2].last_coeffs;
    delete[] look_up_table_1_over_w[1].last_coeffs;
    delete[] look_up_table_1_over_w[0].last_coeffs;

    //purge 1D projection of 1/w
    delete[] look_up_table_phtn_1_o_w_proj[4].last_coeffs;
    delete[] look_up_table_phtn_1_o_w_proj[3].last_coeffs;
    delete[] look_up_table_phtn_1_o_w_proj[2].last_coeffs;
    delete[] look_up_table_phtn_1_o_w_proj[1].last_coeffs;
    delete[] look_up_table_phtn_1_o_w_proj[0].last_coeffs;

    //purge diff crss sctn
    delete[] look_up_table_phtn_diff_crss_sctn[4].last_coeffs;
    delete[] look_up_table_phtn_diff_crss_sctn[3].last_coeffs;
    delete[] look_up_table_phtn_diff_crss_sctn[2].last_coeffs;
    delete[] look_up_table_phtn_diff_crss_sctn[1].last_coeffs;
    delete[] look_up_table_phtn_diff_crss_sctn[0].last_coeffs;
}

/*****************/
/*PAIR PRODUCTION*/
/*****************/

// double asympt_pair_emission_rate_low_k(double k){
//     return (27. * pi * k)/(16. * sqrt(2.)) * exp(-(8/(3 * k))) * (1. - 11./64.*k + 7585./73728.*k*k);
// }

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

//
    look_up_table_pair_prd_rt[5].domain_half_range = 1;
    look_up_table_pair_prd_rt[5].domain_middle_point = 1;
    look_up_table_pair_prd_rt[5].evaluation_order = 0;
//

    std::string emission_name = path_to_coeffs + "pair_rate/pair_production_rate_0-2.txt";
    std::ifstream in_file(emission_name.c_str());
    look_up_table_pair_prd_rt[4] = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_rate/pair_production_rate_2-20.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_prd_rt[3] = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_rate/pair_production_rate_20-80.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_prd_rt[2] = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_rate/pair_production_rate_80-600.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_prd_rt[1] = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_rate/pair_production_rate_600-2000.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_prd_rt[0] = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();


//
    look_up_table_pair_prtl_rt[5].domain_half_range_M = 1;
    look_up_table_pair_prtl_rt[5].domain_half_range_N = 1;
    look_up_table_pair_prtl_rt[5].domain_middle_point_M = 1;
    look_up_table_pair_prtl_rt[5].domain_middle_point_N = 1;
    look_up_table_pair_prtl_rt[5].evaluation_order_M = 0;
    look_up_table_pair_prtl_rt[5].evaluation_order_N = 0;
//

    emission_name = path_to_coeffs + "pair_prtl_rate/pair_partial_rate_0-2.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_prtl_rt[4] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_prtl_rate/pair_partial_rate_2-20.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_prtl_rt[3] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_prtl_rate/pair_partial_rate_20-80.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_prtl_rt[2] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_prtl_rate/pair_partial_rate_80-600.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_prtl_rt[1] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_prtl_rate/pair_partial_rate_600-2000.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_prtl_rt[0] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();


    //photon energy in v
//
    emission_name = path_to_coeffs + "v/pair_nrgs_v_001-03.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs[5] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();
//

    emission_name = path_to_coeffs + "v/pair_nrgs_v_0-2.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs[4] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "v/pair_nrgs_v_2-20.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs[3] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "v/pair_nrgs_v_20-80.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs[2] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "v/pair_nrgs_v_80-600.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs[1] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "v/pair_nrgs_v_600-2000.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs[0] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();


    //photon energy in v high part
//
    emission_name = path_to_coeffs + "v/pair_nrgs_v_high_001-03.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs_high[5] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();
//

    emission_name = path_to_coeffs + "v/pair_nrgs_v_high_0-2.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs_high[4] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "v/pair_nrgs_v_high_2-20.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs_high[3] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "v/pair_nrgs_v_high_20-80.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs_high[2] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "v/pair_nrgs_v_high_80-600.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs_high[1] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "v/pair_nrgs_v_high_600-2000.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs_high[0] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();


    //inverse 1D projection (to be fixed)
//
    look_up_table_pair_v_nrgs_proj[5].a = look_up_table_pair_v_nrgs_high[5].a;
    look_up_table_pair_v_nrgs_proj[5].b = look_up_table_pair_v_nrgs_high[5].b;
    look_up_table_pair_v_nrgs_proj[5].domain_half_range = look_up_table_pair_v_nrgs_high[5].domain_half_range_N;
    look_up_table_pair_v_nrgs_proj[5].domain_middle_point = look_up_table_pair_v_nrgs_high[5].domain_middle_point_N;
    look_up_table_pair_v_nrgs_proj[5].evaluation_order = look_up_table_pair_v_nrgs_high[5].evaluation_order_N;
    look_up_table_pair_v_nrgs_proj[5].last_coeffs = look_up_table_pair_v_nrgs_high[5].evaluate_y(pair_r_exp_limit_001_03);
//

    look_up_table_pair_v_nrgs_proj[4].a = look_up_table_pair_v_nrgs_high[4].a;
    look_up_table_pair_v_nrgs_proj[4].b = look_up_table_pair_v_nrgs_high[4].b;
    look_up_table_pair_v_nrgs_proj[4].domain_half_range = look_up_table_pair_v_nrgs_high[4].domain_half_range_N;
    look_up_table_pair_v_nrgs_proj[4].domain_middle_point = look_up_table_pair_v_nrgs_high[4].domain_middle_point_N;
    look_up_table_pair_v_nrgs_proj[4].evaluation_order = look_up_table_pair_v_nrgs_high[4].evaluation_order_N;
    look_up_table_pair_v_nrgs_proj[4].last_coeffs = look_up_table_pair_v_nrgs_high[4].evaluate_y(pair_r_exp_limit_0_2);

    look_up_table_pair_v_nrgs_proj[3].a = look_up_table_pair_v_nrgs_high[3].a;
    look_up_table_pair_v_nrgs_proj[3].b = look_up_table_pair_v_nrgs_high[3].b;
    look_up_table_pair_v_nrgs_proj[3].domain_half_range = look_up_table_pair_v_nrgs_high[3].domain_half_range_N;
    look_up_table_pair_v_nrgs_proj[3].domain_middle_point = look_up_table_pair_v_nrgs_high[3].domain_middle_point_N;
    look_up_table_pair_v_nrgs_proj[3].evaluation_order = look_up_table_pair_v_nrgs_high[3].evaluation_order_N;
    look_up_table_pair_v_nrgs_proj[3].last_coeffs = look_up_table_pair_v_nrgs_high[3].evaluate_y(pair_r_exp_limit_2_20);

    look_up_table_pair_v_nrgs_proj[2].a = look_up_table_pair_v_nrgs_high[2].a;
    look_up_table_pair_v_nrgs_proj[2].b = look_up_table_pair_v_nrgs_high[2].b;
    look_up_table_pair_v_nrgs_proj[2].domain_half_range = look_up_table_pair_v_nrgs_high[2].domain_half_range_N;
    look_up_table_pair_v_nrgs_proj[2].domain_middle_point = look_up_table_pair_v_nrgs_high[2].domain_middle_point_N;
    look_up_table_pair_v_nrgs_proj[2].evaluation_order = look_up_table_pair_v_nrgs_high[2].evaluation_order_N;
    look_up_table_pair_v_nrgs_proj[2].last_coeffs = look_up_table_pair_v_nrgs_high[2].evaluate_y(pair_r_exp_limit_20_80);

    look_up_table_pair_v_nrgs_proj[1].a = look_up_table_pair_v_nrgs_high[1].a;
    look_up_table_pair_v_nrgs_proj[1].b = look_up_table_pair_v_nrgs_high[1].b;
    look_up_table_pair_v_nrgs_proj[1].domain_half_range = look_up_table_pair_v_nrgs_high[1].domain_half_range_N;
    look_up_table_pair_v_nrgs_proj[1].domain_middle_point = look_up_table_pair_v_nrgs_high[1].domain_middle_point_N;
    look_up_table_pair_v_nrgs_proj[1].evaluation_order = look_up_table_pair_v_nrgs_high[1].evaluation_order_N;
    look_up_table_pair_v_nrgs_proj[1].last_coeffs = look_up_table_pair_v_nrgs_high[1].evaluate_y(pair_r_exp_limit_80_600);

    look_up_table_pair_v_nrgs_proj[0].a = look_up_table_pair_v_nrgs_high[0].a;
    look_up_table_pair_v_nrgs_proj[0].b = look_up_table_pair_v_nrgs_high[0].b;
    look_up_table_pair_v_nrgs_proj[0].domain_half_range = look_up_table_pair_v_nrgs_high[0].domain_half_range_N;
    look_up_table_pair_v_nrgs_proj[0].domain_middle_point = look_up_table_pair_v_nrgs_high[0].domain_middle_point_N;
    look_up_table_pair_v_nrgs_proj[0].evaluation_order = look_up_table_pair_v_nrgs_high[0].evaluation_order_N;
    look_up_table_pair_v_nrgs_proj[0].last_coeffs = look_up_table_pair_v_nrgs_high[0].evaluate_y(pair_r_exp_limit_600_2000);
}

// double  SFQED_Processes::SFQED_PAIR_creation_rate(const double &gamma, const double &chi) const {
    
//     //coefficient for the rate of emission
//     //multiply by \tilde{W_rad} to get the rate of emission
//     double coefW_pair = coef_rate_Wpair / gamma;
                
//     //low impact Branch to calculate \tilde{W_rad} and then the rate
//     //---------------------------------------------------
//     /**/
//     //lookup tables boolean selectors
//     // bool chi_001_03 = chi <= bound_kappa_0th, 
//     //     chi_0_2 = chi <= bound_kappa_1st,
//     //     chi_2_20 = chi <= bound_kappa_2nd,
//     //     chi_20_80 = chi <= bound_kappa_3rd,
//     //     chi_80_600 = chi <= bound_kappa_4th;
//     //     //chi_600_2000 = chi <= 2000.;
                
//     // int lookup_index = chi_001_03*16 + chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600;
//     // lookup_index = (1 - std::signbit(lookup_index-1))*(32 - __builtin_clz(lookup_index));
                

//     int lookup_index;

//     if(chi <= bound_kappa_0th){
//         // lookup_index = 5;
//         return coefW_pair * //asympt_pair_emission_rate_low_k(chi);
//                 (27. * pigreek * chi)/(16. * sqrt(2.)) * exp(-(8/(3 * chi))) * (1. - 11./64.*chi + 7585./73728.*chi*chi);
//     }else if(chi <= bound_kappa_1st){
//         lookup_index = 4;
//     }else if(chi <= bound_kappa_2nd){
//         lookup_index = 3;
//     }else if(chi <= bound_kappa_3rd){
//         lookup_index = 2;
//     }else if(chi <= bound_kappa_4th){
//         lookup_index = 1;
//     }else{
//         lookup_index = 0;
//     }

//     //IMPORTANT: through this constructs we mean
//     //to default every chi > 600 case into the
//     //600 < chi <= 2000 one
//     const Chebyshev_User_1D chebyshev_pair_rate = 
//                 (look_up_table_pair_prd_rt[lookup_index]);                        
                
//     return coefW_pair * chebyshev_pair_rate.evaluate(chi);//
// }

// double  SFQED_Processes::SFQED_PAIR_creation_rate_fast(const double &gamma, const double &chi) const {
    
//     //coefficient for the rate of emission
//     //multiply by \tilde{W_rad} to get the rate of emission
//     double coefW_pair = coef_rate_Wpair / gamma;
                
//     //low impact Branch to calculate \tilde{W_rad} and then the rate
//     //---------------------------------------------------
//     /**/
//     //lookup tables boolean selectors
//     // bool chi_001_03 = chi <= bound_kappa_0th, 
//     //     chi_0_2 = chi <= bound_kappa_1st,
//     //     chi_2_20 = chi <= bound_kappa_2nd,
//     //     chi_20_80 = chi <= bound_kappa_3rd,
//     //     chi_80_600 = chi <= bound_kappa_4th;
//     //     //chi_600_2000 = chi <= 2000.;
                
//     // int lookup_index = chi_001_03*16 + chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600;
//     // lookup_index = (1 - std::signbit(lookup_index-1))*(32 - __builtin_clz(lookup_index));
                

//     int lookup_index;

//     if(chi <= bound_kappa_0th){
//         lookup_index = 5;
//     }else if(chi <= bound_kappa_1st){
//         lookup_index = 4;
//     }else if(chi <= bound_kappa_2nd){
//         lookup_index = 3;
//     }else if(chi <= bound_kappa_3rd){
//         lookup_index = 2;
//     }else if(chi <= bound_kappa_4th){
//         lookup_index = 1;
//     }else{
//         lookup_index = 0;
//     }

//     //IMPORTANT: through this constructs we mean
//     //to default every chi > 600 case into the
//     //600 < chi <= 2000 one
//     const Chebyshev_User_1D chebyshev_pair_rate = 
//                 (look_up_table_pair_prd_rt[lookup_index]);                        
                
//     return coefW_pair * chebyshev_pair_rate.evaluate(chi);//
// }

double SFQED_Processes::SFQED_PAIR_emitted_electron_energy_aux(const int &lookup_index, const double &chi, const double &rescaled_rnd) const{

    const double low_r_limit = (lu_table_pair_lower_r_bounds[lookup_index]),
            high_r_limit = (lu_table_pair_upper_r_bounds[lookup_index]);

    // std::cout << "1\n";

    double v;

    if(rescaled_rnd < low_r_limit){

        // std::cout << "3\n";

        const Chebyshev_User_2D chebyshev_pair_v_nrgs_inv =
                (look_up_table_pair_v_nrgs[lookup_index]);

        // std::cout << "4\n";

        v = chebyshev_pair_v_nrgs_inv.evaluate(chi, rescaled_rnd);

        // std::cout << "6\n";
    }
    else if(rescaled_rnd < high_r_limit){

        // std::cout << "a\n";

        const Chebyshev_User_2D chebyshev_pair_v_nrgs_inv_high =
                (look_up_table_pair_v_nrgs_high[lookup_index]);

        v = chebyshev_pair_v_nrgs_inv_high.evaluate(chi, rescaled_rnd);

        // std::cout << "b\n";
    }
    //in the off chance case that we happen to be beyond the
    //high r threshold, we will assume to be dealing directly
    //with v = 0. Since v = (\epsilon_{\gamma} - 2 \epsilon_e) / \epsilon_{\gamma}
    else{

        // std::cout << "7\n";

        //retrieve the appropriate pair rate set of coefficients
        const Chebyshev_User_1D chebyshev_pair_rate =
                (look_up_table_pair_prd_rt[lookup_index]);

        // std::cout << "8\n";

        //compute the rate of pair emission
        double tildeWrad = chebyshev_pair_rate.evaluate(chi) / 3.;
                            
        double rnd_tildeWrad = rescaled_rnd * tildeWrad;

        // std::cout << "9\n";
        
        //compute the pivot at which the exponential approximation will be applied
        const Chebyshev_User_1D chebyshev_pair_v_nrgs_proj =
                (look_up_table_pair_v_nrgs_proj[lookup_index]);

        // std::cout << "10\n";

         double v0 = chebyshev_pair_v_nrgs_proj.evaluate(chi);
        v = 8.0 / (3.0 * chi * (1.0 - v0*v0));

        // std::cout << "11\n";

        //integrate the pair creation differential probability
        //from the just computed pivot up to 1
        const Chebyshev_User_2D chebyshev_pair_prtl_rate =
                (look_up_table_pair_prtl_rt[lookup_index]);

        // std::cout << "12\n";

        double tail_integral = chebyshev_pair_prtl_rate.evaluate(chi, v0);

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
    // bool chi_001_03 = chi <= bound_kappa_0th, 
    //     chi_0_2 = chi <= bound_kappa_1st,
    //     chi_2_20 = chi <= bound_kappa_2nd,
    //     chi_20_80 = chi <= bound_kappa_3rd,
    //     chi_80_600 = chi <= bound_kappa_4th;
    //     //chi_600_2000 = chi <= 2000.;

    // int lookup_index = chi_001_03*16 + chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600;
    // lookup_index = (1 - std::signbit(lookup_index-1))*(32 - __builtin_clz(lookup_index));


    int lookup_index;

    if(chi <= bound_kappa_0th){
        lookup_index = 5;
    }else if(chi <= bound_kappa_1st){
        lookup_index = 4;
    }else if(chi <= bound_kappa_2nd){
        lookup_index = 3;
    }else if(chi <= bound_kappa_3rd){
        lookup_index = 2;
    }else if(chi <= bound_kappa_4th){
        lookup_index = 1;
    }else{
        lookup_index = 0;
    }

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
    // bool chi_0_2 = chi <= bound_kappa_1st,
    //     chi_2_20 = chi <= bound_kappa_2nd,
    //     chi_20_80 = chi <= bound_kappa_3rd,
    //     chi_80_600 = chi <= bound_kappa_4th;
    //     //chi_600_2000 = chi <= 2000.;

    // int lookup_index = chi_0_2*8 + chi_2_20*4 + chi_20_80*2 + chi_80_600;
    // lookup_index = (1 - std::signbit(lookup_index-1))*(32 - __builtin_clz(lookup_index));

    int lookup_index;

    if(chi <= bound_kappa_1st){
        lookup_index = 4;
    }else if(chi <= bound_kappa_2nd){
        lookup_index = 3;
    }else if(chi <= bound_kappa_3rd){
        lookup_index = 2;
    }else if(chi <= bound_kappa_4th){
        lookup_index = 1;
    }else{
        lookup_index = 0;
    }

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
    //purge all coefficients
    delete[] look_up_table_pair_prd_rt[4].last_coeffs;
    delete[] look_up_table_pair_prd_rt[3].last_coeffs;
    delete[] look_up_table_pair_prd_rt[2].last_coeffs;
    delete[] look_up_table_pair_prd_rt[1].last_coeffs;
    delete[] look_up_table_pair_prd_rt[0].last_coeffs;

    delete[] look_up_table_pair_prtl_rt[4].last_coeffs;
    delete[] look_up_table_pair_prtl_rt[3].last_coeffs;
    delete[] look_up_table_pair_prtl_rt[2].last_coeffs;
    delete[] look_up_table_pair_prtl_rt[1].last_coeffs;
    delete[] look_up_table_pair_prtl_rt[0].last_coeffs;

    delete[] look_up_table_pair_v_nrgs[5].last_coeffs;
    delete[] look_up_table_pair_v_nrgs[4].last_coeffs;
    delete[] look_up_table_pair_v_nrgs[3].last_coeffs;
    delete[] look_up_table_pair_v_nrgs[2].last_coeffs;
    delete[] look_up_table_pair_v_nrgs[1].last_coeffs;
    delete[] look_up_table_pair_v_nrgs[0].last_coeffs;

    delete[] look_up_table_pair_v_nrgs_high[5].last_coeffs;
    delete[] look_up_table_pair_v_nrgs_high[4].last_coeffs;
    delete[] look_up_table_pair_v_nrgs_high[3].last_coeffs;
    delete[] look_up_table_pair_v_nrgs_high[2].last_coeffs;
    delete[] look_up_table_pair_v_nrgs_high[1].last_coeffs;
    delete[] look_up_table_pair_v_nrgs_high[0].last_coeffs;

    delete[] look_up_table_pair_v_nrgs_proj[5].last_coeffs;
    delete[] look_up_table_pair_v_nrgs_proj[4].last_coeffs;
    delete[] look_up_table_pair_v_nrgs_proj[3].last_coeffs;
    delete[] look_up_table_pair_v_nrgs_proj[2].last_coeffs;
    delete[] look_up_table_pair_v_nrgs_proj[1].last_coeffs;
    delete[] look_up_table_pair_v_nrgs_proj[0].last_coeffs;
}