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

#include "Nonlinear_Inverse_Compton_Scattering.h"


#include <cmath>

/********/
/* LCFA */
/********/

void Nonlinear_Inverse_Compton_Scattering::SFQED_init_PHTN_emission(std::string path_to_coeffs){

    // init_phtn_mx_tables();

    //photon emission rate
    std::string emission_name = path_to_coeffs + "phtn_rate/phtn_emission_rate_0-2.txt";
    std::ifstream in_file(emission_name.c_str());
    look_up_table_phtn_mx_rt[4] = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "phtn_rate/phtn_emission_rate_2-20.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_mx_rt[3] = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "phtn_rate/phtn_emission_rate_20-80.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_mx_rt[2] = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "phtn_rate/phtn_emission_rate_80-600.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_mx_rt[1] = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "phtn_rate/phtn_emission_rate_600-2000.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_mx_rt[0] = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();


    //photon partial emission rate g(w,\chi) = \int_0^w{Wrad(\chi,w_{integrated})dw_{integrated}}
    emission_name = path_to_coeffs + "phtn_prtl_rate/phtn_partial_rate_0-2.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_prtl_rt[4] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "phtn_prtl_rate/phtn_partial_rate_2-20.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_prtl_rt[3] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "phtn_prtl_rate/phtn_partial_rate_20-80.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_prtl_rt[2] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "phtn_prtl_rate/phtn_partial_rate_80-600.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_prtl_rt[1] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "phtn_prtl_rate/phtn_partial_rate_600-2000.txt";
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
    // emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_0-2_plus.txt";
    // emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_0-2_extended.txt";
    in_file.open(emission_name.c_str());
    look_up_table_1_over_w[4] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_2-20.txt";
    // emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_2-20_plus.txt";
    // emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_2-20_extended.txt";
    in_file.open(emission_name.c_str());
    look_up_table_1_over_w[3] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_20-80.txt";
    // emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_20-80_plus.txt";
    // emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_20-80_extended.txt";
    in_file.open(emission_name.c_str());
    look_up_table_1_over_w[2] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_80-600.txt";
    // emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_80-600_plus.txt";
    // emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_80-600_extended.txt";
    in_file.open(emission_name.c_str());
    look_up_table_1_over_w[1] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_600-2000.txt";
    // emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_600-2000_plus.txt";
    // emission_name = path_to_coeffs + "1_over_w/phtn_1_over_w_nrg_600-2000_extended.txt";
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

}

//auxiliary inner function used inside the LCFA energy method
//it computes the v-value energy of the emitted photon
double Nonlinear_Inverse_Compton_Scattering::SFQED_LCFA_phtn_nrg_aux(const double &chi, const double &rnd, const int& lookup_index) const {

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

        //std::cout << "3\n" << std::flush;

        //photon emission rate (just the integral without the coefficient \chi / gamma)
        //double tildeWrad = chebyshev_phtn_rate.evaluate(chi);
        double rnd_tildeWrad = rnd * (look_up_table_phtn_mx_rt[lookup_index]).evaluate(chi);
        //////////////////////////////////////

        //std::cout << "4\n" << std::flush;

        //this is w
        v = soft_ph_coef * rnd_tildeWrad;

        //this is v = w^3
        v *= v * v;
    }
    //otherwise we use the inverse equation Chebyshev approximation
    else if (rnd < inverse_r_limit){

        //std::cout << "11\n" << std::flush;

        //std::cout << "12\n" << std::flush;

        //this is w
        //the evaluation of these coefficients results in the
        //inverse of g(\chi, w) - r*W_{rad}(\chi) = 0,
        //where \chi and r are fixed at this stage, will give a w value
        v = (look_up_table_phtn_momw[lookup_index]).evaluate(chi, rnd);

        //this is v = w^3
        v *= v * v;

        //std::cout << "13\n" << std::flush;
    }
    //in case the random number corresponds to some high energy value we resort to the inverse Chebyshev approx
    else if(rnd < high_r_limit){
       
        //this is w
        v = 1. / (look_up_table_1_over_w[lookup_index]).evaluate(chi, rnd);

        //this is v = w^3
        v *= v * v;
    }
    //in the off chance that we happen to be beyond the
    //high energy tail we will employ the exponential 
    //approximation (see notes)
    else{

        //std::cout << "5\n" << std::flush;

        //std::cout << "6\n" << std::flush;

        //photon emission rate (just the integral without the coefficient \chi / gamma)
        double tildeWrad = (look_up_table_phtn_mx_rt[lookup_index]).evaluate(chi);
        double rnd_tildeWrad = rnd * tildeWrad;
        //////////////////////////////////////

        //std::cout << "7\n" << std::flush;

        //std::cout << "8\n" << std::flush;
        
        //pivot w0 from which the exponential tail is computed
        double w0 = 1. / (look_up_table_phtn_1_o_w_proj[lookup_index]).evaluate(chi);
        //this variable is cumbersome but this
        // high tail is statistically never called
        double v0 = w0 * w0 * w0;

        //this evaluation brings us to the partial photon emission rate, by integrating
        //the differential probability from 0 to a certain w
        //(keep chi fixed). Please, notice we are evaluating
        //this on the upper bound of validity of the approximation w0
        double g_of_w_chi = (look_up_table_phtn_prtl_rt[lookup_index]).evaluate(chi, w0);

        //std::cout << "9\n" << std::flush;

        //be careful to when v is null
        v = std::abs((tildeWrad - rnd_tildeWrad) / (tildeWrad - g_of_w_chi));
                
        //this is v
        v = v0 - log(v);

    }

    //return the photon energy disguised as w value
    return v;
}

//auxiliary inner function used inside the BLCFA energy method
//it computes the w-value energy of the emitted photon
double Nonlinear_Inverse_Compton_Scattering::SFQED_BLCFA_phtn_nrg_aux(const double &chi, const double &rnd, const int& lookup_index) const {

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

        //std::cout << "3\n" << std::flush;

        //photon emission rate (just the integral without the coefficient \chi / gamma)
        //double tildeWrad = chebyshev_phtn_rate.evaluate(chi);
        double rnd_tildeWrad = rnd * (look_up_table_phtn_mx_rt[lookup_index]).evaluate(chi);
        //////////////////////////////////////

        //std::cout << "4\n" << std::flush;

        //this is w
        v = soft_ph_coef * rnd_tildeWrad;
    }
    //otherwise we use the inverse equation Chebyshev approximation
    else if (rnd < inverse_r_limit){

        //std::cout << "11\n" << std::flush;

        //std::cout << "12\n" << std::flush;

        //this is w
        //the evaluation of these coefficients results in the
        //inverse of g(\chi, w) - r*W_{rad}(\chi) = 0,
        //where \chi and r are fixed at this stage, will give a w value
        v = (look_up_table_phtn_momw[lookup_index]).evaluate(chi, rnd);

        //std::cout << "13\n" << std::flush;
    }
    //in case the random number corresponds to some high energy value we resort to the inverse Chebyshev approx
    else if(rnd < high_r_limit){
       
        //this is w
        v = 1. / (look_up_table_1_over_w[lookup_index]).evaluate(chi, rnd);

    }
    //in the off chance that we happen to be beyond the
    //high energy tail we will employ the exponential 
    //approximation (see notes)
    else{

        //std::cout << "5\n" << std::flush;

        //std::cout << "6\n" << std::flush;

        //photon emission rate (just the integral without the coefficient \chi / gamma)
        double tildeWrad = (look_up_table_phtn_mx_rt[lookup_index]).evaluate(chi);
        double rnd_tildeWrad = rnd * tildeWrad;
        //////////////////////////////////////

        //std::cout << "7\n" << std::flush;

        //std::cout << "8\n" << std::flush;
        
        //pivot w0 from which the exponential tail is computed
        double w0 = 1. / (look_up_table_phtn_1_o_w_proj[lookup_index]).evaluate(chi);
        //this variable is cumbersome but this
        // high tail is statistically never called
        double v0 = w0 * w0 * w0;

        //this evaluation brings us to the partial photon emission rate, by integrating
        //the differential probability from 0 to a certain w
        //(keep chi fixed). Please, notice we are evaluating
        //this on the upper bound of validity of the approximation w0
        double g_of_w_chi = (look_up_table_phtn_prtl_rt[lookup_index]).evaluate(chi, w0);

        //std::cout << "9\n" << std::flush;

        //be careful to when v is null
        v = std::abs((tildeWrad - rnd_tildeWrad) / (tildeWrad - g_of_w_chi));
                
        //this is w
        v = cbrt(v0 - log(v));

    }

    //return the photon energy disguised as w value
    return v;
} 

//return the emitted photon energy (normalized in units of m_e c^2)
double Nonlinear_Inverse_Compton_Scattering::SFQED_LCFA_emitted_photon_energy(const double &gamma,
                                                                                const double &chi,
                                                                                const double &rnd) const {


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
    
    double v = SFQED_LCFA_phtn_nrg_aux(chi, rnd, lookup_index);

    //CALCULATE THE PHOTON MOMENTUM AND ENERGY
    //in normalized units (i.e., in units of m_e c^2) by using:
    //( \varepsilon_\gamma / \varepsilon ) = 3 \chi v / (2 + 3 \chi v)
    v = 3.0*chi*v;
    v = gamma * v / (2.0 + v);

    return v;
}

void Nonlinear_Inverse_Compton_Scattering::SFQED_finalize_PHTN_emission(){
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
}


/*********/
/* BLCFA */
/*********/

void Nonlinear_Inverse_Compton_Scattering::SFQED_init_PHTN_emission_BLCFA(std::string path_to_coeffs){
    
    //differential cross section
    std::string emission_name = path_to_coeffs + "dP_over_dt_dw/phtn_diff_crss_sctn_0-2.txt";
    std::ifstream in_file(emission_name.c_str());
    look_up_table_phtn_diff_crss_sctn[4] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "dP_over_dt_dw/phtn_diff_crss_sctn_2-20.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_diff_crss_sctn[3] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "dP_over_dt_dw/phtn_diff_crss_sctn_20-80.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_diff_crss_sctn[2] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "dP_over_dt_dw/phtn_diff_crss_sctn_80-600.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_diff_crss_sctn[1] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "dP_over_dt_dw/phtn_diff_crss_sctn_600-2000.txt";
    in_file.open(emission_name.c_str());
    look_up_table_phtn_diff_crss_sctn[0] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();
}

//\gamma here represents the normalized energy of the particle \epsilon,
//that in case of an electron \epsilon = \gamma. In general it is not
//the sheer gamma factor
//if an increasing tail appears at low energies (log scale) please check the conditions
double Nonlinear_Inverse_Compton_Scattering::SFQED_BLCFA_emitted_photon_energy(const double& LCFA_limit,
                                                            const double& gamma, 
                                                            const double& chi, 
                                                            const double& rnd, 
                                                            const double& rnd2) const{

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
    

    //in case the phtn nrg threshold is zero or very close to the particle nrg
    //or it is above the upper bound of the approximated domain
    //we don't want the emission to occur 
    if(LCFA_limit == 0. || LCFA_limit > 0.75 * gamma || LCFA_limit_w > (lu_table_upper_w_bounds[lookup_index])){
        return 0.;
    }

    //compute the w-value energy of the LCFA emitted photon
    double emitted_phtn_LCFA_nrg_w = SFQED_BLCFA_phtn_nrg_aux(chi, rnd, lookup_index),
            diff_probability,
            diff_prob_LCFA,
            rs_dw_over_depsgamma;

    //remember that the assignement operator return the value
    //being assigned


    //BLCFA section (it starts with a comparison)
    if(emitted_phtn_LCFA_nrg_w < LCFA_limit_w){
        

        //compute the threshold value
        rs_dw_over_depsgamma = (1. + 1.5 * chi * LCFA_limit_w * LCFA_limit_w * LCFA_limit_w) / LCFA_limit_w;
        diff_prob_LCFA = (look_up_table_phtn_diff_crss_sctn[lookup_index]).evaluate(chi, LCFA_limit_w) * rs_dw_over_depsgamma * rs_dw_over_depsgamma;

        //compute the differential probability associated to the LCFA supposed phtn nrg
        //(pay attention not to overwrite the w value energy)
        rs_dw_over_depsgamma = (1. + 1.5 * chi * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w) / emitted_phtn_LCFA_nrg_w;
        diff_probability = (look_up_table_phtn_diff_crss_sctn[lookup_index]).evaluate(chi, emitted_phtn_LCFA_nrg_w) * rs_dw_over_depsgamma * rs_dw_over_depsgamma;

        //now I use the rejection sampling method (employing rnd2) to determine if I have
        //to keep or reject the just occurred emission
        if(diff_probability * rnd2 > diff_prob_LCFA){
            return 0.;
        }

    }

    //convert to a proper photon energy!
    emitted_phtn_LCFA_nrg_w = 3.0 * chi * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w * emitted_phtn_LCFA_nrg_w;
    emitted_phtn_LCFA_nrg_w = gamma * emitted_phtn_LCFA_nrg_w / (2.0 + emitted_phtn_LCFA_nrg_w);

    return emitted_phtn_LCFA_nrg_w;
    ///////////////////////////////////////////////
}


void Nonlinear_Inverse_Compton_Scattering::SFQED_finalize_PHTN_emission_BLCFA(){
    //purge diff crss sctn
    delete[] look_up_table_phtn_diff_crss_sctn[4].last_coeffs;
    delete[] look_up_table_phtn_diff_crss_sctn[3].last_coeffs;
    delete[] look_up_table_phtn_diff_crss_sctn[2].last_coeffs;
    delete[] look_up_table_phtn_diff_crss_sctn[1].last_coeffs;
    delete[] look_up_table_phtn_diff_crss_sctn[0].last_coeffs;
}


/***************************/
/* GENERIC ROUTINE SECTION */
/***************************/
void SFQED_collinear_momentum(const double &gamma_out, const double p_in[3], double p_out[3]){
    
    double v = sqrt( (gamma_out*gamma_out - 1.) / (p_in[0]*p_in[0] + p_in[1]*p_in[1] + p_in[2]*p_in[2]) );

    p_out[0] = v*p_in[0];
    p_out[1] = v*p_in[1];
    p_out[2] = v*p_in[2];

    return;
}