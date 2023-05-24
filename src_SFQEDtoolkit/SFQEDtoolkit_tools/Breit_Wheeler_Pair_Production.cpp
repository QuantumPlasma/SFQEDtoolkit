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

#include "Breit_Wheeler_Pair_Production.h"


/***************************************/
/* PAIR PRODUCTION function definition */
/***************************************/

void Breit_Wheeler_Pair_Production::SFQED_init_PAIR_creation(std::string path_to_coeffs){

    // init_pair_crtn_tables();

    //pair creation rate

//
    look_up_table_pair_prd_rt[6].domain_half_range = 1;
    look_up_table_pair_prd_rt[6].domain_middle_point = 1;
    look_up_table_pair_prd_rt[6].evaluation_order = 0;
//

    std::string emission_name = path_to_coeffs + "pair_rate/pair_production_rate_024-04.txt";
    std::ifstream in_file(emission_name.c_str());
    look_up_table_pair_prd_rt[5] = Chebyshev_User_1D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_rate/pair_production_rate_0-2.txt";
    in_file.open(emission_name.c_str());
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
    emission_name = path_to_coeffs + "pair_prod_nrgs/pair_nrgs_v_001-03.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs[5] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();
//

    emission_name = path_to_coeffs + "pair_prod_nrgs/pair_nrgs_v_0-2.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs[4] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_prod_nrgs/pair_nrgs_v_2-20.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs[3] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_prod_nrgs/pair_nrgs_v_20-80.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs[2] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_prod_nrgs/pair_nrgs_v_80-600.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs[1] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_prod_nrgs/pair_nrgs_v_600-2000.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs[0] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();


    //photon energy in v high part
//
    emission_name = path_to_coeffs + "pair_prod_nrgs/pair_nrgs_v_high_001-03.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs_high[5] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();
//

    emission_name = path_to_coeffs + "pair_prod_nrgs/pair_nrgs_v_high_0-2.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs_high[4] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_prod_nrgs/pair_nrgs_v_high_2-20.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs_high[3] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_prod_nrgs/pair_nrgs_v_high_20-80.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs_high[2] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_prod_nrgs/pair_nrgs_v_high_80-600.txt";
    in_file.open(emission_name.c_str());
    look_up_table_pair_v_nrgs_high[1] = Chebyshev_User_2D::init_from_txt_file(in_file);
    in_file.close();
    in_file.clear();

    emission_name = path_to_coeffs + "pair_prod_nrgs/pair_nrgs_v_high_600-2000.txt";
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

double Breit_Wheeler_Pair_Production::SFQED_PAIR_emitted_electron_energy_aux(const int &lookup_index, const double &chi, const double &rescaled_rnd) const {

    const double low_r_limit = (lu_table_pair_lower_r_bounds[lookup_index]),
            high_r_limit = (lu_table_pair_upper_r_bounds[lookup_index]);

    // std::cout << "1\n";

    double v;

    if(rescaled_rnd < low_r_limit){

        // std::cout << "3\n";

        // std::cout << "4\n";

        v = (look_up_table_pair_v_nrgs[lookup_index]).evaluate(chi, rescaled_rnd);

        // std::cout << "6\n";
    }
    else if(rescaled_rnd < high_r_limit){

        // std::cout << "a\n";

        v = (look_up_table_pair_v_nrgs_high[lookup_index]).evaluate(chi, rescaled_rnd);

        // std::cout << "b\n";
    }
    //in the off chance case that we happen to be beyond the
    //high r threshold, we will assume to be dealing directly
    //with v = 0. Since v = (\epsilon_{\gamma} - 2 \epsilon_e) / \epsilon_{\gamma}
    else{

        // std::cout << "7\n";

        //select the correct lookup_index to call the rate
        int lookup_index_rate = lookup_index
                                + (chi < bound_kappa_0th_rate)
                                + (chi >= bound_kappa_0th_nrgs) * (chi < bound_kappa_intermediate_rate);

        // std::cout << "8\n";

        //compute the rate of pair emission
        //introducing a new branching (we are allowed to do this,
        // as this part is called only the 0.0001% of times)
        double tildeWrad = (lookup_index_rate == 6) ?
                                (9. * pigreek * chi)/(16. * sqrt(2.)) * exp(-(8/(3 * chi))) * (1. - 11./64.*chi + 7585./73728.*chi*chi) : 
                                (look_up_table_pair_prd_rt[lookup_index_rate]).evaluate(chi); // / 3.;
                            
        double rnd_tildeWrad = rescaled_rnd * tildeWrad;

        // std::cout << "9\n";

        // std::cout << "10\n";

        //compute the pivot at which the exponential approximation will be applied
        double v0 = (look_up_table_pair_v_nrgs_proj[lookup_index]).evaluate(chi);
        v = 8.0 / (3.0 * chi * (1.0 - v0*v0));

        // std::cout << "11\n";

        // std::cout << "12\n";

        //integrate the pair creation differential probability
        //from the just computed pivot up to 1
        double tail_integral = (look_up_table_pair_prtl_rt[lookup_index]).evaluate(chi, v0);

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
double Breit_Wheeler_Pair_Production::SFQED_PAIR_created_electron_energy(const double &nrg, const double &chi, const double &rnd) const {


    int lookup_index;

    if(chi <= bound_kappa_0th_nrgs){
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

double Breit_Wheeler_Pair_Production::SFQED_PAIR_created_electron_energy_fast(const double &nrg, const double &chi, const double &rnd) const {

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

void Breit_Wheeler_Pair_Production::SFQED_finalize_PAIR_creation(){
    //purge all coefficients
    delete[] look_up_table_pair_prd_rt[5].last_coeffs;
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