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

#include "BLCFA_Object.h"
#include "SFQED_User_Constants.h"

#include <cmath>

//definition of BLCFA_Object's static variables
double BLCFA_Object::one_over_dt, 
        BLCFA_Object::one_over_dt2,
        BLCFA_Object::norm_Compton_time_2;


bool BLCFA_Object::SFQED_BLCFA_update_entities_quantities(const double* const pushed_momentum,
                                                const double* const momentum,
                                                double& delta,
                                                double& part_gamma,
                                                double& part_chi){

    //all the std outputs are for debug purposes, they can easily be canceled!
    // std::cout << "a\n";

    //tmp variables needed to store the lorentz force derivatives
    double Lorentz_F_t_der[3], Lorentz_F_tt_der[3];
    //value to return
    bool returner;                       

    //even though this tmp array variable is meant to store the
    //force derivative value, we will use it to temporarily deal with
    //the actual Lorentz force vector
    Lorentz_F_t_der[0] = BLCFA_Object::one_over_dt * (pushed_momentum[0] - momentum[0]),
    Lorentz_F_t_der[1] = BLCFA_Object::one_over_dt * (pushed_momentum[1] - momentum[1]),
    Lorentz_F_t_der[2] = BLCFA_Object::one_over_dt * (pushed_momentum[2] - momentum[2]);

    // std::cout << "b\n";

    //sum of the pushed and unpushed momenta
    double momentum_sum_x = (pushed_momentum[0] + momentum[0]),
            momentum_sum_y = (pushed_momentum[1] + momentum[1]),
            momentum_sum_z = (pushed_momentum[2] + momentum[2]);

    // std::cout << "c\n";

    //compute the proper gamma factor. we use this to approximate the momentum modulus
    //and avoid possible numerical issue. This is not the actual gamma, but 4 * gamma^2.
    part_gamma = 4. + (momentum_sum_x*momentum_sum_x + momentum_sum_y*momentum_sum_y + momentum_sum_z*momentum_sum_z);

    // std::cout << "d\n";

    //this variable is called chi, but right now it stores the
    //projection factor of the Lorentz force along the 
    //momentum_sum vector direction
    part_chi = (Lorentz_F_t_der[0]*momentum_sum_x + Lorentz_F_t_der[1]*momentum_sum_y + Lorentz_F_t_der[2]*momentum_sum_z) 
                    / part_gamma;

    // std::cout << "d\n";

    //compute the proper gamma [we also update a variable passed from outside]
    part_gamma = std::sqrt(part_gamma*0.25);

    // std::cout << "e\n";

    //now we replace the actual Lorentz force comps with those
    //perpendicular to the momentum
    Lorentz_F_t_der[0] = Lorentz_F_t_der[0] - part_chi * momentum_sum_x,
    Lorentz_F_t_der[1] = Lorentz_F_t_der[1] - part_chi * momentum_sum_y,
    Lorentz_F_t_der[2] = Lorentz_F_t_der[2] - part_chi * momentum_sum_z;

    // std::cout << "f\n";

    //compute the squared modulus of the orthogonal Lorentz force
    //multiplied by the squared gamma factor
    double L_F_Modulus = Lorentz_F_t_der[0]*Lorentz_F_t_der[0] +
                                Lorentz_F_t_der[1]*Lorentz_F_t_der[1] + 
                                Lorentz_F_t_der[2]*Lorentz_F_t_der[2];

    // std::cout << "h\n";

    //we now use the previously defined momentum sum components to store
    //the new orthogonal lorentz force difference
    momentum_sum_x = Lorentz_F_t_der[0] - this->Lorentz_F_Old[0],
    momentum_sum_y = Lorentz_F_t_der[1] - this->Lorentz_F_Old[1],
    momentum_sum_z = Lorentz_F_t_der[2] - this->Lorentz_F_Old[2];

    // std::cout << "i\n";

    //update the old lorentz force quantities. IMPORTANT:
    //before updating this variable be sure to have calculated the
    //new orthogonal Lorentz force difference
    this->Lorentz_F_Old[0] = Lorentz_F_t_der[0];
    this->Lorentz_F_Old[1] = Lorentz_F_t_der[1];
    this->Lorentz_F_Old[2] = Lorentz_F_t_der[2];

    // std::cout << "j\n";

    // //[IMPORTANT!!!]
    // //if the particle has just been created we exit
    // //immediately after we updated the FORCES, GAMMA and CHI! 
    // //and return false (the emission is not possible). But before
    // //actually returning...
    // if(this->just_created){
    //     //... switch the just created particle flag to false
    //     this->just_created = false;

    //     //... set the delta variable (used to extract the LCFA threshold)
    //     // to a pointless value (-1)
    //     delta = -1.;

    //     //... update the old lorentz force difference.
    //     //(remember that the lorentz force has already
    //     //been updated at this point, and we do not need
    //     //to compute the second lorentz derivative).
    //     this->Delta_Lorentz_F_Old[0] = momentum_sum_x;
    //     this->Delta_Lorentz_F_Old[1] = momentum_sum_y;
    //     this->Delta_Lorentz_F_Old[2] = momentum_sum_z;

    //     return false; 
    // }

    //IMPORTANT:
    //if the particles has been emitted at some previous step
    //then just_created is false and we continue to the calculation
    //of the delta variable

    // std::cout << "k\n";
    // std::cout << proc.one_over_dt2 << ' ' << this->Delta_Lorentz_F_Old[0] << ' '
    //             << this->Delta_Lorentz_F_Old[1] << ' ' << this->Delta_Lorentz_F_Old[2] << '\n'; 

    //compute second derivatives
    Lorentz_F_tt_der[0] = BLCFA_Object::one_over_dt2 * (momentum_sum_x - this->Delta_Lorentz_F_Old[0]),
    Lorentz_F_tt_der[1] = BLCFA_Object::one_over_dt2 * (momentum_sum_y - this->Delta_Lorentz_F_Old[1]),
    Lorentz_F_tt_der[2] = BLCFA_Object::one_over_dt2 * (momentum_sum_z - this->Delta_Lorentz_F_Old[2]);

    // std::cout << "l\n";

    //update old lorentz force difference. IMPORTANT:
    //before updating this variable be sure to have calculated the
    //second Lorentz force derivatives
    this->Delta_Lorentz_F_Old[0] = momentum_sum_x;
    this->Delta_Lorentz_F_Old[1] = momentum_sum_y;
    this->Delta_Lorentz_F_Old[2] = momentum_sum_z;

    //at this point all the BLCFA_Object member variables has been updated

    // std::cout << "m\n";

    //compute the actual first derivatives
    Lorentz_F_t_der[0] = BLCFA_Object::one_over_dt * momentum_sum_x;
    Lorentz_F_t_der[1] = BLCFA_Object::one_over_dt * momentum_sum_y;
    Lorentz_F_t_der[2] = BLCFA_Object::one_over_dt * momentum_sum_z;

    // std::cout << "g\n";

    //compute the ultimate chi [we update also a variable passed from outside]
    part_chi = part_gamma * std::sqrt(BLCFA_Object::norm_Compton_time_2 * L_F_Modulus);//* (mass / charge); we decided to provide these methods for positrons and electrons only


    //at this stage all the derivatives are computed and we can proceed
    //to the evaluation of the delta variable and its test

    // std::cout << "n\n";

    //first of all: compute the delta (see article)
    delta = BLCFA_Object::norm_Compton_time_2 * 
                        (Lorentz_F_t_der[0]*Lorentz_F_t_der[0] + Lorentz_F_t_der[1]*Lorentz_F_t_der[1] + Lorentz_F_t_der[2]*Lorentz_F_t_der[2] +
                         std::abs(this->Lorentz_F_Old[0]*Lorentz_F_tt_der[0] + this->Lorentz_F_Old[1]*Lorentz_F_tt_der[1] + this->Lorentz_F_Old[2]*Lorentz_F_tt_der[2]));

    // std::cout << "o\n";

    //this variable stores the quantity (\chi^2)/(\gamma^2)|F_{\perp}|^2
    //which we know to be equal to \tau_C^2 |F_{\perp}|^4 (needed in the control below),
    //since \tau_C = (\chi)/(\gamma) * 1 / |F_{\perp}|
    double chi_2_by_F_2_over_gamma_2 = BLCFA_Object::norm_Compton_time_2 * L_F_Modulus * L_F_Modulus;

    // std::cout << "p\n";

    //the first condition is needed to check if chi is big enough
    //while the second keeps us from encountering numerical issues
    //(like division by 0 or similar) when the fields are constant
    //and the force derivatives are zero.
    returner = (!(this->just_created) && part_chi > zeta && delta / zeta_2 > chi_2_by_F_2_over_gamma_2);
    // returner = (!(this->just_created));
    // std::cout << "q\n";

    this->just_created = false;

    return returner;
}


double BLCFA_Object::SFQED_BLCFA_find_energy_threshold(const double& delta,
                                                        const double& part_gamma,
                                                        const double& part_chi) const {
    
    //remember that this represents the squared modulus
    //of the Lorentz force perpendicular to the momentum
    double L_F_Modulus_Sq = this->Lorentz_F_Old[0]*this->Lorentz_F_Old[0]
                            + this->Lorentz_F_Old[1]*this->Lorentz_F_Old[1]
                            + this->Lorentz_F_Old[2]*this->Lorentz_F_Old[2];

    //this variable stores the quantity (\chi^2)/(\gamma^2)|F_{\perp}|^2
    //which we know to be equal to \tau_C^2 |F_{\perp}|^4 (needed in the control below),
    //since \tau_C = (\chi)/(\gamma) * 1 / |F_{\perp}|
    double formation_time_ratio_mult_chi_over_gamma = twopi * L_F_Modulus_Sq * std::sqrt(BLCFA_Object::norm_Compton_time_2 / delta);

    //debug
    // std::cout << delta << " " << part_gamma << " " << part_chi << " " << formation_time_ratio_mult_chi_over_gamma << " ";

    //this is the normalized LCFA energy threshold
    //[if we want this method to hold for a generic
    // particle of mass m, we should multiply the
    // following expression by m (remember that an electron has m = 1)]
    double LCFA_limit = thresh_factor * part_gamma * part_chi / (part_chi + four_over_threepi * sinh(3. * asinh(formation_time_ratio_mult_chi_over_gamma / 8.)));

    return LCFA_limit;
}


//we expecti this functions to be part of the particle entity
// void BLCFA_Object::init_optical_depth(double tau_0){
//     this->optical_depth_ref = tau_0;
//     this->optical_depth = 0.0;
// }

// void BLCFA_Object::reset_optical_depth(){
//     this->optical_depth = 0.0;
// }

// double BLCFA_Object::increase_optical_depth(double increment){
//     return this->optical_depth += increment;
// }

// bool BLCFA_Object::check_optical_depth_for_emission(){
//     return this->optical_depth >= this->optical_depth_ref;
// }