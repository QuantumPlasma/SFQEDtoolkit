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


#pragma once
#ifndef BLCFA_OBJ
#define BLCFA_OBJ

#include <iostream>

//forward definition of SFQED_Processes
// struct SFQED_Processes;

class BLCFA_Object{
private:

protected:

    // be careful: the following pointer should store the perpendicular
    // component (wrt the direction of the momentum) of the reduced
    // Lorentz force, defined as the time derivative of the reduced
    // momentum pushed by the regular Lorentz pusher
    // \vec{\alpha}
    double Lorentz_F_Old[3];

    // this vector stores the difference between the actual reduced 
    // Lorentz force and the previously processed one
    double Delta_Lorentz_F_Old[3];
    
    bool just_created;

    // double optical_depth_ref, optical_depth;

public:

    //declaration of static member variables
    //(remember you need to also define them outside)
    // the definition of the static variables has been put inside the SFQED_Processes source file
    static double one_over_dt, one_over_dt2, norm_Compton_time_2;

    //useful for debug
    // void* get_F(){
    //     return (void*)Lorentz_F_Old;
    // }

    // void* get_delta_F(){
    //     return (void*)Delta_Lorentz_F_Old;
    // }

    // std::string print_forces(){
    //     return "olf F = " + std::to_string(Lorentz_F_Old[0]) + ", " + std::to_string(Lorentz_F_Old[1]) + ", " + std::to_string(Lorentz_F_Old[2]) + "\n"
    //             + "olf delta_F = " + std::to_string(Delta_Lorentz_F_Old[0]) + ", " + std::to_string(Delta_Lorentz_F_Old[1]) + ", " + std::to_string(Delta_Lorentz_F_Old[2]) + "\n";
    // }

    BLCFA_Object(){
        //the initialization of the pointers has done only once
        Lorentz_F_Old[0] = 0.;
        Lorentz_F_Old[1] = 0.;
        Lorentz_F_Old[2] = 0.;
        
        Delta_Lorentz_F_Old[0] = 0.;
        Delta_Lorentz_F_Old[1] = 0.;
        Delta_Lorentz_F_Old[2] = 0.;
        
        just_created = true;
    }

    ~BLCFA_Object(){
        // delete [] Lorentz_F_Old;
        // delete [] Delta_Lorentz_F_Old;
    }

    bool is_just_created(){
        return just_created;
    }

    double* FL_old(){
        return Lorentz_F_Old;
    }

    double* deltaFL_old(){
        return Delta_Lorentz_F_Old;
    }

    //this function set both the time step related
    // variables and the squared compton time. The 
    // Compton time passed as second argument must be
    // in normalized units.
    static void SFQED_BLCFA_set_time_constants(const double& dt, const double& norm_compton_time){
        one_over_dt = 1. / dt;
        one_over_dt2 = one_over_dt * one_over_dt;
        norm_Compton_time_2 = norm_compton_time * norm_compton_time;
    }

    //this method is used not only to compute and update the Lorentz_F_Old and
    //Delta_Lorentz_F_Old variables, as the name suggests, but also to compute 
    //the delta variable, the gamma and the chi parameter. Its return value
    //is used to infer whether a particle can emit a photon or not. The delta is
    //a variable needed to compute the LCFA validity threshold.
    bool SFQED_BLCFA_update_entities_quantities(const double* const pushed_momentum,
                                                const double* const momentum,
                                                double& delta,
                                                double& part_gamma,
                                                double& part_chi);

    bool SFQED_BLCFA_update_entities_quantities_DEBUG(const double* const pushed_momentum,
                                                const double* const momentum,
                                                double* F,
                                                double* d_F,
                                                double* dd_F,
                                                double& delta,
                                                double& part_gamma,
                                                double& part_chi);

    //The main purpose of this method is to compute the energy threshold
    //below which the LCFA ceases to be valid. It needs the delta variable,
    // gamma and chi parameter computed through the function 
    // SFQED_BLCFA_update_entities_quantities
    double SFQED_BLCFA_find_energy_threshold(const double& delta,
                                                const double& part_gamma, 
                                                const double& part_chi) const;

    //****************************************************************************************
    //the following method are not really needed, as in order
    // to update the BLCFA_Object's members one should resort to  
    // the function SFQED_BLCFA_update_entities_quantitiesupdate

    void update_force(const double* const force){
        Lorentz_F_Old[0] = force[0];
        Lorentz_F_Old[1] = force[1];
        Lorentz_F_Old[2] = force[2];
    }

    void update_delta_force(const double* const delta_force){
        Delta_Lorentz_F_Old[0] = delta_force[0];
        Delta_Lorentz_F_Old[1] = delta_force[1];
        Delta_Lorentz_F_Old[2] = delta_force[2];
    }

};

//////////////////////////////////////////
// FUNCTIONS TO BYPASS THE BLCFA_OBJECT //
//////////////////////////////////////////

//This function offers the possibility to update the Force and
// the difference of the force experienced by the particles at previous timestep
// without resorting to the BLCFA_Object. Be aware that by bypassing the
// BLCFA_Object circuit you are in charge of defining all the variables you need,
// that is two double precision vectors, playing the role of the Force and 
// difference of the force at the previous timestep, and a boolean carrying
// the info about the applicability of the BLCFA method. IMPORTANT: these 3 
// variables, after being used, will be updated with new values everytime
// this function is called. As a byproduct, also the delta variable (which is needed
// to compute the LCFA validity threshold(), the gamma and the chi parameter
// are computed. Its return value is used to infer whether a particle can emit a photon or not.
bool SFQED_BLCFA_update_quantities_raw(const double* const pushed_momentum,
                                                const double* const momentum,
                                                double* const Lorentz_F_Old,
                                                double* const Delta_Lorentz_F_Old,
                                                bool& just_created,
                                                double& delta,
                                                double& part_gamma,
                                                double& part_chi);

//This method returns the energy threshold below which the LCFA ceases to be valid.
// It needs the Lorentz force that's being experienced by the particle, the delta
// variable, the gamma (Lorentz factor) and chi parameter. All these quantities are
// computed through the update function SFQED_BLCFA_update_quantities_raw.
double SFQED_BLCFA_find_energy_threshold_raw(const double* const Lorentz_F,
                                                const double& delta,
                                                const double& part_gamma,
                                                const double& part_chi);

#endif