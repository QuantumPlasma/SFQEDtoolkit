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
struct SFQED_Processes;

class BLCFA_Object{
protected:
    BLCFA_Object(){
        //the initialization of the pointers has done only once
        Lorentz_F_Old = new double[3]{0., 0., 0.};
        Delta_Lorentz_F_Old = new double[3]{0., 0., 0.};
        just_created = true;
    }

    ~BLCFA_Object(){
        delete [] Lorentz_F_Old;
        delete [] Delta_Lorentz_F_Old;
    }


    // be careful: the following pointer should store the perpendicular
    // component (wrt the direction of the momentum) of the reduced
    // Lorentz force, defined as the time derivative of the reduced
    // momentum pushed by the regular Lorentz pusher
    // \vec{\alpha}
    double* Lorentz_F_Old;

    // this vector stores the difference between the actual reduced 
    // Lorentz force and the previously processed one
    double* Delta_Lorentz_F_Old;
    
    bool just_created;

    double optical_depth_ref, optical_depth;

public:

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

    double* FL_old(){
        return Lorentz_F_Old;
    }

    double* deltaFL_old(){
        return Delta_Lorentz_F_Old;
    }

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

    //this method is used not only to compute and update the Lorentz_F_Old and
    //Delta_Lorentz_F_Old variables, as the name suggests, but also to compute 
    //the force derivatives, the gamma and the chi parameter. Its return value
    //is used to infer whether a particle has just been created
    bool SFQED_BLCFA_update_entities_quantities(const SFQED_Processes&, const double* const, const double* const, double*, double*, double&, double&);

    //The main purpose of this method is to compute the energy threshold
    //below which the LCFA ceases to be valid. It needs the 
    double SFQED_BLCFA_find_energy_threshold(const SFQED_Processes&, const double* const, const double* const, const double&, const double&) const;

    void init_optical_depth(double);
    void reset_optical_depth();
    double increase_optical_depth(double);
    bool check_optical_depth_for_emission();

};

#endif