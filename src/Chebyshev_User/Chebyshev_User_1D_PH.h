#pragma once
#ifndef CHEBYSHEV_PH_1D
#define CHEBYSHEV_PH_1D

#include "Chebyshev_User_1D.h"

struct Chebyshev_User_1D_PH: Chebyshev_User_1D {

public:

    //call the super class default constructor (thas won't do nothing)
    Chebyshev_User_1D_PH():Chebyshev_User_1D(){

    }

    double (*func_to_execute)(double);

    virtual double evaluate(double);

};

#endif


