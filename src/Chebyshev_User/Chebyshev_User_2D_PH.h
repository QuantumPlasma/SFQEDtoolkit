#pragma once
#ifndef CHEBYSHEV_PH_2D
#define CHEBYSHEV_PH_2D

#include "Chebyshev_User_2D.h"

struct Chebyshev_User_2D_PH: Chebyshev_User_2D {

public:

    //call the super class default constructor (thas won't do nothing)
    Chebyshev_User_2D_PH():Chebyshev_User_2D(){

    }

    double (*func_to_execute)(double, double);

    virtual double evaluate(double, double);

};

#endif