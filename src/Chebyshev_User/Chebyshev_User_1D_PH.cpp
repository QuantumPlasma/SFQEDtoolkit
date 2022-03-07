#include "Chebyshev_User_1D_PH.h"

double Chebyshev_User_1D_PH::evaluate(double x){
    return (this->func_to_execute)(x);
}

