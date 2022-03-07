#include "Chebyshev_User_2D_PH.h"

double Chebyshev_User_2D_PH::evaluate(double x, double y){
    return (this->func_to_execute)(x, y);
}