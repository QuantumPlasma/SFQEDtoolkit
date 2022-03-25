#include <iostream>
#include <fstream>

#include "Particle3D.h"
#include "PlaneWave3D.h"

using namespace std;

double RelativisticRungeKuttaSolver_LCFA(ofstream&, Particle3D&, PlaneWave3D, double, int, SFQED_Processes&, double (*)(void));

double RelativisticRungeKuttaSolver_BLCFA(ofstream&, Particle3D&, PlaneWave3D, double, int, SFQED_Processes&, double (*)(void));

double RelativisticRungeKuttaSolverPhotons(ofstream&, Particle3D&, PlaneWave3D, double, int, SFQED_Processes&, double (*)(void));