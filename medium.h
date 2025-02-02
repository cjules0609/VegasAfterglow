#ifndef _MEDIUM_
#define _MEDIUM_

#include <iostream>

#include "macros.h"
#include "mesh.h"
class Medium {
   public:
    Medium(UnaryFunc rho, UnaryFunc mass, double cs = 1e-5 * con::c) : rho(rho), mass(mass), cs(cs) {};

    double cs{1e-5 * con::c};  // sound speed
    UnaryFunc rho;
    UnaryFunc mass;
};

Medium createISM(double n_ism, double cs = 1e-5 * con::c);
Medium createWind(double n_r, double r_c, double k, double cs = 1e-5 * con::c);

#endif