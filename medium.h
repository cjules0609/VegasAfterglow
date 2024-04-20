#ifndef _MEDIUM_
#define _MEDIUM_

#include <iostream>

#include "mesh.h"
class Medium {
   public:
    Medium(Profile rho, Profile mass, double eps_e, double eps_B, double xi = 1)
        : rho(rho), mass(mass), eps_e(eps_e), eps_B(eps_B), xi(xi){};
    double eps_e{0.1};
    double eps_B{0.01};
    double zeta{1};
    double xi{1};
    Profile rho;
    Profile mass;
};

Medium create_ISM(double n_ism, double eps_e, double eps_B);

#endif