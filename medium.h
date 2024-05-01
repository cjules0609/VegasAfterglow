#ifndef _MEDIUM_
#define _MEDIUM_

#include <iostream>

#include "macros.h"
#include "mesh.h"
class Medium {
   public:
    Medium(Profile rho, Profile mass, double eps_e, double eps_B, double xi = 1, double cs = 1e-5 * con::c)
        : rho(rho), mass(mass), eps_e(eps_e), eps_B(eps_B), xi(xi), cs(cs){};
    double eps_e{0.1};
    double eps_B{0.01};
    double zeta{1};
    double xi{1};
    double cs{1e-5 * con::c};  // sound speed
    Profile rho;
    Profile mass;
};

Medium create_ISM(double n_ism, double eps_e, double eps_B, double xi = 1, double cs = 1e-5 * con::c);

#endif