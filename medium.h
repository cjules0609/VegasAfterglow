#ifndef _MEDIUM_
#define _MEDIUM_

#include <iostream>

#include "mesh.h"

class Medium {
   public:
    Medium(Profile rho, Profile mass, double eps_e, double eps_B, double pel, double xi = 1)
        : rho(rho), mass(mass), eps_e(eps_e), eps_B(eps_B), pel(pel), xi(xi) {
        if (pel < 1) {
            std::cout << "pel should be larger than 1" << std::endl;
            exit(1);
        }
    };
    double eta_rad{1};
    double eps_e{0.1};
    double eps_B{0.01};
    double zeta{1};
    double xi{1};
    double pel{2.3};
    Profile rho;
    Profile mass;
};

Medium create_ISM(double n_ism, double eps_e, double eps_B, double pel);

double calc_B_field(double eps_B, double Gamma, double rho);

double shock_width_com(double r, double Gamma);
#endif