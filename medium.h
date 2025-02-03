#ifndef _MEDIUM_
#define _MEDIUM_

#include <iostream>

#include "macros.h"
#include "mesh.h"
#include "utilities.h"
class Medium {
   public:
    Medium(double k, double n_c, double r_c);
    Medium(double n_c);

    double rho(double r) const;
    double mass(double r) const;

   private:
    double const k;
    double const n_c;
    double const r_c;
};

inline double Medium::rho(double r) const {
    if (k == 0) {
        return n_c * con::mp;
    } else {
        return n_c * fastPow(r / r_c, -k) * con::mp;
    }
}

inline double Medium::mass(double r) const {
    if (k == 0) {
        return 4 * con::pi / 3 * r * r * r * n_c * con::mp;
    } else if (k == 3) {
        return 4 * con::pi * n_c * r_c * r_c * r_c * fastLog(r / con::cm) * con::mp;
    } else {
        return 4 * con::pi / (3 - k) * r * r * r * n_c * fastPow(r / r_c, -k) * con::mp;
    }
}

Medium createISM(double n_ism);

#endif