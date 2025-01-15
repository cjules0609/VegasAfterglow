#include "medium.h"

#include <cmath>

#include "macros.h"

Medium create_ISM(double n_ism, double cs) {
    auto rho = [=](double r) { return n_ism * con::mp; };
    auto m = [=](double r) { return 4 * con::pi / 3 * r * r * r * n_ism * con::mp; };
    return Medium{rho, m, cs};
}

Medium create_wind(double n_r, double r_c, double k, double cs) {
    auto rho = [=](double r) { return n_r * pow(r / r_c, -k) * con::mp; };
    auto m = [=](double r) { return 4 * con::pi / (3 - k) * r * r * r * n_r * pow(r / r_c, -k) * con::mp; };
    return Medium{rho, m, cs};
}