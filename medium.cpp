#include "medium.h"

#include <cmath>

#include "macros.h"
double calc_B_field(double eps_B, double Gamma, double rho) {
    return sqrt(8 * con::pi * eps_B * rho * (4 * Gamma * Gamma - 4 * Gamma)) * con::c;
}

double shock_width_com(double r, double Gamma) { return r / Gamma / 12; }

Medium create_ISM(double n_ism, double eps_e, double eps_B, double pel) {
    auto rho = [=](double r) { return n_ism * con::mp; };
    auto m = [=](double r) { return 4 * con::pi / 3 * r * r * r * n_ism * con::mp; };
    return Medium{rho, m, eps_e, eps_B, pel};
}