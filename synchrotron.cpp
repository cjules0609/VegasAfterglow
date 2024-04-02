
#include "synchrotron.h"

#include <cmath>

#include "macros.h"
#include "medium.h"
#include "mesh.h"
#include "utilities.h"

SynRadMesh createSynRadGrid(size_t theta_size, size_t r_size, SynRad val = {0, 0, 0, 0, 0}) {
    return SynRadMesh(theta_size, SynRadArray(r_size, val));
}

double SynRad::I_nu(double nu, double pel) const { return I_nu_peak * I_nu_(nu, pel); }

double SynRad::I_nu_(double nu, double pel) const {
    if (nu_a < nu_m && nu_m < nu_c) {
        if (nu <= nu_a) {
            return pow(nu_a / nu_m, 1.0 / 3) * pow(nu / nu_a, 2);
        } else if (nu <= nu_m) {
            return pow(nu / nu_m, 1.0 / 3);
        } else if (nu <= nu_c) {
            return pow(nu / nu_m, -(pel - 1) / 2);
        } else if (nu <= nu_M) {
            return pow(nu_c / nu_m, -(pel - 1) / 2) * pow(nu / nu_c, -pel / 2);
        } else {
            return pow(nu_c / nu_m, -(pel - 1) / 2) * pow(nu_M / nu_c, -pel / 2) * sqrt(nu / nu_M) * exp(-nu / nu_M);
        }
    } else if (nu_m < nu_a && nu_a < nu_c) {
        if (nu <= nu_m) {
            return pow(nu_m / nu_a, (pel + 4) / 2) * pow(nu / nu_m, 2);
        } else if (nu <= nu_a) {
            return pow(nu_a / nu_m, -(pel - 1) / 2) * pow(nu / nu_a, 5.0 / 2);
        } else if (nu <= nu_c) {
            return pow(nu / nu_m, -(pel - 1) / 2);
        } else if (nu <= nu_M) {
            return pow(nu_c / nu_m, -(pel - 1) / 2) * pow(nu / nu_c, -pel / 2);
        } else {
            return pow(nu_c / nu_m, -(pel - 1) / 2) * pow(nu_M / nu_c, -pel / 2) * sqrt(nu / nu_M) * exp(-nu / nu_M);
        }
    } else if (nu_a < nu_c && nu_c < nu_m) {
        if (nu <= nu_a) {
            return pow(nu_a / nu_c, 1.0 / 3) * pow(nu / nu_a, 2);
        } else if (nu <= nu_c) {
            return pow(nu / nu_c, 1.0 / 3);
        } else if (nu <= nu_m) {
            return pow(nu / nu_c, -1 / 2);
        } else if (nu <= nu_M) {
            return pow(nu_m / nu_c, -1 / 2) * pow(nu / nu_m, -pel / 2);
        } else {
            return pow(nu_m / nu_c, -1 / 2) * pow(nu_M / nu_m, -pel / 2) * sqrt(nu / nu_M) * exp(-nu / nu_M);
        }
    } else if (nu_c < nu_a && nu_a < nu_m) {
        if (nu <= nu_a) {
            return pow(nu / nu_a, 2);
        } else if (nu <= nu_m) {
            return pow(nu / nu_a, -1.0 / 2) * pow(nu_c / nu_a, 1.0 / 2) / 3;
        } else if (nu <= nu_M) {
            return pow(nu_m / nu_a, -1.0 / 2) * pow(nu / nu_m, -pel / 2) * pow(nu_c / nu_a, 1.0 / 2) / 3;
        } else {
            return pow(nu_m / nu_a, -1.0 / 2) * pow(nu_M / nu_m, -pel / 2) * pow(nu_c / nu_a, 1.0 / 2) / 3 *
                   sqrt(nu / nu_M) * exp(-nu / nu_M);
        }
    } else if (nu_m < nu_c && nu_c < nu_a) {
        if (nu <= nu_a) {
            return pow(nu / nu_a, 2);
        } else if (nu <= nu_M) {
            return (pel - 1) / 3 * pow(nu / nu_a, -pel / 2) * pow(nu_m / nu_a, (pel - 1) / 2) *
                   pow(nu_c / nu_a, 1.0 / 2);
        } else {
            return (pel - 1) / 3 * pow(nu_M / nu_a, -pel / 2) * pow(nu_m / nu_a, (pel - 1) / 2) *
                   pow(nu_c / nu_a, 1.0 / 2) * sqrt(nu / nu_M) * exp(-nu / nu_M);
        }
    } else if (nu_c < nu_m && nu_m < nu_a) {
        if (nu <= nu_a) {
            return pow(nu / nu_a, 2);
        } else if (nu <= nu_M) {
            return pow(nu / nu_a, -pel / 2) * pow(nu_m / nu_a, (pel - 1) / 2) * pow(nu_c / nu_a, 1.0 / 2) / 3;
        } else {
            return pow(nu_M / nu_a, -pel / 2) * pow(nu_m / nu_a, (pel - 1) / 2) * pow(nu_c / nu_a, 1.0 / 2) / 3 *
                   sqrt(nu / nu_M) * exp(-nu / nu_M);
        }
    }
}

inline double shock_width_com(double r, double Gamma) { return r / Gamma / 12; }

inline double calc_B_field(double eps_B, double Gamma, double rho) {
    return sqrt(8 * con::pi * eps_B * rho * (4 * Gamma * Gamma - 4 * Gamma)) * con::c;
}

inline double calc_syn_I_nu_peak(double r, double Gamma, double Bprime, double rho, double xi, double pel) {
    return (pel - 1) / 2 * sqrt(3) * con::e * con::e * con::e * Bprime / (con::me * con::c2) * 4 * Gamma *
           (rho / con::mp) * xi * shock_width_com(r, Gamma);
}

inline double syn_nu(double gamma, double B) {
    double nu = 3 / (4 * con::pi) * con::e * B / (con::me * con::c) * gamma * gamma;
    return nu;
}

inline double calc_syn_gamma_M(double Bprime, double zeta, double Y_tilt) {
    return sqrt(6 * con::pi * con::e / (con::sigmaT * Bprime * zeta * (1 + Y_tilt)));
}

inline double calc_syn_gamma_m(double Gamma, double gamma_M, double eps_e, double xi, double pel) {
    double gamma_m = 1;
    if (pel > 2) {
        gamma_m = (pel - 2) / (pel - 1) * eps_e * (Gamma - 1) * 1836 / xi + 1;
    } else if (pel < 2) {
        gamma_m =
            pow((2 - pel) / (pel - 1) * eps_e * (Gamma - 1) * 1836 / xi * pow(gamma_M, pel - 1), 1 / (pel - 1)) + 1;
    } else {
        gamma_m = root_bisection(
            [=](double x) -> double {
                return (x * log(gamma_M) - (x + 1) * log(x) - eps_e * (Gamma - 1) * 1836 / xi - log(gamma_M));
            },
            1, gamma_M);
    }
    return gamma_m;
}

inline double calc_syn_gamma_c(double Gamma, double t_com, double Bprime, double Y_tilt) {
    double beta = sqrt(Gamma * Gamma - 1) / Gamma;
    /*double gamma_c = 6 * consts::pi * consts::me * consts::c /
                     (consts::sigmaT * t_com * beta * beta * Bprime * Bprime * (1 + Y_tilt));*/

    // work around to avoid gamma_c < 1
    double gamma_c = 6 * con::pi * con::me * con::c / (con::sigmaT * t_com * Bprime * Bprime * (1 + Y_tilt)) + 1;

    return gamma_c;
}

inline double calc_syn_nu_a(double Gamma, double Bprime, double I_syn_peak, double gamma_m, double gamma_c,
                            double gamma_M) {
    double gamma_peak = std::min(gamma_m, gamma_c);
    double nu_peak = syn_nu(gamma_peak, Bprime);

    double gamma_eos = (4 * Gamma + 1) / (3 * Gamma);  // adiabatic index
    double kT = (gamma_peak - 1) * con::me * con::c2 * (gamma_eos - 1);
    double nu_a = pow(I_syn_peak * con::c2 / pow(nu_peak, 1.0 / 3) / kT / 2,
                      3.0 / 5);  // 2kT(nv_a/c)^2 = I_peak*(nu_a/nu_peak)^(1/3)

    // the nu_peak is not the real peak, peak at nu_a; kT = (gamma_a-1) * me *c^2*(gamma_eos-1), I_syn = I_peak;
    if (nu_a > nu_peak) {
        /*nu_a = pow(I_syn_peak / con::me / 2 / (gamma_eos - 1) /
                       sqrt(4 * con::pi / 3 * con::me * con::c / con::e / Bprime),
                   2.0 / 5);*/ //this works only for gamma >> 1
        double nu_M = syn_nu(gamma_M, Bprime);
        double A = sqrt(4 * con::pi / 3 * con::me * con::c / con::e / Bprime);
        double B = I_syn_peak / (2 * con::me * (gamma_eos - 1));
        nu_a = root_bisection([=](double x) -> double { return A * x * x * x * x * x - x * x * x * x - B; },
                              sqrt(nu_peak), sqrt(nu_M));
        nu_a *= nu_a;
    }
    return nu_a;
}

SynRad syn_radiation(double r, double Gamma, double t_com, Medium const& medium) {
    double rho = medium.rho(r);
    double Bprime = calc_B_field(medium.eps_B, Gamma, rho);
    double I_nu_peak = calc_syn_I_nu_peak(r, Gamma, Bprime, rho, medium.xi, medium.pel);

    double gamma_M = calc_syn_gamma_M(Bprime, medium.zeta, medium.Y_tilt);
    double nu_M = syn_nu(gamma_M, Bprime);

    double gamma_m = calc_syn_gamma_m(Gamma, gamma_M, medium.eps_e, medium.xi, medium.pel);
    double nu_m = syn_nu(gamma_m, Bprime);

    double gamma_c = calc_syn_gamma_c(Gamma, t_com, Bprime, medium.Y_tilt);
    double nu_c = syn_nu(gamma_c, Bprime);

    double nu_a = calc_syn_nu_a(Gamma, Bprime, I_nu_peak, gamma_m, gamma_c, gamma_M);
    return {I_nu_peak, nu_m, nu_c, nu_a, nu_M};
}

SynRadMesh calc_syn_radiation(Coord const& coord, MeshGrid const& Gamma, MeshGrid const& t_com, Medium const& medium) {
    SynRadMesh syn_rad = createSynRadGrid(coord.theta.size(), coord.r.size());
    for (size_t j = 0; j < coord.theta.size(); ++j) {
        for (size_t k = 0; k < coord.r.size(); ++k) {
            syn_rad[j][k] = syn_radiation(coord.r[k], Gamma[j][k], t_com[j][k], medium);
        }
    }
    return syn_rad;
}