#include "inverse-compton.h"

#include <cmath>
#include <iostream>

#include "macros.h"
#include "synchrotron.h"
#include "utilities.h"
ICRadMesh createICRadGrid(size_t theta_size, size_t r_size, ICRad val) {
    return ICRadMesh(theta_size, ICRadArray(r_size, val));
}
inline const double IC_x0 = sqrt(2) / 3;

double ICRad::I_nu(double nu) const { return I_nu_peak * I_nu_(nu); }

inline bool order(double a, double b, double c) { return a < b && b < c; };

double compton_sigma(double nu) {
    double x = con::h * nu / (con::me * con::c2);
    return 0.75 * con::sigmaT *
           ((1 + x) / (x * x * x) * (2 * x * (1 + x) / (1 + 2 * x) - log(1 + 2 * x)) + log(1 + 2 * x) / (2 * x) -
            (1 + 3 * x) / (1 + 2 * x) / (1 + 2 * x));
}

double ICRad::I_nu_(double nu) const {
    if (order(nu_ma, nu_mm, nu_mc)) {
        if (nu < nu_ma) {
            return 5.0 / 2 * (pel - 1) / (pel + 1) * pow(nu_ma / nu_mm, 1.0 / 3) * (nu / nu_ma);
        } else if (nu < nu_mm) {
            return 3.0 / 2 * (pel - 1) / (pel - 1.0 / 3) * pow(nu / nu_mm, 1.0 / 3);
        } else if (nu < nu_mc) {
            return (pel - 1) / (pel + 1) * pow(nu / nu_mm, (1 - pel) / 2) *
                   (4 * (pel + 1.0 / 3) / (pel + 1) / (pel - 1.0 / 3) + log(nu / nu_mm));
        } else if (nu < nu_cc) {
            return (pel - 1) / (pel + 1) * pow(nu / nu_mm, (1 - pel) / 2) *
                   (2 * (2 * pel + 3) / (pel + 2) - 2 / (pel + 1) / (pel + 2) + log(nu_cc / nu));
        } else {
            return (pel - 1) / (pel + 1) * pow(nu / nu_mm, -pel / 2) * (nu_mc / nu_mm) *
                   (2 * (2 * pel + 3) / (pel + 2) - 2 / (pel + 2) / (pel + 2) +
                    (pel + 1) / (pel + 2) * log(nu / nu_cc));
        }
    } else if (order(nu_mm, nu_ma, nu_mc)) {
        if (nu < nu_ma) {
            return 2 * (pel + 4) * (pel - 1) / 3 / (pel + 1) / (pel + 1) * pow(nu_mm / nu_ma, (pel + 1) / 2) *
                   (nu / nu_mm);
        } else if (nu < nu_mc) {
            return (pel - 1) / (pel + 1) * pow(nu / nu_mm, (1 - pel) / 2) *
                   (2 * (2 * pel + 5) / (pel + 1) / (pel + 4) + log(nu / nu_ma));
        } else if (nu < nu_ca) {
            return (pel - 1) / (pel + 1) * pow(nu / nu_mm, (1 - pel) / 2) * (2 + 2 / (pel + 4) + log(nu / nu_mc));
        } else if (nu < nu_cc) {
        } else {
        }
    } else if (order(nu_ma, nu_mc, nu_mm)) {
        if (nu < nu_ca) {
        } else if (nu < nu_cc) {
        } else if (nu < nu_cm) {
        } else if (nu < nu_mm) {
        } else {
        }
    } else if (order(nu_mc, nu_ma, nu_mm)) {
        if (nu < nu_aa) {
        } else if (nu < nu_am) {
        } else if (nu < nu_mm) {
        } else {
        }
    } else if (order(nu_mm, nu_mc, nu_ma)) {
        if (nu < nu_aa) {
        } else {
        }
    } else if (order(nu_mc, nu_mm, nu_ma)) {
        if (nu < nu_aa) {
        } else {
        }
    }
}

inline double calc_IC_I_nu_peak(double I_syn_peak, double r, double Gamma, double rho, double xi, double x0) {
    double tau_es = con::sigmaT * (rho / con::mp) * xi * shock_width_com(r, Gamma) * Gamma;

    return tau_es * I_syn_peak * x0;
}

double IC_nu(double gamma, double nu, double x0) { return 4 * gamma * gamma * nu * x0; }

double IC_nu_E_peak(ICRad const& rad) {
    if (order(rad.nu_ma, rad.nu_mm, rad.nu_mc)) {
        return rad.nu_cc;
    } else if (order(rad.nu_mm, rad.nu_ma, rad.nu_mc)) {
        return rad.nu_cc;
    } else if (order(rad.nu_ma, rad.nu_mc, rad.nu_mm)) {
        return rad.nu_mm;
    } else if (order(rad.nu_mc, rad.nu_ma, rad.nu_mm)) {
        if (rad.nu_ma * rad.nu_ma < rad.nu_mm * rad.nu_mc) {
            return rad.nu_mm;
        } else {
            return rad.nu_aa;
        }
    } else if (order(rad.nu_mm, rad.nu_mc, rad.nu_ma)) {
        return rad.nu_aa;
    } else if (order(rad.nu_mc, rad.nu_mm, rad.nu_ma)) {
        return rad.nu_aa;
    }
}

ICRad IC_radiation(double r, double Gamma, double B, SynRad const& electron, SynRad const& photon,
                   Medium const& medium) {
    double x0 = IC_x0;
    double rho = medium.rho(r);
    double gamma_m = syn_gamma(electron.nu_m, B);
    double gamma_c = syn_gamma(electron.nu_c, B);
    double gamma_a = syn_gamma(electron.nu_a, B);

    ICRad rad;

    rad.I_nu_peak = calc_IC_I_nu_peak(photon.I_nu_peak, r, Gamma, rho, medium.xi, x0);

    rad.nu_mm = IC_nu(gamma_m, photon.nu_m, x0);
    rad.nu_mc = IC_nu(gamma_m, photon.nu_c, x0);
    rad.nu_ma = IC_nu(gamma_m, photon.nu_a, x0);

    rad.nu_cm = IC_nu(gamma_c, photon.nu_m, x0);
    rad.nu_cc = IC_nu(gamma_c, photon.nu_c, x0);
    rad.nu_ca = IC_nu(gamma_c, photon.nu_a, x0);

    rad.nu_am = IC_nu(gamma_a, photon.nu_m, x0);
    rad.nu_ac = IC_nu(gamma_a, photon.nu_c, x0);
    rad.nu_aa = IC_nu(gamma_a, photon.nu_a, x0);

    rad.nu_E_peak = IC_nu_E_peak(rad);

    rad.pel = electron.pel;

    return rad;
}

ICRadMesh calc_IC_radiation(Coord const& coord, MeshGrid const& Gamma, MeshGrid B, SynRadMesh const& electrons,
                            SynRadMesh const& photons, Medium const& medium) {
    ICRadMesh IC_rad = createICRadGrid(coord.theta.size(), coord.r.size());
    for (size_t j = 0; j < coord.theta.size(); ++j) {
        for (size_t k = 0; k < coord.r.size(); ++k) {
            IC_rad[j][k] = IC_radiation(coord.r[k], Gamma[j][k], B[j][k], electrons[j][k], photons[j][k], medium);
        }
    }
    return IC_rad;
}

inline double eta_rad(double nu_m, double nu_c, double pel) { return nu_c < nu_m ? 1 : pow(nu_c / nu_m, 2 - pel); }

double IC_Y_tilt(double b, size_t order) {
    switch (order) {
        case 0:
            return 0;
        case 1:
            return (sqrt(1 + 4 * b) - 1) / 2;
        case 2: {
            double b2 = b * b;
            double x1 = pow(5.19615242270663 * sqrt(27 * b2 + 14 * b + 3) + 27 * b + 7, 1.0 / 3);
            double Y = 0.26456684199469993 * x1 - 0.8399473665965821 / x1 - 1.0 / 3;
            return Y + Y * Y;
        }
        default:
            return 0;
    }
}

double KN_factor(double gamma_e, double nu_p) {
    double x = gamma_e * con::me * con::c2 / (con::h * nu_p);
    return std::min(1.0, x * x);
}

auto calc_Y_nu_c(size_t order, double Gamma, double B, double t_com, double eps_e, double eps_B, SynRad const& electron,
                 ICRad const& photon) {
    double nu_c = electron.nu_c;
    double nu_m = electron.nu_m;
    double pel = electron.pel;
    double eta_e = eta_rad(nu_m, nu_c, pel);
    double b = eta_e * eps_e / eps_B;
    double Y0 = IC_Y_tilt(b, order);
    double Y1 = 2 * Y0;
    for (; fabs((Y1 - Y0) / Y0) > 1e-6;) {
        Y1 = Y0;
        double gamma_c = calc_syn_gamma_c(Gamma, t_com, B, Y1);
        nu_c = syn_nu(gamma_c, B);
        eta_e = eta_rad(nu_m, nu_c, pel);
        b = eta_e * eps_e / eps_B;
        Y0 = IC_Y_tilt(b, order);
    }
    return std::make_pair(Y0, nu_c);
}

auto calc_Y_nu_c_KN(size_t order, double Gamma, double B, double t_com, double eps_e, double eps_B,
                    SynRad const& electron, ICRad const& photon) {
    double nu_c = electron.nu_c;
    double nu_m = electron.nu_m;
    double nu_M = electron.nu_M;
    double gamma_M = syn_gamma(nu_M, B);
    double pel = electron.pel;
    double eta_e = eta_rad(nu_m, nu_c, pel);
    double b = eta_e * eps_e / eps_B;
    double Y0 = IC_Y_tilt(b, order);

    double Y1 = 2 * Y0;
    for (; fabs((Y1 - Y0) / Y0) > 1e-6;) {
        Y1 = Y0;
        double gamma_c = calc_syn_gamma_c(Gamma, t_com, B, Y1);
        nu_c = syn_nu(gamma_c, B);
        eta_e = eta_rad(nu_m, nu_c, pel);
        b = eta_e * eps_e / eps_B * KN_factor(gamma_c, photon.nu_E_peak);
        Y0 = IC_Y_tilt(b, order);
    }
    return std::make_pair(Y0, nu_c);
}

MeshGrid IC_cooling_noKN(MeshGrid const& Gamma, MeshGrid const& t_com, MeshGrid const& B, SynRadMesh& electrons,
                         ICRadMesh const& photons, Medium const& medium, size_t order) {
    MeshGrid Y_tilt = createGrid_like(Gamma);

    for (size_t j = 0; j < electrons.size(); ++j) {
        for (size_t k = 0; k < electrons[j].size(); ++k) {
            auto [Y_, nu_c] = calc_Y_nu_c(order, Gamma[j][k], B[j][k], t_com[j][k], medium.eps_e, medium.eps_B,
                                          electrons[j][k], photons[j][k]);
            electrons[j][k].nu_c = nu_c;
            electrons[j][k].nu_M = syn_nu(calc_syn_gamma_M(B[j][k], medium.zeta, Y_), B[j][k]);
            electrons[j][k].nu_E_peak = syn_nu_E_peak(electrons[j][k]);
            Y_tilt[j][k] = Y_;
        }
    }
    return Y_tilt;
}

MeshGrid IC_cooling(MeshGrid const& Gamma, MeshGrid const& t_com, MeshGrid const& B, SynRadMesh& electrons,
                    ICRadMesh const& photons, Medium const& medium, size_t order) {
    MeshGrid Y_tilt = createGrid_like(Gamma);

    for (size_t j = 0; j < electrons.size(); ++j) {
        for (size_t k = 0; k < electrons[j].size(); ++k) {
            auto [Y_, nu_c] = calc_Y_nu_c_KN(order, Gamma[j][k], B[j][k], t_com[j][k], medium.eps_e, medium.eps_B,
                                             electrons[j][k], photons[j][k]);
            electrons[j][k].nu_c = nu_c;
            electrons[j][k].nu_M = syn_nu(calc_syn_gamma_M(B[j][k], medium.zeta, Y_), B[j][k]);
            electrons[j][k].nu_E_peak = syn_nu_E_peak(electrons[j][k]);
            Y_tilt[j][k] = Y_;
        }
    }
    return Y_tilt;
}

// nu : photon frequency at electron rest frame
