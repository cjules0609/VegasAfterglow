#include "inverse-compton.h"

#include <cmath>
#include <iostream>

#include "macros.h"
#include "utilities.h"
ICPhotonMesh create_IC_photon_grid(size_t theta_size, size_t r_size) {
    return ICPhotonMesh(theta_size, ICPhotonArray(r_size));
}
inline const double IC_x0 = sqrt(2) / 3;

double ICPhoton::I_nu(double nu) const { return I_nu_peak * I_nu_(nu); }

inline bool order(double a, double b, double c) { return a < b && b < c; };

double compton_sigma(double nu) {
    double x = con::h * nu / (con::me * con::c2);
    return 0.75 * con::sigmaT *
           ((1 + x) / (x * x * x) * (2 * x * (1 + x) / (1 + 2 * x) - log(1 + 2 * x)) + log(1 + 2 * x) / (2 * x) -
            (1 + 3 * x) / (1 + 2 * x) / (1 + 2 * x));
}

double ICPhoton::I_nu_(double nu) const {
    if (order(nu_ma, nu_mm, nu_mc)) {
        if (nu < nu_ma) {
            return 5. / 2 * (pel - 1) / (pel + 1) * pow(nu_ma / nu_mm, 1. / 3) * (nu / nu_ma);
        } else if (nu < nu_mm) {
            return 3. / 2 * (pel - 1) / (pel - 1. / 3) * pow(nu / nu_mm, 1. / 3);
        } else if (nu < nu_mc) {
            return (pel - 1) / (pel + 1) * pow(nu / nu_mm, (1 - pel) / 2) *
                   (4 * (pel + 1. / 3) / (pel + 1) / (pel - 1. / 3) + log(nu / nu_mm));
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

double calc_IC_I_nu_peak(double syn_I_nu_peak, double r, double Gamma, double D_com, double rho, double xi) {
    double tau_es = con::sigmaT * (rho / con::mp) * xi * D_com * Gamma;
    return tau_es * syn_I_nu_peak * IC_x0;
}

double IC_nu(double gamma, double nu) { return 4 * gamma * gamma * nu * IC_x0; }

double IC_nu_E_peak(ICPhoton const& ph) {
    if (order(ph.nu_ma, ph.nu_mm, ph.nu_mc)) {
        return ph.nu_cc;
    } else if (order(ph.nu_mm, ph.nu_ma, ph.nu_mc)) {
        return ph.nu_cc;
    } else if (order(ph.nu_ma, ph.nu_mc, ph.nu_mm)) {
        return ph.nu_mm;
    } else if (order(ph.nu_mc, ph.nu_ma, ph.nu_mm)) {
        if (ph.nu_ma * ph.nu_ma < ph.nu_mm * ph.nu_mc) {
            return ph.nu_mm;
        } else {
            return ph.nu_aa;
        }
    } else if (order(ph.nu_mm, ph.nu_mc, ph.nu_ma)) {
        return ph.nu_aa;
    } else if (order(ph.nu_mc, ph.nu_mm, ph.nu_ma)) {
        return ph.nu_aa;
    }
}

inline double eta_rad(double gamma_m, double gamma_c, double pel) {
    return gamma_c < gamma_m ? 1 : pow(gamma_c / gamma_m, (2 - pel));
}

double IC_Y_tilt(double b) { return (sqrt(1 + 4 * b) - 1) / 2; }

double KN_correct(double gamma_e, double nu_p) {
    double x = gamma_e * con::me * con::c2 / (con::h * nu_p);
    return std::min(1.0, x * x);
}

double eff_Y_IC_Thomson(double Gamma, double B, double t_com, double eps_e, double eps_B, SynElectrons const& e) {
    double eta_e = eta_rad(e.gamma_m, e.gamma_c, e.p);
    double b = eta_e * eps_e / eps_B;
    double Y0 = IC_Y_tilt(b);
    double Y1 = 2 * Y0;
    for (; fabs((Y1 - Y0) / Y0) > 1e-6;) {
        Y1 = Y0;
        double gamma_c = syn_gamma_c(Gamma, t_com, B, Y1);
        eta_e = eta_rad(e.gamma_m, gamma_c, e.p);
        b = eta_e * eps_e / eps_B;
        Y0 = IC_Y_tilt(b);
    }
    return Y0;
}

double eff_Y_IC_KN(double Gamma, double B, double t_com, double eps_e, double eps_B, SynElectrons const& e, double nu) {
    double eta_e = eta_rad(e.gamma_m, e.gamma_c, e.p);
    double b = eta_e * eps_e / eps_B;
    double Y0 = IC_Y_tilt(b);
    double Y1 = 2 * Y0;
    for (; fabs((Y1 - Y0) / Y0) > 1e-6;) {
        Y1 = Y0;
        double gamma_c = syn_gamma_c(Gamma, t_com, B, Y1);
        double gamma_N_peak = syn_gamma_N_peak(e.gamma_a, e.gamma_m, gamma_c);
        eta_e = eta_rad(e.gamma_m, gamma_c, e.p);
        b = eta_e * eps_e / eps_B * KN_correct(gamma_N_peak, nu);
        Y0 = IC_Y_tilt(b);
    }
    return Y0;
}

MeshGrid IC_cooling_Thomson(Shock const& shock, SynElectronsMesh& e, Medium const& medium) {
    MeshGrid Y_eff = create_grid_like(shock.Gamma);

    for (size_t j = 0; j < e.size(); ++j) {
        for (size_t k = 0; k < e[j].size(); ++k) {
            Y_eff[j][k] = eff_Y_IC_Thomson(shock.Gamma[j][k], shock.B[j][k], shock.t_com[j][k], medium.eps_e,
                                           medium.eps_B, e[j][k]);
            e[j][k].gamma_c = syn_gamma_c(shock.Gamma[j][k], shock.t_com[j][k], shock.B[j][k], Y_eff[j][k]);
            e[j][k].gamma_M = syn_gamma_M(shock.B[j][k], medium.zeta, Y_eff[j][k]);
        }
    }
    return Y_eff;
}

MeshGrid IC_cooling_KN(Shock const& shock, SynElectronsMesh& e, ICPhotonMesh const& ph, Medium const& medium) {
    MeshGrid Y_eff = create_grid_like(shock.Gamma);

    for (size_t j = 0; j < e.size(); ++j) {
        for (size_t k = 0; k < e[j].size(); ++k) {
            Y_eff[j][k] = eff_Y_IC_KN(shock.Gamma[j][k], shock.B[j][k], shock.t_com[j][k], medium.eps_e, medium.eps_B,
                                      e[j][k], ph[j][k].nu_E_peak);
            e[j][k].gamma_c = syn_gamma_c(shock.Gamma[j][k], shock.t_com[j][k], shock.B[j][k], Y_eff[j][k]);
            e[j][k].gamma_M = syn_gamma_M(shock.B[j][k], medium.zeta, Y_eff[j][k]);
        }
    }
    return Y_eff;
}

ICPhotonMesh gen_IC_photons(Coord const& coord, Shock const& shock, SynElectronsMesh& e, SynPhotonsMesh const& ph,
                            Medium const& medium) {
    IC_cooling_Thomson(shock, e, medium);

    ICPhotonMesh IC_ph = create_IC_photon_grid(coord.theta.size(), coord.r.size());
    for (size_t j = 0; j < coord.theta.size(); ++j) {
        for (size_t k = 0; k < coord.r.size(); ++k) {
            double r = coord.r[k];
            double rho = medium.rho(r);
            double Gamma = shock.Gamma[j][k];
            double B = shock.B[j][k];
            double D_com = shock.D_com[j][k];
            double syn_I_nu_peak = ph[j][k].I_nu_peak;

            IC_ph[j][k].I_nu_peak = calc_IC_I_nu_peak(syn_I_nu_peak, r, Gamma, D_com, rho, medium.xi);
            IC_ph[j][k].nu_mm = IC_nu(e[j][k].gamma_m, ph[j][k].nu_m);
            IC_ph[j][k].nu_mc = IC_nu(e[j][k].gamma_m, ph[j][k].nu_c);
            IC_ph[j][k].nu_ma = IC_nu(e[j][k].gamma_m, ph[j][k].nu_a);

            IC_ph[j][k].nu_cm = IC_nu(e[j][k].gamma_c, ph[j][k].nu_m);
            IC_ph[j][k].nu_cc = IC_nu(e[j][k].gamma_c, ph[j][k].nu_c);
            IC_ph[j][k].nu_ca = IC_nu(e[j][k].gamma_c, ph[j][k].nu_a);

            IC_ph[j][k].nu_am = IC_nu(e[j][k].gamma_a, ph[j][k].nu_m);
            IC_ph[j][k].nu_ac = IC_nu(e[j][k].gamma_a, ph[j][k].nu_c);
            IC_ph[j][k].nu_aa = IC_nu(e[j][k].gamma_a, ph[j][k].nu_a);

            IC_ph[j][k].nu_E_peak = IC_nu_E_peak(IC_ph[j][k]);
            IC_ph[j][k].pel = e[j][k].p;
        }
    }
    return IC_ph;
}
