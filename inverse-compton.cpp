#include "inverse-compton.h"

#include <cmath>
#include <iostream>

#include "macros.h"
#include "utilities.h"
ICPhotonMesh create_IC_photon_grid(size_t theta_size, size_t r_size) {
    return ICPhotonMesh(theta_size, ICPhotonArray(r_size));
}

inline bool order(double a, double b, double c) { return a < b && b < c; };

double ICPhoton::j_nu(double nu) const {
    if (nu >= nu_max) {
        return 0.0;
    } else {
        return interp_log_extra_lo(nu, nu_IC_, j_nu_);
    }
}

/*
double ICPhoton::I_nu(double nu) const { return I_nu_peak * I_nu_(nu); }

double ICPhoton::I_nu_(double nu) const {
    if (order(nu_ma, nu_mm, nu_mc)) {
        if (nu < nu_ma) {
            return 5. / 2 * (p - 1) / (p + 1) * pow(nu_ma / nu_mm, 1. / 3) * (nu / nu_ma);
        } else if (nu < nu_mm) {
            return 3. / 2 * (p - 1) / (p - 1. / 3) * pow(nu / nu_mm, 1. / 3);
        } else if (nu < nu_mc) {
            return (p - 1) / (p + 1) * pow(nu / nu_mm, (1 - p) / 2) *
                   (4 * (p + 1. / 3) / (p + 1) / (p - 1. / 3) + log(nu / nu_mm));
        } else if (nu < nu_cc) {
            return (p - 1) / (p + 1) * pow(nu / nu_mm, (1 - p) / 2) *
                   (2 * (2 * p + 3) / (p + 2) - 2 / (p + 1) / (p + 2) + log(nu_cc / nu));
        } else {
            return (p - 1) / (p + 1) * pow(nu / nu_mm, -p / 2) * (nu_mc / nu_mm) *
                   (2 * (2 * p + 3) / (p + 2) - 2 / (p + 2) / (p + 2) + (p + 1) / (p + 2) * log(nu / nu_cc));
        }
    } else if (order(nu_mm, nu_ma, nu_mc)) {
        if (nu < nu_ma) {
            return 2 * (p + 4) * (p - 1) / 3 / (p + 1) / (p + 1) * pow(nu_mm / nu_ma, (p + 1) / 2) * (nu / nu_mm);
        } else if (nu < nu_mc) {
            return (p - 1) / (p + 1) * pow(nu / nu_mm, (1 - p) / 2) *
                   (2 * (2 * p + 5) / (p + 1) / (p + 4) + log(nu / nu_ma));
        } else if (nu < nu_ca) {
            return (p - 1) / (p + 1) * pow(nu / nu_mm, (1 - p) / 2) * (2 + 2 / (p + 4) + log(nu / nu_mc));
        } else if (nu < nu_cc) {
            return (p - 1) / (p + 1) * pow(nu / nu_mm, (1 - p) / 2) * (2 * (2 * p + 1) / (p + 1) + log(nu_cc / nu));
        } else {
            return (p - 1) / (p + 1) * (nu_mc / nu_mm) * pow(nu / nu_mm, -p / 2) *
                   (2 * (2 * p + 5) / (p + 2) + log(nu / nu_cc));
        }
    } else if (order(nu_ma, nu_mc, nu_mm)) {
        if (nu < nu_ca) {
            return 5. / 6 * pow(nu_ma / nu_mc, 1. / 3) * (nu / nu_ca);
        } else if (nu < nu_cc) {
            return 9. / 10 * pow(nu / nu_cc, 1. / 3);
        } else if (nu < nu_cm) {
            return 1. / 3 * pow(nu / nu_cc, -1. / 2) * (28. / 15 + log(nu / nu_ca));
        } else if (nu < nu_mm) {
            return 1. / 3 * pow(nu / nu_cc, -1. / 2) *
                   (2 * (p + 5) / (p + 2) / (p - 1) - 2. / 3 * (p - 1) / (p + 2) + log(nu_mm / nu));
        } else {
            return 1. / (p + 2) * (nu_mc / nu_mm) * pow(nu / nu_mm, -p / 2) *
                   (2. / 3 * (p + 5) / (p - 1) - 2. / 3 * (p - 1) / (p + 2) + log(nu / nu_mm));
        }
    } else if (order(nu_mc, nu_ma, nu_mm)) {
        double R = pow(nu_mc / nu_ma, 1. / 2);
        if (nu < nu_aa) {
            return (R / 2 + 1) * (R + 4) * (nu / nu_aa);
        } else if (nu < nu_am) {
            return R * pow(nu / nu_aa, -1. / 2) * (R / 6 + 0.9 + R / 4 * log(nu / nu_aa));
        } else if (nu < nu_mm) {
            return R * R * pow(nu / nu_aa, -1. / 2) * (3 / (p - 1) - 1. / 2 + 0.75 * log(nu_mm / nu));
        } else {
            return 9 * R * R / 2 / (p + 2) * (nu_ma / nu_mm) * pow(nu / nu_mm, -p / 2) *
                   (4 / (p + 3) * pow(nu_am / nu_mm, (p - 1) / 2) * sqrt(nu_am / nu_cm) +
                    3 * (p + 1) / (p - 1) / (p + 2) + 0.5 * log(nu / nu_mm));
        }
    } else if (order(nu_mm, nu_mc, nu_ma)) {
        double R = (p - 1) / 3 * pow(nu_mc / nu_ma, 1. / 2) * pow(nu_mm / nu_ma, (p - 1) / 2);
        if (nu < nu_aa) {
            return (3 * R / 2 / (p + 2) + 1) * (3 * R / (p + 2) + 4) * (nu / nu_aa);
        } else {
            return pow(nu / nu_aa, -p / 2) / (p + 2) *
                   (6 * R / (p + 3) + R * (9 * R / 2 / (p + 2) + 1) + 9 * R * R / 4 * log(nu / nu_aa));
        }
    } else if (order(nu_mc, nu_mm, nu_ma)) {
        double R = 1. / 3 * pow(nu_mc / nu_ma, 1. / 2) * pow(nu_mm / nu_ma, (p - 1) / 2);
        if (nu < nu_aa) {
            return (3 * R / 2 / (p + 2) + 1) * (3 * R / (p + 2) + 4) * (nu / nu_aa);
        } else {
            return pow(nu / nu_aa, -p / 2) / (p + 2) *
                   (6 * R / (p + 3) + R * (9 * R / 2 / (p + 2) + 1) + 9 * R * R / 4 * log(nu / nu_aa));
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
*/
inline double eta_rad(double gamma_m, double gamma_c, double p) {
    return gamma_c < gamma_m ? 1 : pow(gamma_c / gamma_m, (2 - p));
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
        double gamma_c = syn_gamma_c(t_com, B, Y1);
        eta_e = eta_rad(e.gamma_m, gamma_c, e.p);
        b = eta_e * eps_e / eps_B;
        Y0 = IC_Y_tilt(b);
    }
    return Y0;
}

double eff_Y_IC_KN(double Gamma, double B, double t_com, double eps_e, double eps_B, SynElectrons const& e) {
    double eta_e = eta_rad(e.gamma_m, e.gamma_c, e.p);
    double b = eta_e * eps_e / eps_B;
    double Y0 = IC_Y_tilt(b);
    double Y1 = 2 * Y0;
    for (; fabs((Y1 - Y0) / Y0) > 1e-6;) {
        Y1 = Y0;
        double gamma_c = syn_gamma_c(t_com, B, Y1);
        eta_e = eta_rad(e.gamma_m, gamma_c, e.p);
        double nu_c = syn_nu(gamma_c, B);
        double gamma_N_peak = syn_gamma_N_peak(e.gamma_a, e.gamma_m, gamma_c);
        b = eta_e * eps_e / eps_B * compton_sigma(nu_c / gamma_N_peak) / con::sigmaT;
        Y0 = IC_Y_tilt(b);
    }
    return Y0;
}

MeshGrid solve_IC_Y_Thomson(Shock const& shock, SynElectronsMesh const& e, Medium const& medium) {
    MeshGrid Y_eff = create_grid_like(shock.Gamma);

    for (size_t j = 0; j < e.size(); ++j) {
        for (size_t k = 0; k < e[j].size(); ++k) {
            Y_eff[j][k] = eff_Y_IC_Thomson(shock.Gamma[j][k], shock.B[j][k], shock.t_com[j][k], medium.eps_e,
                                           medium.eps_B, e[j][k]);
        }
    }
    return Y_eff;
}

MeshGrid solve_IC_Y_KN(Shock const& shock, SynElectronsMesh const& e, Medium const& medium) {
    MeshGrid Y_eff = create_grid_like(shock.Gamma);

    for (size_t j = 0; j < e.size(); ++j) {
        for (size_t k = 0; k < e[j].size(); ++k) {
            Y_eff[j][k] =
                eff_Y_IC_KN(shock.Gamma[j][k], shock.B[j][k], shock.t_com[j][k], medium.eps_e, medium.eps_B, e[j][k]);
        }
    }
    return Y_eff;
}

ICPhotonMesh gen_IC_photons(Coord const& coord, Shock const& shock, SynElectronsMesh const& e, SynPhotonsMesh const& ph,
                            Medium const& medium) {
    solve_IC_Y_Thomson(shock, e, medium);
    ICPhotonMesh IC_ph = create_IC_photon_grid(coord.theta.size(), coord.r.size());
    for (size_t j = 0; j < coord.theta.size(); ++j) {
        for (size_t k = 0; k < coord.r.size(); ++k) {
            IC_ph[j][k] = ICPhoton(e[j][k], ph[j][k]);
        }
    }
    return IC_ph;
}
/*
ICPhotonMesh gen_IC_photons(Coord const& coord, Shock const& shock, SynElectronsMesh const& e, SynPhotonsMesh const& ph,
                            Medium const& medium) {
    solve_IC_Y_Thomson(shock, e, medium);

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
            IC_ph[j][k].p = e[j][k].p;
        }
    }
    return IC_ph;
}
*/