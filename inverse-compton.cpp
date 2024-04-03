#include "inverse-compton.h"

#include <cmath>
#include <iostream>

#include "macros.h"
#include "utilities.h"
ICPhotonMesh create_IC_photon_grid(size_t theta_size, size_t r_size, ICPhoton val) {
    return ICPhotonMesh(theta_size, ICPhotonArray(r_size, val));
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

double calc_IC_I_nu_peak(double I_nu_peak, double r, double Gamma, double rho, double xi, double x0) {
    double tau_es = con::sigmaT * (rho / con::mp) * xi * co_moving_shock_width(r, Gamma) * Gamma;

    return tau_es * I_nu_peak * x0;
}

double IC_nu(double gamma, double nu, double x0) { return 4 * gamma * gamma * nu * x0; }

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

ICPhoton gen_IC_photon(double r, double Gamma, double B, SynElectron const& electron, SynElectron const& photon,
                       Medium const& medium) {
    double x0 = IC_x0;
    double rho = medium.rho(r);
    double gamma_m = syn_gamma(electron.nu_m, B);
    double gamma_c = syn_gamma(electron.nu_c, B);
    double gamma_a = syn_gamma(electron.nu_a, B);

    ICPhoton IC_ph;

    IC_ph.I_nu_peak = calc_IC_I_nu_peak(photon.I_nu_peak, r, Gamma, rho, medium.xi, x0);

    IC_ph.nu_mm = IC_nu(gamma_m, photon.nu_m, x0);
    IC_ph.nu_mc = IC_nu(gamma_m, photon.nu_c, x0);
    IC_ph.nu_ma = IC_nu(gamma_m, photon.nu_a, x0);

    IC_ph.nu_cm = IC_nu(gamma_c, photon.nu_m, x0);
    IC_ph.nu_cc = IC_nu(gamma_c, photon.nu_c, x0);
    IC_ph.nu_ca = IC_nu(gamma_c, photon.nu_a, x0);

    IC_ph.nu_am = IC_nu(gamma_a, photon.nu_m, x0);
    IC_ph.nu_ac = IC_nu(gamma_a, photon.nu_c, x0);
    IC_ph.nu_aa = IC_nu(gamma_a, photon.nu_a, x0);

    IC_ph.nu_E_peak = IC_nu_E_peak(IC_ph);

    IC_ph.pel = electron.pel;

    return IC_ph;
}

inline double eta_rad(double nu_m, double nu_c, double pel) { return nu_c < nu_m ? 1 : pow(nu_c / nu_m, 2 - pel); }

double IC_Y_tilt(double b) { return (sqrt(1 + 4 * b) - 1) / 2; }

double KN_correct(double gamma_e, double nu_p) {
    double x = gamma_e * con::me * con::c2 / (con::h * nu_p);
    return std::min(1.0, x * x);
}

double eff_Y_IC_Thomson(double Gamma, double B, double t_com, double eps_e, double eps_B, SynElectron const& electron) {
    double nu_c = electron.nu_c;
    double nu_m = electron.nu_m;
    double pel = electron.pel;
    double eta_e = eta_rad(nu_m, nu_c, pel);
    double b = eta_e * eps_e / eps_B;
    double Y0 = IC_Y_tilt(b);
    double Y1 = 2 * Y0;
    for (; fabs((Y1 - Y0) / Y0) > 1e-6;) {
        Y1 = Y0;
        double gamma_c = syn_gamma_c(Gamma, t_com, B, Y1);
        nu_c = syn_nu(gamma_c, B);
        eta_e = eta_rad(nu_m, nu_c, pel);
        b = eta_e * eps_e / eps_B;
        Y0 = IC_Y_tilt(b);
    }
    return Y0;
}

double eff_Y_IC_KN(double Gamma, double B, double t_com, double eps_e, double eps_B, SynElectron const& electron,
                   ICPhoton const& photon) {
    double nu_c = electron.nu_c;
    double nu_m = electron.nu_m;
    double nu_M = electron.nu_M;
    double gamma_M = syn_gamma(nu_M, B);
    double pel = electron.pel;
    double eta_e = eta_rad(nu_m, nu_c, pel);
    double b = eta_e * eps_e / eps_B;
    double Y0 = IC_Y_tilt(b);

    double Y1 = 2 * Y0;
    for (; fabs((Y1 - Y0) / Y0) > 1e-6;) {
        Y1 = Y0;
        double gamma_c = syn_gamma_c(Gamma, t_com, B, Y1);
        nu_c = syn_nu(gamma_c, B);
        eta_e = eta_rad(nu_m, nu_c, pel);
        b = eta_e * eps_e / eps_B * KN_correct(gamma_c, photon.nu_E_peak);
        Y0 = IC_Y_tilt(b);
    }
    return Y0;
}

ICPhotonMesh gen_IC_photons(Coord const& coord, Shock const& shock, SynElectronMesh const& electrons,
                            SynElectronMesh const& photons, Medium const& medium) {
    ICPhotonMesh IC_photons = create_IC_photon_grid(coord.theta.size(), coord.r.size());
    for (size_t j = 0; j < coord.theta.size(); ++j) {
        for (size_t k = 0; k < coord.r.size(); ++k) {
            IC_photons[j][k] =
                gen_IC_photon(coord.r[k], shock.Gamma[j][k], shock.B[j][k], electrons[j][k], photons[j][k], medium);
        }
    }
    return IC_photons;
}

MeshGrid IC_cooling_Thomson(Shock const& shock, SynElectronMesh& e, ICPhotonMesh const& ph, Medium const& medium) {
    MeshGrid Y_eff = create_grid_like(shock.Gamma);

    for (size_t j = 0; j < e.size(); ++j) {
        for (size_t k = 0; k < e[j].size(); ++k) {
            Y_eff[j][k] = eff_Y_IC_Thomson(shock.Gamma[j][k], shock.B[j][k], shock.t_com[j][k], medium.eps_e,
                                           medium.eps_B, e[j][k]);
            double gamma_c = syn_gamma_c(shock.Gamma[j][k], shock.t_com[j][k], shock.B[j][k], Y_eff[j][k]);
            double gamma_M = syn_gamma_M(shock.B[j][k], medium.zeta, Y_eff[j][k]);
            e[j][k].nu_c = syn_nu(gamma_c, shock.B[j][k]);
            e[j][k].nu_M = syn_nu(gamma_M, shock.B[j][k]);
            e[j][k].nu_E_peak = syn_nu_E_peak(e[j][k]);
        }
    }
    return Y_eff;
}

MeshGrid IC_cooling_KN(Shock const& shock, SynElectronMesh& e, ICPhotonMesh const& ph, Medium const& medium) {
    MeshGrid Y_eff = create_grid_like(shock.Gamma);

    for (size_t j = 0; j < e.size(); ++j) {
        for (size_t k = 0; k < e[j].size(); ++k) {
            Y_eff[j][k] = eff_Y_IC_KN(shock.Gamma[j][k], shock.B[j][k], shock.t_com[j][k], medium.eps_e, medium.eps_B,
                                      e[j][k], ph[j][k]);
            double gamma_c = syn_gamma_c(shock.Gamma[j][k], shock.t_com[j][k], shock.B[j][k], Y_eff[j][k]);
            double gamma_M = syn_gamma_M(shock.B[j][k], medium.zeta, Y_eff[j][k]);
            e[j][k].nu_c = syn_nu(gamma_c, shock.B[j][k]);
            e[j][k].nu_M = syn_nu(gamma_M, shock.B[j][k]);
            e[j][k].nu_E_peak = syn_nu_E_peak(e[j][k]);
        }
    }
    return Y_eff;
}
