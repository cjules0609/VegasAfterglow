#include "inverse-compton.h"

#include <cmath>
#include <iostream>
#include <thread>

#include "macros.h"
#include "utilities.h"
ICPhotonMesh create_IC_photon_grid(size_t theta_size, size_t r_size) {
    return ICPhotonMesh(theta_size, ICPhotonArray(r_size));
}

inline bool order(double a, double b, double c) { return a < b && b < c; };

double ICPhoton::L_nu(double nu) const { return loglog_interp_eq_spaced(nu, this->nu_IC_, this->j_nu_, true, true); }

double ICPhoton::E_nu(double nu) const { return L_nu(nu) * dt_com_; }

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
        double gamma_c = syn_gamma_c(t_com, B, e.Ys, e.p);
        eta_e = eta_rad(e.gamma_m, gamma_c, e.p);
        b = eta_e * eps_e / eps_B;
        Y0 = IC_Y_tilt(b);
    }
    return Y0;
}

/*
double eff_Y_IC_KN(double Gamma, double B, double t_com, double eps_e, double eps_B, SynElectrons& e) {
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
}*/

MeshGrid solve_SSC_Y_Thomson(SynElectronsMesh const& e, Shock const& shock) {
    MeshGrid Y_eff = create_grid_like(shock.Gamma);

    for (size_t j = 0; j < e.size(); ++j) {
        for (size_t k = 0; k < e[j].size(); ++k) {
            Y_eff[j][k] = eff_Y_IC_Thomson(shock.Gamma[j][k], shock.B[j][k], shock.t_com[j][k], shock.eps_e,
                                           shock.eps_B, e[j][k]);
        }
    }
    return Y_eff;
}

/*
MeshGrid solve_IC_Y_KN(SynElectronsMesh& e, Shock const& shock) {
    MeshGrid Y_eff = create_grid_like(shock.Gamma);

    for (size_t j = 0; j < e.size(); ++j) {
        for (size_t k = 0; k < e[j].size(); ++k) {
            Y_eff[j][k] =
                eff_Y_IC_KN(shock.Gamma[j][k], shock.B[j][k], shock.t_com[j][k], shock.eps_e, shock.eps_B, e[j][k]);
        }
    }
    return Y_eff;
}*/

void gen_IC_photons_(ICPhotonArray& IC_ph, SynElectronsArray const& e, SynPhotonsArray const& ph, Array const& D_com) {
    for (size_t k = 0; k < IC_ph.size(); ++k) {
        IC_ph[k].gen(e[k], ph[k], D_com[k]);
    }
}

ICPhotonMesh gen_IC_photons(SynElectronsMesh const& e, SynPhotonsMesh const& ph, Shock const& shock) {
    ICPhotonMesh IC_ph = create_IC_photon_grid(shock.width_eff.size(), shock.width_eff[0].size());

    // multithreading acceleration
    std::vector<std::thread> thread;
    thread.reserve(IC_ph.size());
    for (size_t j = 0; j < IC_ph.size(); ++j) {
        thread.emplace_back(gen_IC_photons_, std::ref(IC_ph[j]), std::cref(e[j]), std::cref(ph[j]),
                            std::cref(shock.width_eff[j]));
    }

    for (size_t j = 0; j < thread.size(); ++j) {
        thread[j].join();
    }

    /*for (size_t j = 0; j < IC_ph.size(); ++j) {
        for (size_t k = 0; k < IC_ph[0].size(); ++k) {
            IC_ph[j][k].gen(e[j][k], ph[j][k], shock.D_com[j][k]);
        }
    }*/
    return IC_ph;
}

void IC_cooling_Thomson(SynElectronsMesh& e, SynPhotonsMesh const& ph, Shock const& shock) {
    for (size_t j = 0; j < e.size(); ++j) {
        for (size_t k = 0; k < e[j].size(); ++k) {
            double Y_T = eff_Y_IC_Thomson(shock.Gamma[j][k], shock.B[j][k], shock.t_com[j][k], shock.eps_e, shock.eps_B,
                                          e[j][k]);

            e[j][k].Ys.clear();
            e[j][k].Ys.emplace_back(Y_T);
        }
    }
    update_electrons_4_Y(e, shock);
}

void IC_cooling_KN(SynElectronsMesh& e, SynPhotonsMesh const& ph, Shock const& shock) {
    for (size_t j = 0; j < e.size(); ++j) {
        for (size_t k = 0; k < e[j].size(); ++k) {
            double Y_T = eff_Y_IC_Thomson(shock.Gamma[j][k], shock.B[j][k], shock.t_com[j][k], shock.eps_e, shock.eps_B,
                                          e[j][k]);
            e[j][k].Ys.clear();
            e[j][k].Ys.emplace_back(ph[j][k].nu_m, ph[j][k].nu_c, shock.B[j][k], Y_T);
        }
    }
    update_electrons_4_Y(e, shock);
}