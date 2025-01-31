#include "inverse-compton.h"

#include <cmath>
#include <iostream>
#include <thread>

#include "macros.h"
#include "utilities.h"
ICPhotonMesh createICPhotonGrid(size_t theta_size, size_t r_size) {
    return ICPhotonMesh(theta_size, ICPhotonArray(r_size));
}

inline bool order(double a, double b, double c) { return a < b && b < c; };

double ICPhoton::I_nu(double nu) const { return loglogInterpEqSpaced(nu, this->nu_IC_, this->j_nu_, true, true); }

inline double eta_rad(double gamma_m, double gamma_c, double p) {
    return gamma_c < gamma_m ? 1 : std::pow(gamma_c / gamma_m, (2 - p));
}

double effectiveYThomson(double Gamma, double B, double t_com, double eps_e, double eps_B, SynElectrons const& e) {
    double eta_e = eta_rad(e.gamma_m, e.gamma_c, e.p);
    double b = eta_e * eps_e / eps_B;
    double Y0 = (std::sqrt(1 + 4 * b) - 1) / 2;
    double Y1 = 2 * Y0;
    for (; std::fabs((Y1 - Y0) / Y0) > 1e-5;) {
        Y1 = Y0;
        double gamma_c = syn_gamma_c(t_com, B, e.Ys, e.p);
        eta_e = eta_rad(e.gamma_m, gamma_c, e.p);
        b = eta_e * eps_e / eps_B;
        Y0 = (std::sqrt(1 + 4 * b) - 1) / 2;
    }
    return Y0;
}

ICPhotonMesh genICPhotons(SynElectronsMesh const& e, SynPhotonsMesh const& ph, Shock const& shock) {
    ICPhotonMesh IC_ph = createICPhotonGrid(shock.width_eff.size(), shock.width_eff[0].size());

    for (size_t j = 0; j < IC_ph.size(); ++j) {
        for (size_t k = 0; k < IC_ph[0].size(); ++k) {
            IC_ph[j][k].gen(e[j][k], ph[j][k], shock.width_eff[j][k]);
        }
    }
    return IC_ph;
}

void eCoolingThomson(SynElectronsMesh& e, SynPhotonsMesh const& ph, Shock const& shock) {
    for (size_t j = 0; j < e.size(); ++j) {
        for (size_t k = 0; k < e[j].size(); ++k) {
            double Y_T = effectiveYThomson(shock.Gamma[j][k], shock.B[j][k], shock.t_com[j][k], shock.eps_e,
                                           shock.eps_B, e[j][k]);

            e[j][k].Ys.clear();
            e[j][k].Ys.emplace_back(Y_T);
        }
    }
    updateElectrons4Y(e, shock);
}

void eCoolingKleinNishina(SynElectronsMesh& e, SynPhotonsMesh const& ph, Shock const& shock) {
    for (size_t j = 0; j < e.size(); ++j) {
        for (size_t k = 0; k < e[j].size(); ++k) {
            double Y_T = effectiveYThomson(shock.Gamma[j][k], shock.B[j][k], shock.t_com[j][k], shock.eps_e,
                                           shock.eps_B, e[j][k]);
            e[j][k].Ys.clear();
            e[j][k].Ys.emplace_back(ph[j][k].nu_m, ph[j][k].nu_c, shock.B[j][k], Y_T);
        }
    }
    updateElectrons4Y(e, shock);
}