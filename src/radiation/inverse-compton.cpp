//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "inverse-compton.h"

#include <cmath>
#include <iostream>
#include <thread>

#include "macros.h"
#include "utilities.h"

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Computes the radiative efficiency parameter (ηₑ) given minimum electron Lorentz factors.
 * @details If gamma_c is less than gamma_m, it returns 1; otherwise, it returns (gamma_c/gamma_m)^(2-p).
 * @param gamma_m Minimum electron Lorentz factor
 * @param gamma_c Cooling electron Lorentz factor
 * @param p Electron distribution power-law index
 * @return The radiative efficiency parameter
 * <!-- ************************************************************************************** -->
 */
inline Real eta_rad(Real gamma_m, Real gamma_c, Real p) {
    return gamma_c < gamma_m ? 1 : std::pow(gamma_c / gamma_m, (2 - p));
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Computes the effective Compton Y parameter in the Thomson regime.
 * @details Iteratively solves for Y until convergence using the relation:
 *          Y0 = (sqrt(1+4b) - 1)/2, where b = (ηₑ * eps_e / eps_B).
 *          The electron cooling parameters are updated during each iteration.
 * @param B Magnetic field strength
 * @param t_com Comoving time
 * @param eps_e Electron energy fraction
 * @param eps_B Magnetic energy fraction
 * @param e Synchrotron electron properties
 * @return The effective Thomson Y parameter
 * <!-- ************************************************************************************** -->
 */
Real effectiveYThomson(Real B, Real t_com, Real eps_e, Real eps_B, SynElectrons const& e) {
    Real eta_e = eta_rad(e.gamma_m, e.gamma_c, e.p);
    Real b = eta_e * eps_e / eps_B;
    Real Y0 = (std::sqrt(1 + 4 * b) - 1) / 2;
    Real Y1 = 2 * Y0;
    for (; std::fabs((Y1 - Y0) / Y0) > 1e-5;) {
        Y1 = Y0;
        Real gamma_c = compute_gamma_c(t_com, B, e.Ys, e.p);
        eta_e = eta_rad(e.gamma_m, gamma_c, e.p);
        b = eta_e * eps_e / eps_B;
        Y0 = (std::sqrt(1 + 4 * b) - 1) / 2;
    }
    return Y0;
}

Real ICPhoton::I_nu(Real nu) const { return eq_space_loglog_interp(nu, this->nu_IC_, this->j_nu_, true, true); }

Real compton_sigma(Real nu) {
    Real x = con::h * nu / (con::me * con::c2);
    if (x <= 1) {
        return con::sigmaT;
    } else {
        return 0;
    }
    /*
    if (x < 1e-2) {
         return con::sigmaT * (1 - 2 * x);
     } else if (x > 1e2) {
         return 3. / 8 * con::sigmaT * (log(2 * x) + 1.0 / 2) / x;
     } else {
         return 0.75 * con::sigmaT *
                ((1 + x) / (x * x * x) * (2 * x * (1 + x) / (1 + 2 * x) - log(1 + 2 * x)) + log(1 + 2 * x) / (2 * x) -
                 (1 + 3 * x) / (1 + 2 * x) / (1 + 2 * x));
     }
    */
}

ICPhotonGrid gen_IC_photons(SynElectronGrid const& e, SynPhotonGrid const& ph) {
    size_t phi_size = e.shape()[0];
    size_t theta_size = e.shape()[1];
    size_t r_size = e.shape()[2];
    ICPhotonGrid IC_ph({phi_size, theta_size, r_size});

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < r_size; ++k) {
                // Generate the IC photon spectrum for each grid cell.
                IC_ph(i, j, k).gen(e(i, j, k), ph(i, j, k));
            }
        }
    }
    return IC_ph;
}

void Thomson_cooling(SynElectronGrid& e, SynPhotonGrid const& ph, Shock const& shock) {
    size_t phi_size = e.shape()[0];
    size_t theta_size = e.shape()[1];
    size_t r_size = e.shape()[2];

    for (size_t i = 0; i < phi_size; i++) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < r_size; ++k) {
                Real Y_T =
                    effectiveYThomson(shock.B(i, j, k), shock.t_comv(i, j, k), shock.eps_e, shock.eps_B, e(i, j, k));

                e(i, j, k).Ys = InverseComptonY(Y_T);
            }
        }
    }
    update_electrons_4Y(e, shock);
}

void KN_cooling(SynElectronGrid& e, SynPhotonGrid const& ph, Shock const& shock) {
    size_t phi_size = e.shape()[0];
    size_t theta_size = e.shape()[1];
    size_t r_size = e.shape()[2];
    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < r_size; ++k) {
                Real Y_T =
                    effectiveYThomson(shock.B(i, j, k), shock.t_comv(i, j, k), shock.eps_e, shock.eps_B, e(i, j, k));
                // Clear existing Ys and emplace a new InverseComptonY with additional synchrotron frequency parameters.
                // e[i][j][k].Ys.clear();
                // e[i][j][k].Ys.emplace_back(ph[i][j][k].nu_m, ph[i][j][k].nu_c, shock.B[i][j][k], Y_T);
                e(i, j, k).Ys = InverseComptonY(ph(i, j, k).nu_m, ph(i, j, k).nu_c, shock.B(i, j, k), Y_T);
            }
        }
    }
    update_electrons_4Y(e, shock);
}