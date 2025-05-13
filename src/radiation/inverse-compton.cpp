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

#include "IO.h"
#include "macros.h"
#include "utilities.h"

void ICPhoton::fill_integration_grid(IntegratorGrid& grid, bool KN) noexcept {
    // For each (nu0, gamma) pair, compute differential contributions and fill in I0
    for (size_t i = 0; i < grid.num; ++i) {
        Real dnu = grid.nu_edge(i + 1) - grid.nu_edge(i);
        for (size_t j = 0; j < grid.num; ++j) {
            Real dgamma = grid.gamma_edge(j + 1) - grid.gamma_edge(j);
            Real dS = std::fabs(dnu * dgamma);
            Real gamma_nu = grid.nu(i) * grid.gamma(j);
            Real factor = 4 * gamma_nu * gamma_nu;
            grid.I0(i, j) = KN && (grid.nu(i) / grid.gamma(j) * (con::h / con::mec2) > 1)
                                ? 0
                                : grid.column_den(j) * grid.I_nu_syn(i) * dS / factor * con::sigmaT;
            /*if (KN) {
                grid.I0(i, j) = grid.column_den(j) * grid.I_nu_syn(i) * dS / factor *
                                compton_cross_section(grid.nu(i) / grid.gamma(j));
            } else {
                grid.I0(i, j) = grid.column_den(j) * grid.I_nu_syn(i) * dS / factor * con::sigmaT;
            }*/
        }
    }
}

void ICPhoton::integrate_IC_spectrum(IntegratorGrid const& grid) noexcept {
    // Calculate output frequency range
    Real gamma_min = grid.gamma.front();
    Real gamma_max = grid.gamma.back();
    Real nu0_min = grid.nu.front();
    Real nu0_max = grid.nu.back();

    Real nu_min = 4 * gamma_min * gamma_min * nu0_min * IC_x0 * 0.999;  // 0.999 to avoid the nu_IC_>max(max_freq)
    Real nu_max = 4 * gamma_max * gamma_max * nu0_max * IC_x0 * 0.999;

    nu_IC_ = xt::logspace(std::log10(nu_min), std::log10(nu_max), spectrum_resol);
    I_nu_IC_ = xt::zeros<Real>({spectrum_resol});

    // Integrate over the grid to compute the final IC photon spectrum
    for (size_t i = 0; i < grid.num; ++i) {
        for (size_t j = 0; j < grid.num; ++j) {
            Real max_freq = 4 * grid.gamma(j) * grid.gamma(j) * grid.nu(i) * IC_x0;
            for (size_t k = 0; k < nu_IC_.size(); ++k) {
                if (nu_IC_(k) <= max_freq) {
                    I_nu_IC_(k) += grid.I0(i, j) * nu_IC_(k);
                } else {
                    break;
                }
            }
        }
    }
}

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
Real compute_Thomson_Y(Real B, Real t_com, Real eps_e, Real eps_B, SynElectrons const& e) {
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

Real ICPhoton::compute_I_nu(Real nu) const noexcept {
    return eq_space_loglog_interp(nu, this->nu_IC_, this->I_nu_IC_, true, true);
}

Real ICPhoton::compute_log2_I_nu(Real log2_nu) const noexcept {
    Real dnu = log2_nu_IC_(1) - log2_nu_IC_(0);
    size_t idx = static_cast<size_t>((log2_nu - log2_nu_IC_(0)) / dnu + 1);
    if (idx < 1) {
        idx = 1;
    } else if (idx > log2_nu_IC_.size() - 1) {
        idx = log2_nu_IC_.size() - 1;
    }

    Real slope = (log2_I_nu_(idx) - log2_I_nu_(idx - 1)) / dnu;

    return log2_I_nu_(idx - 1) + slope * (log2_nu - log2_nu_IC_(idx - 1));
}

Real compton_cross_section(Real nu) {
    Real x = con::h / (con::me * con::c2) * nu;
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

ICPhotonGrid generate_IC_photons(SynElectronGrid const& e, SynPhotonGrid const& ph, bool KN) noexcept {
    size_t phi_size = e.shape()[0];
    size_t theta_size = e.shape()[1];
    size_t t_size = e.shape()[2];
    ICPhotonGrid IC_ph({phi_size, theta_size, t_size});

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < t_size; ++k) {
                // Generate the IC photon spectrum for each grid cell.
                IC_ph(i, j, k).compute_IC_spectrum(e(i, j, k), ph(i, j, k), KN);
            }
        }
    }
    return IC_ph;
}

void Thomson_cooling(SynElectronGrid& e, SynPhotonGrid& ph, Shock const& shock) {
    size_t phi_size = e.shape()[0];
    size_t theta_size = e.shape()[1];
    size_t t_size = e.shape()[2];

    for (size_t i = 0; i < phi_size; i++) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < t_size; ++k) {
                Real Y_T =
                    compute_Thomson_Y(shock.B(i, j, k), shock.t_comv(i, j, k), shock.eps_e, shock.eps_B, e(i, j, k));
                e(i, j, k).Ys = InverseComptonY(Y_T);
            }
        }
    }
    update_electrons_4Y(e, shock);
    generate_syn_photons(ph, shock, e);
}

void KN_cooling(SynElectronGrid& e, SynPhotonGrid& ph, Shock const& shock) {
    size_t phi_size = e.shape()[0];
    size_t theta_size = e.shape()[1];
    size_t r_size = e.shape()[2];
    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < r_size; ++k) {
                Real Y_T =
                    compute_Thomson_Y(shock.B(i, j, k), shock.t_comv(i, j, k), shock.eps_e, shock.eps_B, e(i, j, k));
                // Clear existing Ys and emplace a new InverseComptonY with additional synchrotron frequency parameters.
                e(i, j, k).Ys = InverseComptonY(ph(i, j, k).nu_m, ph(i, j, k).nu_c, shock.B(i, j, k), Y_T);
            }
        }
    }
    update_electrons_4Y(e, shock);
    generate_syn_photons(ph, shock, e);
}
