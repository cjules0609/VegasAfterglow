//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "observer.h"

#include <cmath>

#include "macros.h"
#include "physics.h"
#include "utilities.h"

/********************************************************************************************************************
 * METHOD: LogScaleInterp::interpLuminosity
 * DESCRIPTION: Interpolates the luminosity at a given observation time (t) using linear
 *              interpolation in log-space. The result is returned in linear space.
 ********************************************************************************************************************/

Real LogScaleInterp::interpLuminosity(Real t_obs) const {
    Real log_t = fastLog(t_obs / t_obs_lo);
    return L_lo * fastExp(log_t * log_L_ratio / log_t_ratio);
}

/********************************************************************************************************************
 * TEMPLATE METHOD: Observer::changeViewingAngle
 * DESCRIPTION: Updates the Observer's viewing angle
 ********************************************************************************************************************/
void Observer::changeViewingAngle(Real theta_view) {
    // Set effective phi grid size based on the observation angle and jet dimensionality.
    if (theta_view == 0 && interp.jet_3d == 0) {
        eff_phi_size = 1;
    } else {
        eff_phi_size = coord.phi.size();
    }
    theta_obs = theta_view;
    calcSolidAngle();
    calcObsTimeGrid();
}

/********************************************************************************************************************
 * METHOD: Observer::calcSolidAngle
 * DESCRIPTION: Calculates the solid angle (dOmega) for each effective phi and theta grid point.
 *              The solid angle is computed as the product of the differential cosine of theta and either 2π (if
 *              the effective phi size is 1) or the differential phi value.
 ********************************************************************************************************************/
void Observer::calcSolidAngle() {
    for (size_t i = 0; i < coord.phi.size(); ++i) {
        Real phi_lo = 0;
        Real phi_hi = 0;
        if (eff_phi_size == 1) {
            phi_lo = 0;
            phi_hi = 2 * con::pi;
        } else if (i == 0) {  // note this also implys phi.size() > 1
            phi_lo = coord.phi[0];
            phi_hi = coord.phi[1];
        } else if (i == coord.phi.size() - 1) {
            phi_lo = coord.phi[i - 1];
            phi_hi = coord.phi[i];
        } else {
            phi_lo = coord.phi[i - 1];
            phi_hi = coord.phi[i + 1];
        }
        Real dphi = std::abs(phi_hi - phi_lo) / 2;
        size_t i_eff = i * interp.jet_3d;
        for (size_t j = 0; j < coord.theta.size(); ++j) {
            for (size_t k = 0; k < coord.t.size(); ++k) {
                Real theta_lo = 0;
                Real theta_hi = 0;
                if (j == 0) {
                    theta_lo = theta_grid[i_eff][0][k];
                    theta_hi = 0.5 * (theta_grid[i_eff][1][k] + theta_grid[i_eff][0][k]);
                } else if (j == coord.theta.size() - 1) {
                    theta_lo = 0.5 * (theta_grid[i_eff][j - 1][k] + theta_grid[i_eff][j][k]);
                    theta_hi = theta_grid[i_eff][j][k];
                } else {
                    theta_lo = 0.5 * (theta_grid[i_eff][j][k] + theta_grid[i_eff][j - 1][k]);
                    theta_hi = 0.5 * (theta_grid[i_eff][j][k] + theta_grid[i_eff][j + 1][k]);
                }
                Real dcos = std::abs(std::cos(theta_hi) - std::cos(theta_lo));
                dOmega[i][j][k] = dcos * dphi;
            }
        }
    }

    /*for (size_t i = 0; i < eff_phi_size; ++i) {
        for (size_t j = 0; j < coord.theta.size(); ++j) {
            // If there is only one effective phi point, use 2π; otherwise, use the corresponding differential phi.
            dOmega[i][j] = coord.dcos[j] * (eff_phi_size == 1 ? 2 * con::pi : coord.dphi[i]);
        }
    }*/
}

/********************************************************************************************************************
 * METHOD: Observer::calcObsTimeGrid
 * DESCRIPTION: Calculates the observation time grid (t_obs_grid) and updates the doppler factor grid based on the
 *              Gamma (Lorentz factor) and engine time (t) array.
 *              For each grid point, the Doppler factor is computed and the observed time is calculated taking
 *              redshift into account.
 ********************************************************************************************************************/
void Observer::calcObsTimeGrid() {
    auto [phi_size, theta_size, t_size] = coord.shape();
    Real cos_obs = std::cos(theta_obs);
    Real sin_obs = std::sin(theta_obs);
    for (size_t i = 0; i < eff_phi_size; ++i) {
        Real cos_phi = std::cos(coord.phi[i]);
        for (size_t j = 0; j < theta_size; ++j) {
            // Compute the cosine of the angle between the local velocity vector and the observer's line of sight.
            Real cos_v = std::sin(coord.theta[j]) * cos_phi * sin_obs + std::cos(coord.theta[j]) * cos_obs;
            for (size_t k = 0; k < t_size; ++k) {
                Real gamma_ = Gamma[i * interp.jet_3d][j][k];  // Get Gamma at the grid point.
                Real r = r_grid[i * interp.jet_3d][j][k];
                Real t_eng_ = coord.t[k];         // Get engine time at the grid point.
                Real beta = gammaTobeta(gamma_);  // Convert Gamma to beta.
                // Compute the Doppler factor: D = 1 / [Gamma * (1 - beta * cos_v)]
                doppler[i][j][k] = 1 / (gamma_ * (1 - beta * cos_v));
                // Compute the observed time: t_obs = [t_eng + (1 - cos_v) * r / c] * (1 + z)
                t_obs_grid[i][j][k] = (t_eng_ + (1 - cos_v) * r / con::c) * (1 + z);
            }
        }
    }
}