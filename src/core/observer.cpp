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
 * METHOD: LogScaleInterp::interpRadius
 * DESCRIPTION: Interpolates the radius at a given logarithmic observation time (log_t) using linear interpolation
 *              in log-space. The result is returned in linear space.
 ********************************************************************************************************************/
Real LogScaleInterp::interpRadius(Real log_t) const {
    // Linearly interpolate between the lower and upper log radius boundaries, then exponentiate.
    return fastExp(log_r_lo + (log_r_hi - log_r_lo) * (log_t - log_t_lo) / (log_t_hi - log_t_lo));
}

/********************************************************************************************************************
 * METHOD: LogScaleInterp::interpIntensity
 * DESCRIPTION: Interpolates the intensity at a given logarithmic observation time (log_t) using linear interpolation
 *              in log-space. The result is returned in linear space.
 ********************************************************************************************************************/
Real LogScaleInterp::interpIntensity(Real log_t) const {
    // Linearly interpolate between the lower and upper log intensity boundaries, then exponentiate.
    return fastExp(log_I_lo + (log_I_hi - log_I_lo) * (log_t - log_t_lo) / (log_t_hi - log_t_lo));
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
    calcSolidAngle();
    calcObsTimeGrid();
}

/********************************************************************************************************************
 * METHOD: LogScaleInterp::interpDoppler
 * DESCRIPTION: Interpolates the Doppler factor at a given logarithmic observation time (log_t) using linear
 *              interpolation in log-space. The result is returned in linear space.
 ********************************************************************************************************************/
Real LogScaleInterp::interpDoppler(Real log_t) const {
    // Linearly interpolate between the lower and upper log Doppler boundaries, then exponentiate.
    return fastExp(log_d_lo + (log_d_hi - log_d_lo) * (log_t - log_t_lo) / (log_t_hi - log_t_lo));
}

/********************************************************************************************************************
 * METHOD: Observer::calcSolidAngle
 * DESCRIPTION: Calculates the solid angle (dOmega) for each effective phi and theta grid point.
 *              The solid angle is computed as the product of the differential cosine of theta and either 2π (if
 *              the effective phi size is 1) or the differential phi value.
 ********************************************************************************************************************/
void Observer::calcSolidAngle() {
    for (size_t i = 0; i < eff_phi_size; ++i) {
        for (size_t j = 0; j < coord.theta.size(); ++j) {
            // If there is only one effective phi point, use 2π; otherwise, use the corresponding differential phi.
            dOmega[i][j] = coord.dcos[j] * (eff_phi_size == 1 ? 2 * con::pi : coord.dphi[i]);
        }
    }
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