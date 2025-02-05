//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "observer.h"

#include <boost/numeric/odeint.hpp>
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
 * CONSTRUCTOR: Observer::Observer
 * DESCRIPTION: Constructs an Observer object with the given coordinate grid.
 *              It initializes the doppler and t_obs_grid as 3D grids, dOmega as a 2D grid, and log_r as a 1D array.
 *              The logarithm of each radius value is computed and stored in log_r. Finally, the solid angle is
 *calculated.
 ********************************************************************************************************************/
Observer::Observer(Coord const& coord)
    : t_obs_grid(boost::extents[coord.phi.size()][coord.theta.size()][coord.r.size()]),
      doppler(boost::extents[coord.phi.size()][coord.theta.size()][coord.r.size()]),
      theta_obs(0),
      lumi_dist(1),
      z(0),
      dOmega(boost::extents[coord.phi.size()][coord.theta.size()]),
      log_r(boost::extents[coord.r.size()]),
      interp(),
      coord(coord),
      eff_phi_size(1) {
    // Compute logarithm (using fastLog) for each radius value in the coordinate grid.
    for (size_t i = 0; i < coord.r.size(); ++i) {
        log_r[i] = fastLog(coord.r[i]);
    }

    // Calculate the solid angle grid.
    calcSolidAngle();
}

/********************************************************************************************************************
 * METHOD: Observer::calcObsTimeGrid
 * DESCRIPTION: Calculates the observation time grid (t_obs_grid) and updates the doppler factor grid based on the
 *              provided Gamma (Lorentz factor) and engine time (t_eng) grids.
 *              For each grid point, the Doppler factor is computed and the observed time is calculated taking
 *              redshift into account.
 ********************************************************************************************************************/
void Observer::calcObsTimeGrid(MeshGrid3d const& Gamma, MeshGrid3d const& t_eng) {
    auto [phi_size, theta_size, r_size] = coord.shape();
    Real cos_obs = std::cos(theta_obs);
    Real sin_obs = std::sin(theta_obs);
    for (size_t i = 0; i < eff_phi_size; ++i) {
        Real cos_phi = std::cos(coord.phi[i]);
        for (size_t j = 0; j < theta_size; ++j) {
            // Compute the cosine of the angle between the local velocity vector and the observer's line of sight.
            Real cos_v = std::sin(coord.theta[j]) * cos_phi * sin_obs + std::cos(coord.theta[j]) * cos_obs;
            for (size_t k = 0; k < r_size; ++k) {
                Real gamma_ = Gamma[i * interp.jet_3d][j][k];  // Get Gamma at the grid point.
                Real t_eng_ = t_eng[i * interp.jet_3d][j][k];  // Get engine time at the grid point.
                Real beta = gammaTobeta(gamma_);               // Convert Gamma to beta.
                // Compute the Doppler factor: D = 1 / [Gamma * (1 - beta * cos_v)]
                doppler[i][j][k] = 1 / (gamma_ * (1 - beta * cos_v));
                if (gamma_ == 1) {
                    // For non-relativistic case, set observed time to infinity.
                    t_obs_grid[i][j][k] = std::numeric_limits<Real>::infinity();
                } else {
                    // Compute the observed time: t_obs = [t_eng + (1 - cos_v) * r / c] * (1 + z)
                    t_obs_grid[i][j][k] = (t_eng_ + (1 - cos_v) * coord.r[k] / con::c) * (1 + z);
                }
            }
        }
    }
}