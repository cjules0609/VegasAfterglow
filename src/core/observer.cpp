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
    Real log_t = fastLog2(t_obs / t_obs_lo);
    return surface_lo * fastExp2(log_I_lo + log_t * slope);
}

/********************************************************************************************************************
 * METHOD: Observer::calcEmissionSurface
 * DESCRIPTION: Calculates the observe frame effective emission surface for each effective phi and theta grid point.
 *              The emission surface is computed as the product of the differential cosine of theta and either 2π (if
 *              the effective phi size is 1) or the differential phi value.
 ********************************************************************************************************************/
void Observer::calcEmissionSurface(Coord const& coord, Shock const& shock) {
    Array dphi({eff_phi_size}, 0);

    for (size_t i = 0; i < eff_phi_size; ++i) {
        Real phi_lo = 0;
        Real phi_hi = 0;
        if (eff_phi_size == 1) {
            phi_lo = 0;
            phi_hi = 4 * con::pi;
        } else if (i == 0) {  // note this also implys phi.size() > 1
            phi_lo = coord.phi(0);
            phi_hi = coord.phi(1);
        } else if (i == coord.phi.size() - 1) {
            phi_lo = coord.phi(i - 1);
            phi_hi = coord.phi(i);
        } else [[likely]] {
            phi_lo = coord.phi(i - 1);
            phi_hi = coord.phi(i + 1);
        }
        dphi(i) = std::abs(phi_hi - phi_lo) / 2;
    }

    for (size_t i = 0; i < eff_phi_size; ++i) {
        size_t i_eff = i * interp.jet_3d;
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < t_size; ++k) {
                if (shock.required(i_eff, j, k) == 0) {
                    continue;
                }
                Real theta_lo = 0;
                Real theta_hi = 0;
                if (j == 0) {
                    theta_lo = shock.theta(i_eff, 0, k);
                    theta_hi = 0.5 * (shock.theta(i_eff, 1, k) + shock.theta(i_eff, 0, k));
                } else if (j == coord.theta.size() - 1) {
                    theta_lo = 0.5 * (shock.theta(i_eff, j - 1, k) + shock.theta(i_eff, j, k));
                    theta_hi = shock.theta(i_eff, j, k);
                } else [[likely]] {
                    theta_lo = 0.5 * (shock.theta(i_eff, j, k) + shock.theta(i_eff, j - 1, k));
                    theta_hi = 0.5 * (shock.theta(i_eff, j, k) + shock.theta(i_eff, j + 1, k));
                }
                Real dcos = std::abs(std::cos(theta_hi) - std::cos(theta_lo));
                Real r = shock.r(i_eff, j, k);
                Real D = doppler(i, j, k);
                surface(i, j, k) = dcos * dphi(i) * r * r * D * D * D;
            }
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
void Observer::calcObsTimeGrid(Coord const& coord, MeshGrid3d const& Gamma, MeshGrid3d const& r_grid, Real theta_obs) {
    Real cos_obs = std::cos(theta_obs);
    Real sin_obs = std::sin(theta_obs);
    for (size_t i = 0; i < eff_phi_size; ++i) {
        Real cos_phi = std::cos(coord.phi[i]);
        for (size_t j = 0; j < theta_size; ++j) {
            // Compute the cosine of the angle between the local velocity vector and the observer's line of sight.
            Real cos_v = std::sin(coord.theta[j]) * cos_phi * sin_obs + std::cos(coord.theta[j]) * cos_obs;
            for (size_t k = 0; k < t_size; ++k) {
                Real gamma_ = Gamma(i * interp.jet_3d, j, k);  // Get Gamma at the grid point.
                Real r = r_grid(i * interp.jet_3d, j, k);
                Real t_eng_ = coord.t[k];  // Get engine time at the grid point.
                // Compute the Doppler factor: D = 1 / [Gamma * (1 - beta * cos_v)]
                doppler(i, j, k) = 1 / (gamma_ - std::sqrt(gamma_ * gamma_ - 1) * cos_v);
                // Compute the observed time: t_obs = [t_eng + (1 - cos_v) * r / c] * (1 + z)
                t_obs_grid(i, j, k) = (t_eng_ + (1 - cos_v) * r / con::c) * (1 + z);
            }
        }
    }
}

void Observer::updateRequired(MaskGrid& required, Array const& t_obs) {
    size_t t_obs_size = t_obs.size();

    // Loop over effective phi and theta grid points.
    for (size_t i = 0; i < eff_phi_size; i++) {
        size_t i_eff = i * interp.jet_3d;
        for (size_t j = 0; j < theta_size; j++) {
            // Skip observation times that are below the grid's start time
            size_t t_idx = iterate_to(t_obs_grid(i, j, 0), t_obs, 0);

            // find the grid points that are required for the interpolation.
            for (size_t k = 0; k < t_size - 1 && t_idx < t_obs_size; k++) {
                Real const t_lo = t_obs_grid(i, j, k);
                Real const t_hi = t_obs_grid(i, j, k + 1);

                if (t_lo <= t_obs(t_idx) && t_obs(t_idx) < t_hi) {
                    required(i_eff, j, k) = 1;
                    required(i_eff, j, k + 1) = 1;
                }

                t_idx = iterate_to(t_hi, t_obs, t_idx);
            }
        }
    }
}

void Observer::buildObsTimeGrid(Coord const& coord, Shock const& shock, Real theta_view, Real luminosity_dist,
                                Real redshift) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    // Determine if the jet is 3D (more than one phi value)
    interp.jet_3d = static_cast<size_t>((phi_size > 1));
    // Set effective phi grid size based on the observation angle and jet dimensionality.
    if (theta_view == 0 && interp.jet_3d == 0) {
        this->eff_phi_size = 1;  // optimize for on-axis observer
    } else {
        this->eff_phi_size = coord.phi.size();
    }
    this->theta_size = theta_size;
    this->t_size = t_size;
    this->lumi_dist = luminosity_dist;
    this->z = redshift;

    t_obs_grid.resize({eff_phi_size, theta_size, t_size});
    doppler.resize({eff_phi_size, theta_size, t_size});
    surface.resize({eff_phi_size, theta_size, t_size});

    // Calculate the solid angle grid and observation time grid.
    calcObsTimeGrid(coord, shock.Gamma, shock.r, theta_view);
}

/********************************************************************************************************************
 * CONSTRUCTOR: Observer::Observer
 * DESCRIPTION: Constructs an Observer object with the given coordinate grid, shocks, observation angle, luminosity
 *              and redshift. It initializes the observation time and Doppler factor grids, as well as the
 *              interpolation object.
 ********************************************************************************************************************/
void Observer::observe(Coord const& coord, Shock const& shock, Real theta_view, Real luminosity_dist, Real redshift) {
    buildObsTimeGrid(coord, shock, theta_view, luminosity_dist, redshift);
    calcEmissionSurface(coord, shock);
}

void Observer::observeAt(Array const& t_obs, Coord const& coord, Shock& shock, Real theta_view, Real luminosity_dist,
                         Real redshift) {
    buildObsTimeGrid(coord, shock, theta_view, luminosity_dist, redshift);

    xt::noalias(shock.required) = 0;

    updateRequired(shock.required, t_obs);

    calcEmissionSurface(coord, shock);
}