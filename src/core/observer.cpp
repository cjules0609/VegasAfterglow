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
 * METHOD: Observer::interpolate
 * DESCRIPTION: Interpolates the luminosity at a given observation time (t) using linear
 *              interpolation in log-space. The result is returned in linear space.
 ********************************************************************************************************************/

Real Observer::interpolate(InterpState const& state, size_t i, size_t j, size_t k, Real t_obs) const noexcept {
    Real log_t = fast_log2(t_obs / t_obs_grid(i, j, k));
    return surface(i, j, k) * fast_exp2(state.log_I_lo + log_t * state.slope);
}

/********************************************************************************************************************
 * METHOD: Observer::calc_emission_surface
 * DESCRIPTION: Calculates the observe frame effective emission surface for each effective phi and theta grid point.
 *              The emission surface is computed as the product of the differential cosine of theta and either 2π (if
 *              the effective phi size is 1) or the differential phi value.
 ********************************************************************************************************************/
void Observer::calc_emission_surface(Coord const& coord, Shock const& shock) {
    Array dphi({eff_phi_size}, 0);

    if (eff_phi_size == 1) {
        dphi(0) = 2 * con::pi;
    } else {
        int last = eff_phi_size - 1;
        for (int i = 0; i < eff_phi_size; ++i) {
            dphi(i) = 0.5 * (coord.phi(std::min(i + 1, last)) - coord.phi(std::max(i - 1, 0)));
        }
    }

    // precompute the dcos to avoid recomputing in axisymmetric jet
    static thread_local MeshGrid3d dcos;
    if (dcos.shape() != shock.theta.shape()) {
        dcos.resize(shock.theta.shape());
    }

    int last = theta_size - 1;
    size_t shock_phi_size = shock.theta.shape(0);
    for (size_t i = 0; i < shock_phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            size_t j_m1 = (j == 0) ? 0 : (j - 1);
            size_t j_p1 = (j == last) ? last : (j + 1);
            for (size_t k = 0; k < t_size; ++k) {
                Real theta_lo = 0.5 * (shock.theta(i, j, k) + shock.theta(i, j_m1, k));
                Real theta_hi = 0.5 * (shock.theta(i, j, k) + shock.theta(i, j_p1, k));
                dcos(i, j, k) = std::cos(theta_hi) - std::cos(theta_lo);
            }
        }
    }

    for (size_t i = 0; i < eff_phi_size; ++i) {
        size_t i_eff = i * jet_3d;
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < t_size; ++k) {
                if (shock.required(i_eff, j, k) == 0) {
                    continue;
                }

                Real dOmega = std::fabs(dcos(i_eff, j, k) * dphi(i));
                Real r = shock.r(i_eff, j, k);
                Real D = doppler(i, j, k);
                surface(i, j, k) = dOmega * r * r * D * D * D;
            }
        }
    }
}

/********************************************************************************************************************
 * METHOD: Observer::calc_t_obs_grid
 * DESCRIPTION: Calculates the observation time grid (t_obs_grid) and updates the doppler factor grid based on the
 *              Gamma (Lorentz factor) and engine time (t) array.
 *              For each grid point, the Doppler factor is computed and the observed time is calculated taking
 *              redshift into account.
 ********************************************************************************************************************/
void Observer::calc_t_obs_grid(Coord const& coord, Shock const& shock) {
    Real cos_obs = std::cos(coord.theta_view);
    Real sin_obs = std::sin(coord.theta_view);

    static thread_local MeshGrid3d cos_theta, sin_theta;
    if (cos_theta.shape() != shock.theta.shape()) {
        cos_theta.resize(shock.theta.shape());
        sin_theta.resize(shock.theta.shape());
    }

    size_t shock_phi_size = shock.theta.shape(0);
    for (size_t i = 0; i < shock_phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < t_size; ++k) {
                cos_theta(i, j, k) = std::cos(shock.theta(i, j, k));
                sin_theta(i, j, k) = std::sin(shock.theta(i, j, k));
            }
        }
    }

    for (size_t i = 0; i < eff_phi_size; ++i) {
        Real cos_phi = std::cos(coord.phi[i]);
        size_t i_eff = i * jet_3d;
        for (size_t j = 0; j < theta_size; ++j) {
            // Compute the cosine of the angle between the local velocity vector and the observer's line of sight.
            // Real cos_v = std::sin(coord.theta[j]) * cos_phi * sin_obs + std::cos(coord.theta[j]) * cos_obs;
            for (size_t k = 0; k < t_size; ++k) {
                Real gamma_ = shock.Gamma(i_eff, j, k);  // Get Gamma at the grid point.
                Real r = shock.r(i_eff, j, k);
                Real t_eng_ = coord.t(i_eff, j, k);  // Get engine time at the grid point.
                Real cos_v = sin_theta(i_eff, j, k) * cos_phi * sin_obs + cos_theta(i_eff, j, k) * cos_obs;
                // Compute the Doppler factor: D = 1 / [Gamma * (1 - beta * cos_v)]
                doppler(i, j, k) = 1 / (gamma_ - std::sqrt(gamma_ * gamma_ - 1) * cos_v);
                // Compute the observed time: t_obs = [t_eng + (1 - cos_v) * r / c] * (1 + z)
                t_obs_grid(i, j, k) = (t_eng_ + (1 - cos_v) * r / con::c) * one_plus_z;
            }
        }
    }
}

void Observer::update_required(MaskGrid& required, Array const& t_obs) {
    size_t t_obs_size = t_obs.size();

    // Loop over effective phi and theta grid points.
    for (size_t i = 0; i < eff_phi_size; i++) {
        size_t i_eff = i * jet_3d;
        for (size_t j = 0; j < theta_size; j++) {
            // Skip observation times that are below the grid's start time
            size_t t_idx = 0;
            iterate_to(t_obs_grid(i, j, 0), t_obs, t_idx);

            // find the grid points that are required for the interpolation.
            for (size_t k = 0; k < t_size - 1 && t_idx < t_obs_size; k++) {
                Real const t_lo = t_obs_grid(i, j, k);
                Real const t_hi = t_obs_grid(i, j, k + 1);

                if (t_lo <= t_obs(t_idx) && t_obs(t_idx) < t_hi) {
                    required(i_eff, j, k) = 1;
                    required(i_eff, j, k + 1) = 1;
                }

                iterate_to(t_hi, t_obs, t_idx);
            }
        }
    }
}

void Observer::build_t_obs_grid(Coord const& coord, Shock const& shock, Real luminosity_dist, Real redshift) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    // Determine if the jet is 3D (more than one phi value)
    jet_3d = static_cast<size_t>((phi_size > 1));
    // Set effective phi grid size based on the observation angle and jet dimensionality.
    if (coord.theta_view == 0 && jet_3d == 0) {
        this->eff_phi_size = 1;  // optimize for on-axis observer
    } else {
        this->eff_phi_size = coord.phi.size();
    }
    this->theta_size = theta_size;
    this->t_size = t_size;
    this->lumi_dist = luminosity_dist;
    this->one_plus_z = 1 + redshift;

    t_obs_grid.resize({eff_phi_size, theta_size, t_size});
    doppler.resize({eff_phi_size, theta_size, t_size});
    surface.resize({eff_phi_size, theta_size, t_size});

    // Calculate the solid angle grid and observation time grid.
    calc_t_obs_grid(coord, shock);
}

/********************************************************************************************************************
 * METHOD: Observer::observe
 * DESCRIPTION: Sets up the Observer with the given coordinate grid, shocks, observation angle, luminosity
 *              and redshift. It initializes the observation time and Doppler factor grids, as well as the
 *              emission surface.
 ********************************************************************************************************************/
void Observer::observe(Coord const& coord, Shock const& shock, Real luminosity_dist, Real redshift) {
    build_t_obs_grid(coord, shock, luminosity_dist, redshift);
    calc_emission_surface(coord, shock);
}

void Observer::observe_at(Array const& t_obs, Coord const& coord, Shock& shock, Real luminosity_dist, Real redshift) {
    build_t_obs_grid(coord, shock, luminosity_dist, redshift);

    xt::noalias(shock.required) = 0;

    update_required(shock.required, t_obs);

    calc_emission_surface(coord, shock);
}