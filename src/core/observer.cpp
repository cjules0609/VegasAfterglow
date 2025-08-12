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

Real Observer::interpolate(InterpState const& state, size_t i, size_t j, size_t k, Real lg2_t_obs) const noexcept {
    Real dlg2_t = lg2_t_obs - lg2_t(i, j, k);
    return fast_exp2(state.lg2_P_nu_lo + dlg2_t * state.slope);
}

void Observer::calc_obs_solid_angle(Coord const& coord, Shock const& shock) {
    Array dphi({eff_phi_grid}, 0);

    if (eff_phi_grid == 1) {
        dphi(0) = 2 * con::pi;
    } else {
        int last = eff_phi_grid - 1;
        for (int i = 0; i < eff_phi_grid; ++i) {
            dphi(i) = 0.5 * (coord.phi(std::min(i + 1, last)) - coord.phi(std::max(i - 1, 0)));
        }
    }

    // precompute the dcos to avoid recomputing in axisymmetric jet
    static thread_local MeshGrid3d dcos;
    if (dcos.shape() != shock.theta.shape()) {
        dcos.resize(shock.theta.shape());
    }

    int last = theta_grid - 1;
    size_t shock_phi_size = shock.theta.shape(0);
    for (size_t i = 0; i < shock_phi_size; ++i) {
        for (size_t j = 0; j < theta_grid; ++j) {
            size_t j_p1 = (j == last) ? last : (j + 1);
            for (size_t k = 0; k < t_grid; ++k) {
                Real theta_lo = (j == 0) ? 0.0 : 0.5 * (shock.theta(i, j, k) + shock.theta(i, j - 1, k));
                Real theta_hi = 0.5 * (shock.theta(i, j, k) + shock.theta(i, j_p1, k));
                dcos(i, j, k) = std::cos(theta_hi) - std::cos(theta_lo);
            }
        }
    }

    for (size_t i = 0; i < eff_phi_grid; ++i) {
        size_t i_eff = i * jet_3d;
        for (size_t j = 0; j < theta_grid; ++j) {
            for (size_t k = 0; k < t_grid; ++k) {
                // if (shock.required(i_eff, j, k) == 0) {
                //     continue;
                // }  // maybe remove this inner branch harm to vectorization
                Real dOmega = std::fabs(dcos(i_eff, j, k) * dphi(i));
                lg2_Omega(i, j, k) = std::log2(dOmega);
            }
        }
    }
}

void Observer::calc_t_obs(Coord const& coord, Shock const& shock) {
    Real cos_obs = std::cos(coord.theta_view);
    Real sin_obs = std::sin(coord.theta_view);

    static thread_local MeshGrid3d cos_theta, sin_theta;
    if (cos_theta.shape() != shock.theta.shape()) {
        cos_theta.resize(shock.theta.shape());
        sin_theta.resize(shock.theta.shape());
    }

    size_t shock_phi_size = shock.theta.shape(0);
    for (size_t i = 0; i < shock_phi_size; ++i) {
        for (size_t j = 0; j < theta_grid; ++j) {
            for (size_t k = 0; k < t_grid; ++k) {
                cos_theta(i, j, k) = std::cos(shock.theta(i, j, k));
                sin_theta(i, j, k) = std::sin(shock.theta(i, j, k));
            }
        }
    }

    for (size_t i = 0; i < eff_phi_grid; ++i) {
        Real cos_phi = std::cos(coord.phi[i] - coord.phi_view);
        size_t i_eff = i * jet_3d;
        for (size_t j = 0; j < theta_grid; ++j) {
            // Compute the cosine of the angle between the local velocity vector and the observer's line of sight.
            for (size_t k = 0; k < t_grid; ++k) {
                Real gamma_ = shock.Gamma(i_eff, j, k);
                Real r = shock.r(i_eff, j, k);
                Real t_eng_ = coord.t(i_eff, j, k);
                Real cos_v = sin_theta(i_eff, j, k) * cos_phi * sin_obs + cos_theta(i_eff, j, k) * cos_obs;
                // Compute the Doppler factor: D = 1 / [Gamma * (1 - beta * cos_v)]
                lg2_doppler(i, j, k) = -std::log2((gamma_ - std::sqrt(gamma_ * gamma_ - 1) * cos_v));
                // Compute the observed time: t_obs = [t_eng + (1 - cos_v) * r / c] * (1 + z)
                time(i, j, k) = (t_eng_ + (1 - cos_v) * r / con::c) * one_plus_z;
            }
        }
    }
    lg2_t = xt::log2(time);
}

void Observer::update_required(MaskGrid& required, Array const& t_obs) {
    size_t t_obs_size = t_obs.size();
    Array log2_t_obs = xt::log2(t_obs);
    // Loop over effective phi and theta grid points.
    for (size_t i = 0; i < eff_phi_grid; i++) {
        size_t i_eff = i * jet_3d;
        for (size_t j = 0; j < theta_grid; j++) {
            // Skip observation times that are below the grid's start time
            size_t t_idx = 0;
            iterate_to(lg2_t(i, j, 0), log2_t_obs, t_idx);
            // find the grid points that are required for the interpolation.
            for (size_t k = 0; k < t_grid - 1 && t_idx < t_obs_size; k++) {
                Real const t_lo = lg2_t(i, j, k);
                Real const t_hi = lg2_t(i, j, k + 1);

                if (t_lo < log2_t_obs(t_idx) && log2_t_obs(t_idx) <= t_hi) {
                    required(i_eff, j, k) = 1;
                    required(i_eff, j, k + 1) = 1;
                }
                iterate_to(t_hi, log2_t_obs, t_idx);
            }
        }
    }
}

void Observer::build_time_grid(Coord const& coord, Shock const& shock, Real luminosity_dist, Real redshift) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    // Determine if the jet is 3D (more than one phi value)
    jet_3d = static_cast<size_t>((phi_size > 1));
    // Set effective phi grid size based on the observation angle and jet dimensionality.
    if (coord.theta_view == 0 && jet_3d == 0) {
        this->eff_phi_grid = 1;  // optimize for on-axis observer
    } else {
        this->eff_phi_grid = coord.phi.size();
    }
    this->theta_grid = theta_size;
    this->t_grid = t_size;
    this->lumi_dist = luminosity_dist;
    this->one_plus_z = 1 + redshift;
    this->lg2_one_plus_z = std::log2(one_plus_z);

    time.resize({eff_phi_grid, theta_size, t_size});
    lg2_t.resize({eff_phi_grid, theta_size, t_size});
    lg2_doppler.resize({eff_phi_grid, theta_size, t_size});
    lg2_Omega.resize({eff_phi_grid, theta_size, t_size});

    // Calculate the solid angle grid and observation time grid.
    calc_t_obs(coord, shock);
}

void Observer::observe(Coord const& coord, Shock const& shock, Real luminosity_dist, Real redshift) {
    build_time_grid(coord, shock, luminosity_dist, redshift);
    calc_obs_solid_angle(coord, shock);
}

void Observer::observe_at(Array const& t_obs, Coord const& coord, Shock& shock, Real luminosity_dist, Real redshift) {
    build_time_grid(coord, shock, luminosity_dist, redshift);
    xt::view(shock.required, xt::all(), xt::all(), xt::all()) = 0;
    update_required(shock.required, t_obs);
    calc_obs_solid_angle(coord, shock);
}