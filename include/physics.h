//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <cmath>

#include "jet.h"
#include "macros.h"
#include "medium.h"
#include "mesh.h"
/********************************************************************************************************************
 * FUNCTION PROTOTYPES: Cosmological and Shock Distance Calculations
 * DESCRIPTION: These functions compute various distances used in shock and cosmological calculations.
 ********************************************************************************************************************/
// Computes the deceleration radius, maximum of thin and thick shell cases
Real dec_radius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura);
// Computes deceleration radius for the thin shell case
Real thin_shell_dec_radius(Real E_iso, Real n_ism, Real Gamma0);
// Computes deceleration radius for the thick shell case
Real thick_shell_dec_radius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura);
// Computes the radius at which shell spreading becomes significant
Real shell_spreading_radius(Real Gamma0, Real engine_dura);
// Computes the radius at which reverse shock transitions
Real RS_transition_radius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura);
// Computes the shell thickness parameter (xi), characterizing shell geometry
Real shell_thickness_param(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura);
// Calculates engine duration based on other physical parameters and shell geometry
Real calc_engine_duration(Real E_iso, Real n_ism, Real Gamma0, Real xi);
/********************************************************************************************************************
 * INLINE FUNCTIONS: Gamma Conversion and Adiabatic Index
 * DESCRIPTION: Helper inline functions to convert a Lorentz factor (gamma) to the corresponding beta value
 *              and to compute the adiabatic index as a function of gamma.
 ********************************************************************************************************************/
// Converts Lorentz factor (gamma) to velocity fraction (beta)
inline Real gamma_to_beta(Real gamma) {
    return std::sqrt(1 - 1 / (gamma * gamma));  // Convert Lorentz factor to velocity fraction (beta).
}

// Computes adiabatic index as a function of Lorentz factor
inline Real adiabatic_idx(Real gamma) {
    return (4 * gamma + 1) / (3 * gamma);  // Compute adiabatic index.
}

/********************************************************************************************************************
 * INLINE FUNCTION: SedovLength
 * DESCRIPTION: Computes the Sedov length—a characteristic scale for blast wave deceleration—given the
 *              isotropic equivalent energy (E_iso) and the ISM number density (n_ism).
 ********************************************************************************************************************/
inline Real sedov_length(Real E_iso, Real n_ism) {
    // return std::cbrt(E_iso / (n_ism * con::mp * con::c2));
    return std::cbrt(E_iso / (4 * con::pi / 3 * n_ism * con::mp * con::c2));
}

/********************************************************************************************************************
 * INLINE FUNCTION: RShockCrossingRadius
 * DESCRIPTION: Returns the radius at which the reverse shock crosses, defined as the thick shell deceleration radius.
 ********************************************************************************************************************/
inline Real RS_crossing_radius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura) {
    size_t l = sedov_length(E_iso, n_ism);
    return std::sqrt(std::sqrt(l * l * l * con::c * engine_dura));
}

/********************************************************************************************************************
 * TEMPLATE FUNCTION: find_jet_edge
 * DESCRIPTION: Determines the edge of the jet based on a given gamma cut-off (gamma_cut) using binary search.
 *              It returns the angle (in radians) at which the jet's Lorentz factor drops to gamma_cut.
 ********************************************************************************************************************/
template <typename Ejecta>
Real find_jet_edge(Ejecta const& jet, Real gamma_cut) {
    if (jet.Gamma0(0, con::pi / 2) >= gamma_cut) {
        return con::pi / 2;  // If the Lorentz factor at pi/2 is above the cut, the jet extends to pi/2.
    }
    Real low = 0;
    Real hi = con::pi / 2;
    Real eps = 1e-9;
    for (; hi - low > eps;) {
        Real mid = 0.5 * (low + hi);
        if (jet.Gamma0(0, mid) > gamma_cut) {
            low = mid;
        } else {
            hi = mid;
        }
    }
    return low;
}

/********************************************************************************************************************
 * TEMPLATE FUNCTION: jet_spreading_edge
 * DESCRIPTION: Determines the edge of the jet where the spreading is strongest
 *              (dp/dtheta \propto d((Gamma-1)Gamma rho)/dtheta is largest).
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
Real jet_spreading_edge(Ejecta const& jet, Medium const& medium, Real phi, Real theta_min, Real theta_max, Real t0) {
    Real step = (theta_max - theta_min) / 256;
    Real theta_s = theta_min;
    Real dp_min = 0;

    for (Real theta = theta_min; theta <= theta_max; theta += step) {
        Real G = jet.Gamma0(phi, theta);
        Real beta0 = gamma_to_beta(G);
        Real r0 = beta0 * con::c * t0 / (1 - beta0);
        Real rho = medium.rho(phi, theta, 0);
        Real th_lo = std::max(theta - step, theta_min);
        Real th_hi = std::min(theta + step, theta_max);
        Real dG = (jet.Gamma0(phi, th_hi) - jet.Gamma0(phi, th_lo)) / (th_hi - th_lo);
        Real drho = (medium.rho(phi, th_hi, r0) - medium.rho(phi, th_lo, r0)) / (th_hi - th_lo);
        Real dp = dG;  //(2 * G - 1) * rho * dG + (G - 1) * G * drho;

        if (dp < dp_min) {
            dp_min = dp;
            theta_s = theta;
        }
    }
    if (dp_min == 0) {
        theta_s = theta_max;
    }

    return theta_s;
}

/********************************************************************************************************************
 * TEMPLATE FUNCTION: auto_grid
 * DESCRIPTION: Constructs a coordinate grid (Coord) for shock evolution. The grid is based on the
 *              observation times (t_obs), a maximum theta value, and specified numbers of grid points in phi,
 *              theta, and t. The radial grid is logarithmically spaced between t_min and t_max, and the theta grid
 *              is generated linearly.
 ********************************************************************************************************************/
template <typename Ejecta>
Coord auto_grid(Ejecta const& jet, Array const& t_obs, Real theta_cut, Real theta_view, Real z, size_t phi_num = 32,
                size_t theta_num = 32, size_t t_num = 32, bool is_axisymmetric = true) {
    Coord coord;
    coord.theta_view = theta_view;
    coord.phi = xt::linspace(0., 2 * con::pi, phi_num);  // Generate phi grid linearly spaced.
    Real jet_edge = find_jet_edge(jet, con::Gamma_cut);  // Determine the jet edge angle.
    // Array theta = uniform_cos(0, std::min(jet_edge, theta_cut), theta_num);  // Generate theta grid uniformly in
    // cosine.
    coord.theta = xt::linspace(1e-4_r, std::min(jet_edge, theta_cut), theta_num);  // Generate theta grid uniformly

    Real t_max = *std::max_element(t_obs.begin(), t_obs.end());  // Maximum observation time.
    Real t_min = *std::min_element(t_obs.begin(), t_obs.end());  // Minimum observation time.
    size_t phi_size_needed = is_axisymmetric ? 1 : phi_num;
    coord.t = xt::zeros<Real>({phi_size_needed, theta_num, t_num});
    for (size_t i = 0; i < phi_size_needed; ++i) {
        for (size_t j = 0; j < theta_num; ++j) {
            Real b = gamma_to_beta(jet.Gamma0(coord.phi(i), coord.theta(j)));
            Real theta_max = coord.theta(j) + theta_view;

            Real t_start = 0.99 * t_min * (1 - b) / (1 - std::cos(theta_max) * b) / (1 + z);
            Real t_end = 1.01 * t_max / (1 + z);
            xt::view(coord.t, i, j, xt::all()) = xt::logspace(std::log10(t_start), std::log10(t_end), t_num);
        }
    }

    return coord;  // Construct coordinate object.
}
