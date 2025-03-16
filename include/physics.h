//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _RELATIVITY_H_
#define _RELATIVITY_H_
#include <cmath>

#include "jet.h"
#include "macros.h"
#include "medium.h"
#include "mesh.h"
/********************************************************************************************************************
 * FUNCTION PROTOTYPES: Cosmological and Shock Distance Calculations
 * DESCRIPTION: These functions compute various distances used in shock and cosmological calculations.
 ********************************************************************************************************************/
Real decRadius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura);  // Deceleration radius.
Real thinShellDecRadius(Real E_iso, Real n_ism, Real Gamma0);           // Thin shell deceleration radius.
Real thickShellDecRadius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura);
Real shellSpreadingRadius(Real Gamma0, Real engine_dura);  // Shell spreading radius.
Real RSTransitionRadius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura);
Real shellRegimeParam(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura);
Real calcEngineDuration(Real E_iso, Real n_ism, Real Gamma0, Real xi);
/********************************************************************************************************************
 * INLINE FUNCTIONS: Gamma Conversion and Adiabatic Index
 * DESCRIPTION: Helper inline functions to convert a Lorentz factor (gamma) to the corresponding beta value
 *              and to compute the adiabatic index as a function of gamma.
 ********************************************************************************************************************/
inline Real gammaTobeta(Real gamma) {
    return std::sqrt(gamma * gamma - 1) / gamma;  // Convert Lorentz factor to velocity fraction (beta).
}
inline Real adiabaticIndex(Real gamma) {
    return (4 * gamma + 1) / (3 * gamma);  // Compute adiabatic index.
}

/********************************************************************************************************************
 * INLINE FUNCTION: SedovLength
 * DESCRIPTION: Computes the Sedov length—a characteristic scale for blast wave deceleration—given the
 *              isotropic equivalent energy (E_iso) and the ISM number density (n_ism).
 ********************************************************************************************************************/
inline Real SedovLength(Real E_iso, Real n_ism) {
    return std::cbrt(E_iso / (n_ism * con::mp * con::c2));
    // return std::cbrt(E_iso / (4 * con::pi / 3 * n_ism * con::mp * con::c2));
}

/********************************************************************************************************************
 * INLINE FUNCTION: RShockCrossingRadius
 * DESCRIPTION: Returns the radius at which the reverse shock crosses, defined as the thick shell deceleration radius.
 ********************************************************************************************************************/
inline Real RShockCrossingRadius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura) {
    size_t l = SedovLength(E_iso, n_ism);
    return std::sqrt(std::sqrt(l * l * l * con::c * engine_dura));
}

/********************************************************************************************************************
 * TEMPLATE FUNCTION: jetEdge
 * DESCRIPTION: Determines the edge of the jet based on a given gamma cut-off (gamma_cut) using binary search.
 *              It returns the angle (in radians) at which the jet's Lorentz factor drops to gamma_cut.
 ********************************************************************************************************************/
template <typename Ejecta>
Real jetEdge(Ejecta const& jet, Real gamma_cut) {
    if (jet.Gamma(0, con::pi / 2, 0) >= gamma_cut) {
        return con::pi / 2;  // If the Lorentz factor at pi/2 is above the cut, the jet extends to pi/2.
    }
    Real low = 0;
    Real hi = con::pi / 2;
    Real eps = 1e-9;
    for (; hi - low > eps;) {
        Real mid = 0.5 * (low + hi);
        if (jet.Gamma(0, mid, 0) > gamma_cut) {
            low = mid;
        } else {
            hi = mid;
        }
    }
    return low;
}

/********************************************************************************************************************
 * TEMPLATE FUNCTION: autoGrid
 * DESCRIPTION: Constructs a coordinate grid (Coord) for shock evolution. The grid is based on the
 *              observation times (t_obs), a maximum theta value, and specified numbers of grid points in phi,
 *              theta, and t. The radial grid is logarithmically spaced between t_min and t_max, and the theta grid
 *              is generated linearly.
 ********************************************************************************************************************/
template <typename Ejecta>
Coord autoGrid(Ejecta const& jet, Array const& t_obs, Real theta_cut, size_t phi_num = 32, size_t theta_num = 32,
               size_t t_num = 32, double theta_view_max = con::pi / 2) {
    Array phi = linspace(0, 2 * con::pi, phi_num);  // Generate phi grid linearly spaced.
    Real jet_edge = jetEdge(jet, con::Gamma_cut);   // Determine the jet edge angle.
    // Array theta = uniform_cos(0, std::min(jet_edge, theta_cut), theta_num);  // Generate theta grid uniformly in
    // cosine.
    Array theta = linspace(1e-6, std::min(jet_edge, theta_cut), theta_num);  // Generate theta grid uniformly
    Real theta_max = theta[theta_num - 1];                                   // Maximum theta value.
    Real t_max = *std::max_element(t_obs.begin(), t_obs.end());              // Maximum observation time.
    Real t_min = *std::min_element(t_obs.begin(), t_obs.end());              // Minimum observation time.
    Real b_max = gammaTobeta(jet.Gamma(0, 0, 0));                            // Maximum beta value.
    Real t_start =
        t_min * (1 - b_max) / (1 - std::cos(theta_max + theta_view_max) * b_max);  // Start time for the grid.
    t_start = std::min(t_start, 10 * con::sec);
    Real t_end = t_max;
    Array t = logspace(t_start, t_end, t_num);  // Generate logarithmically spaced radial grid.

    return Coord(phi, theta, t);  // Construct coordinate object.
}
#endif  // _RELATIVITY_H_
