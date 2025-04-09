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
    return std::sqrt(1 - 1 / (gamma * gamma));  // Convert Lorentz factor to velocity fraction (beta).
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
    // return std::cbrt(E_iso / (n_ism * con::mp * con::c2));
    return std::cbrt(E_iso / (4 * con::pi / 3 * n_ism * con::mp * con::c2));
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
 * TEMPLATE FUNCTION: jetSpreadingEdge
 * DESCRIPTION: Determines the edge of the jet where the spreading is strongest
 *              (dp/dtheta \propto d((Gamma-1)Gamma rho)/dtheta is largest).
 ********************************************************************************************************************/
template <typename Ejecta>
Real jetSpreadingEdge(Ejecta const& jet, Medium const& medium, Real phi, Real theta_min, Real theta_max, Real t0) {
    Real step = (theta_max - theta_min) / 256;
    Real theta_s = theta_min;
    Real dp_min = 0;

    for (double theta = theta_min; theta <= theta_max; theta += step) {
        Real G = jet.Gamma0(phi, theta);
        Real beta0 = gammaTobeta(G);
        Real r0 = beta0 * con::c * t0 / (1 - beta0);
        Real rho = medium.rho(phi, theta, 0);
        Real th_lo = std::max(theta - step, theta_min);
        Real th_hi = std::min(theta + step, theta_max);
        Real dG = (jet.Gamma0(phi, th_hi) - jet.Gamma0(phi, th_lo)) / (th_hi - th_lo);
        Real drho = (medium.rho(phi, th_hi, r0) - medium.rho(phi, th_lo, r0)) / (th_hi - th_lo);
        Real dp = (2 * G - 1) * rho * dG + (G - 1) * G * drho;

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

template <typename Ejecta>
Real jetSpreadingEdge(Ejecta const& jet, Real phi, Real theta_min, Real theta_max) {
    Real step = (theta_max - theta_min) / 256;
    Real theta_s = theta_min;
    Real dG_min = 0;

    for (double theta = theta_min; theta <= theta_max; theta += step) {
        Real th_lo = std::max(theta - step, theta_min);
        Real th_hi = std::min(theta + step, theta_max);
        Real dG = (jet.Gamma(phi, th_hi, 0) - jet.Gamma(phi, th_lo, 0)) / (th_hi - th_lo);

        if (dG < dG_min) {
            dG_min = dG;
            theta_s = theta;
        }
    }

    if (dG_min == 0) {
        theta_s = theta_max;
    }

    return theta_s;
}
/********************************************************************************************************************
 * TEMPLATE FUNCTION: autoGrid
 * DESCRIPTION: Constructs a coordinate grid (Coord) for shock evolution. The grid is based on the
 *              observation times (t_obs), a maximum theta value, and specified numbers of grid points in phi,
 *              theta, and t. The radial grid is logarithmically spaced between t_min and t_max, and the theta grid
 *              is generated linearly.
 ********************************************************************************************************************/
template <typename Ejecta>
Coord autoGrid(Ejecta const& jet, Array const& t_obs, Real theta_cut, Real theta_view_max, size_t phi_num = 32,
               size_t theta_num = 32, size_t t_num = 32) {
    Array phi = linspace(0, 2 * con::pi, phi_num);  // Generate phi grid linearly spaced.
    Real jet_edge = jetEdge(jet, con::Gamma_cut);   // Determine the jet edge angle.
    // Array theta = uniform_cos(0, std::min(jet_edge, theta_cut), theta_num);  // Generate theta grid uniformly in
    // cosine.
    Array theta = linspace(1e-4, std::min(jet_edge, theta_cut), theta_num);  // Generate theta grid uniformly
    Real theta_max = theta[theta_num - 1];                                   // Maximum theta value.
    Real t_max = *std::max_element(t_obs.begin(), t_obs.end());              // Maximum observation time.
    Real t_min = *std::min_element(t_obs.begin(), t_obs.end());              // Minimum observation time.
    t_min = std::min(t_min, 10 * con::sec);
    Real b_max = gammaTobeta(jet.Gamma0(0, 0));  // Maximum beta value.
    Real t_start =
        t_min * (1 - b_max) / (1 - std::cos(theta_max + theta_view_max) * b_max);  // Start time for the grid.

    /*Real R0 = 1e10 * con::cm;
    Real Rc = std::max(R0, jet.T0 * con::c);
    Real gamma_max = jet.Gamma0(0, 0);
    Real b_max = gammaTobeta(gamma_max);
    Real t_start = Rc * gamma_max * (1 - b_max) / (b_max * con::c);*/
    Real t_end = t_max;
    Array t = logspace(t_start, t_end, t_num);  // Generate logarithmically spaced radial grid.

    return Coord(phi, theta, t);  // Construct coordinate object.
}
#endif  // _RELATIVITY_H_
