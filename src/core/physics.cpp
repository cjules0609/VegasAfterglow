//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "physics.h"

#include "mesh.h"
#include "shock.h"
#include "utilities.h"
/********************************************************************************************************************
 * FUNCTION: decRadius
 * DESCRIPTION: Computes the deceleration radius of the shock.
 *              For a given isotropic energy E_iso, ISM density n_ism, initial Lorentz factor Gamma0,
 *              and engine duration, the deceleration radius is the maximum of the thin shell and thick shell
 *              deceleration radii.
 ********************************************************************************************************************/
Real decRadius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura) {
    return std::max(thinShellDecRadius(E_iso, n_ism, Gamma0), thickShellDecRadius(E_iso, n_ism, Gamma0, engine_dura));
}

/********************************************************************************************************************
 * FUNCTION: thinShellDecRadius
 * DESCRIPTION: Computes the deceleration radius for the thin shell case using the formula:
 *                  R_dec = [3E_iso / (4π n_ism mp c^2 Gamma0^2)]^(1/3)
 ********************************************************************************************************************/
Real thinShellDecRadius(Real E_iso, Real n_ism, Real Gamma0) {
    return std::cbrt(3 * E_iso / (4 * con::pi * n_ism * con::mp * con::c2 * Gamma0 * Gamma0));
}

/********************************************************************************************************************
 * FUNCTION: thickShellDecRadius
 * DESCRIPTION: Computes the deceleration radius for the thick shell case using the formula:
 *                  R_dec = [3 E_iso engine_dura c / (4π n_ism mp c^2)]^(1/4)
 ********************************************************************************************************************/
Real thickShellDecRadius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura) {
    return std::sqrt(std::sqrt(3 * E_iso * engine_dura * con::c / (4 * con::pi * n_ism * con::mp * con::c2)));
}

/********************************************************************************************************************
 * FUNCTION: shellSpreadingRadius
 * DESCRIPTION: Computes the radius at which shell spreading becomes significant.
 *              The formula is: R_spread = Gamma0^2 * c * engine_dura.
 ********************************************************************************************************************/
Real shellSpreadingRadius(Real Gamma0, Real engine_dura) { return Gamma0 * Gamma0 * con::c * engine_dura; }

/********************************************************************************************************************
 * FUNCTION: RSTransitionRadius
 * DESCRIPTION: Computes the radius at which the reverse shock transitions, based on the Sedov length,
 *              engine duration, and initial Lorentz factor.
 *              The formula is: R_RS = (SedovLength^(1.5)) / (sqrt(c * engine_dura) * Gamma0^2)
 ********************************************************************************************************************/
Real RSTransitionRadius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura) {
    return std::pow(SedovLength(E_iso, n_ism), 1.5) / std::sqrt(con::c * engine_dura) / Gamma0 / Gamma0;
}

/********************************************************************************************************************
 * FUNCTION: shellRegimeParam
 * DESCRIPTION:
 ********************************************************************************************************************/
Real shellRegimeParam(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura) {
    Real Sedov_l = SedovLength(E_iso, n_ism);
    Real shell_width = con::c * engine_dura;
    return std::sqrt(Sedov_l / shell_width) * std::pow(Gamma0, -4. / 3);
}

/********************************************************************************************************************
 * FUNCTION: calcEngineDuration
 * DESCRIPTION:
 ********************************************************************************************************************/
Real calcEngineDuration(Real E_iso, Real n_ism, Real Gamma0, Real xi) {
    Real Sedov_l = SedovLength(E_iso, n_ism);
    return Sedov_l / (xi * xi * std::pow(Gamma0, 8. / 3) * con::c);
}

/********************************************************************************************************************
 * TEMPLATE FUNCTION: jetEdge
 * DESCRIPTION: Determines the edge of the jet based on a given gamma cut-off (gamma_cut) using binary search.
 *              It returns the angle (in radians) at which the jet's Lorentz factor drops to gamma_cut.
 ********************************************************************************************************************/
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
 * TEMPLATE FUNCTION: autoGrid
 * DESCRIPTION: Constructs a coordinate grid (Coord) for shock evolution. The grid is based on the
 *              observation times (t_obs), a maximum theta value, and specified numbers of grid points in phi,
 *              theta, and t. The radial grid is logarithmically spaced between t_min and t_max, and the theta grid
 *              is generated linearly.
 ********************************************************************************************************************/
Coord autoGrid(Ejecta const& jet, Array const& t_obs, Real theta_cut, Real theta_view_max, size_t phi_num,
               size_t theta_num, size_t t_num) {
    Array phi = xt::linspace(0., 2 * con::pi, phi_num);  // Generate phi grid linearly spaced.
    Real jet_edge = jetEdge(jet, con::Gamma_cut);        // Determine the jet edge angle.
    // Array theta = uniform_cos(0, std::min(jet_edge, theta_cut), theta_num);  // Generate theta grid uniformly in
    // cosine.
    Array theta = xt::linspace(1e-4, std::min(jet_edge, theta_cut), theta_num);  // Generate theta grid uniformly
    Real theta_max = theta.back();                                               // Maximum theta value.
    Real t_max = *std::max_element(t_obs.begin(), t_obs.end());                  // Maximum observation time.
    Real t_min = *std::min_element(t_obs.begin(), t_obs.end());                  // Minimum observation time.
    t_min = std::min(t_min, 10 * con::sec);
    Real b_max = gammaTobeta(jet.Gamma0(0, 0));  // Maximum beta value.
    Real t_start =
        t_min * (1 - b_max) / (1 - std::cos(theta_max + theta_view_max) * b_max);  // Start time for the grid.

    Real t_end = t_max;
    Array t =
        xt::logspace(std::log10(t_start), std::log10(t_end), t_num);  // Generate logarithmically spaced radial grid.

    return Coord(phi, theta, t);  // Construct coordinate object.
}

void autoGrid(Coord& coord, Ejecta const& jet, Array const& t_obs, Real theta_cut, Real theta_view_max, size_t phi_num,
              size_t theta_num, size_t t_num) {
    coord.phi = xt::linspace(0., 2 * con::pi, phi_num);  // Generate phi grid linearly spaced.
    Real jet_edge = jetEdge(jet, con::Gamma_cut);        // Determine the jet edge angle.
    // Array theta = uniform_cos(0, std::min(jet_edge, theta_cut), theta_num);  // Generate theta grid uniformly in
    // cosine.
    coord.theta = xt::linspace(1e-4, std::min(jet_edge, theta_cut), theta_num);  // Generate theta grid uniformly
    Real theta_max = coord.theta.back();                                         // Maximum theta value.
    Real t_max = *std::max_element(t_obs.begin(), t_obs.end());                  // Maximum observation time.
    Real t_min = *std::min_element(t_obs.begin(), t_obs.end());                  // Minimum observation time.
    t_min = std::min(t_min, 10 * con::sec);
    Real b_max = gammaTobeta(jet.Gamma0(0, 0));  // Maximum beta value.
    Real t_start =
        t_min * (1 - b_max) / (1 - std::cos(theta_max + theta_view_max) * b_max);  // Start time for the grid.

    Real t_end = t_max;
    coord.t =
        xt::logspace(std::log10(t_start), std::log10(t_end), t_num);  // Generate logarithmically spaced radial grid.
}
