//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/


#pragma once
#include <cmath>

#include "macros.h"

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the deceleration radius of the shock.
 * @details For a given isotropic energy E_iso, ISM density n_ism, initial Lorentz factor Gamma0,
 *          and engine duration, the deceleration radius is the maximum of the thin shell and thick shell
 *          deceleration radii.
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param Gamma0 Initial Lorentz factor
 * @param engine_dura Engine duration
 * @return The deceleration radius
 * <!-- ************************************************************************************** -->
 */
Real dec_radius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura);

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the deceleration radius for the thin shell case.
 * @details Uses the formula: R_dec = [3E_iso / (4π n_ism mp c^2 Gamma0^2)]^(1/3)
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param Gamma0 Initial Lorentz factor
 * @return The thin shell deceleration radius
 * <!-- ************************************************************************************** -->
 */
Real thin_shell_dec_radius(Real E_iso, Real n_ism, Real Gamma0);

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the deceleration radius for the thick shell case.
 * @details Uses the formula: R_dec = [3 E_iso engine_dura c / (4π n_ism mp c^2)]^(1/4)
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param engine_dura Engine duration
 * @return The thick shell deceleration radius
 * <!-- ************************************************************************************** -->
 */
Real thick_shell_dec_radius(Real E_iso, Real n_ism, Real engine_dura);

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the radius at which shell spreading becomes significant.
 * @details Uses the formula: R_spread = Gamma0^2 * c * engine_dura
 * @param Gamma0 Initial Lorentz factor
 * @param engine_dura Engine duration
 * @return The shell spreading radius
 * <!-- ************************************************************************************** -->
 */
Real shell_spreading_radius(Real Gamma0, Real engine_dura);

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the radius at which the reverse shock transitions.
 * @details Based on the Sedov length, engine duration, and initial Lorentz factor.
 *          Uses the formula: R_RS = (SedovLength^(1.5)) / (sqrt(c * engine_dura) * Gamma0^2)
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param Gamma0 Initial Lorentz factor
 * @param engine_dura Engine duration
 * @return The reverse shock transition radius
 * <!-- ************************************************************************************** -->
 */
Real RS_transition_radius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura);

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the dimensionless parameter (ξ) that characterizes the shell geometry.
 * @details This parameter helps determine whether the shell behaves as thick or thin.
 *          Uses the formula: ξ = sqrt(Sedov_length / shell_width) * Gamma0^(-4/3)
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param Gamma0 Initial Lorentz factor
 * @param engine_dura Engine duration
 * @return The shell thickness parameter ξ
 * <!-- ************************************************************************************** -->
 */
Real shell_thickness_param(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura);

/**
 * <!-- ************************************************************************************** -->
 * @brief Calculates the engine duration needed to achieve a specific shell thickness parameter.
 * @details Uses the formula: T_engine = Sedov_l / (ξ^2 * Gamma0^(8/3) * c)
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param Gamma0 Initial Lorentz factor
 * @param xi Target shell thickness parameter
 * @return The required engine duration
 * <!-- ************************************************************************************** -->
 */
Real calc_engine_duration(Real E_iso, Real n_ism, Real Gamma0, Real xi);

/**
 * <!-- ************************************************************************************** -->
 * @defgroup GammaConversions Gamma Conversion and Adiabatic Index Functions
 * @brief Helper functions for Lorentz factor conversions and adiabatic index calculations
 * <!-- ************************************************************************************** -->
 */

/**
 * <!-- ************************************************************************************** -->
 * @brief Converts Lorentz factor (gamma) to velocity fraction (beta)
 * @param gamma Lorentz factor
 * @return Velocity fraction (beta = v/c)
 * <!-- ************************************************************************************** -->
 */
inline Real gamma_to_beta(Real gamma) { return std::sqrt(1 - 1 / (gamma * gamma)); }

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes adiabatic index as a function of Lorentz factor
 * @param gamma Lorentz factor
 * @return Adiabatic index
 * <!-- ************************************************************************************** -->
 */
inline Real adiabatic_idx(Real gamma) { return (4 * gamma + 1) / (3 * gamma); }

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the Sedov length—a characteristic scale for blast wave deceleration
 * @param E_iso Isotropic equivalent energy
 * @param n_ism ISM number density
 * @return Sedov length
 * @details The Sedov length is a characteristic scale defined as the cube root of (E_iso / (ρc²)),
 *          where ρ is the ambient medium mass density
 * <!-- ************************************************************************************** -->
 */
inline Real sedov_length(Real E_iso, Real n_ism) {
    return std::cbrt(E_iso / (4 * con::pi / 3 * n_ism * con::mp * con::c2));
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Returns the radius at which the reverse shock crosses, defined as the thick shell deceleration radius
 * @param E_iso Isotropic equivalent energy
 * @param n_ism ISM number density
 * @param engine_dura Engine duration
 * @return Reverse shock crossing radius
 * <!-- ************************************************************************************** -->
 */
inline Real RS_crossing_radius(Real E_iso, Real n_ism, Real engine_dura) {
    Real l = sedov_length(E_iso, n_ism);
    return std::sqrt(std::sqrt(l * l * l * con::c * engine_dura));
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Determines the edge of the jet based on a given gamma cut-off using binary search
 * @tparam Ejecta Type of the jet/ejecta class
 * @param jet The jet/ejecta object
 * @param gamma_cut Lorentz factor cutoff value
 * @return Angle (in radians) at which the jet's Lorentz factor drops to gamma_cut
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta>
Real find_jet_edge(Ejecta const& jet, Real gamma_cut);

/**
 * <!-- ************************************************************************************** -->
 * @brief Determines the edge of the jet where the spreading is strongest
 * @tparam Ejecta Type of the jet/ejecta class
 * @tparam Medium Type of the ambient medium
 * @param jet The jet/ejecta object
 * @param medium The ambient medium object
 * @param phi Azimuthal angle
 * @param theta_min Minimum polar angle to consider
 * @param theta_max Maximum polar angle to consider
 * @param t0 Initial time
 * @return Angle (in radians) where the spreading is strongest
 * @details The spreading strength is measured by the derivative of the pressure with respect to theta,
 *          which is proportional to d((Gamma-1)Gamma rho)/dtheta
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium>
Real jet_spreading_edge(Ejecta const& jet, Medium const& medium, Real phi, Real theta_min, Real theta_max, Real t0);

//========================================================================================================
//                                  template function implementation
//========================================================================================================

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

template <typename Ejecta, typename Medium>
Real jet_spreading_edge(Ejecta const& jet, Medium const& medium, Real phi, Real theta_min, Real theta_max, Real t0) {
    Real step = (theta_max - theta_min) / 256;
    Real theta_s = theta_min;
    Real dp_min = 0;

    for (Real theta = theta_min; theta <= theta_max; theta += step) {
        // Real G = jet.Gamma0(phi, theta);
        // Real beta0 = gamma_to_beta(G);
        // Real r0 = beta0 * con::c * t0 / (1 - beta0);
        // Real rho = medium.rho(phi, theta, 0);
        Real th_lo = std::max(theta - step, theta_min);
        Real th_hi = std::min(theta + step, theta_max);
        Real dG = (jet.Gamma0(phi, th_hi) - jet.Gamma0(phi, th_lo)) / (th_hi - th_lo);
        // Real drho = (medium.rho(phi, th_hi, r0) - medium.rho(phi, th_lo, r0)) / (th_hi - th_lo);
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

/**
 * <!-- ************************************************************************************** -->
 * @brief Parameters for radiation transport
 * @details Parameters for radiation transport
 * <!-- ************************************************************************************** -->
 */
struct RadParams {
    Real eps_e{0.1};   ///< Electron energy fraction
    Real eps_B{0.01};  ///< Magnetic field energy fraction
    Real p{2.3};       ///< Electron energy distribution index
    Real xi_e{1};      ///< Electron self-absorption parameter
};