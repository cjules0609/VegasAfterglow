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
 * FUNCTION: dec_radius
 * DESCRIPTION: Computes the deceleration radius of the shock.
 *              For a given isotropic energy E_iso, ISM density n_ism, initial Lorentz factor Gamma0,
 *              and engine duration, the deceleration radius is the maximum of the thin shell and thick shell
 *              deceleration radii.
 ********************************************************************************************************************/
Real dec_radius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura) {
    return std::max(thin_shell_dec_radius(E_iso, n_ism, Gamma0), thick_shell_dec_radius(E_iso, n_ism, engine_dura));
}

/********************************************************************************************************************
 * FUNCTION: thin_shell_dec_radius
 * DESCRIPTION: Computes the deceleration radius for the thin shell case using the formula:
 *                  R_dec = [3E_iso / (4π n_ism mp c^2 Gamma0^2)]^(1/3)
 ********************************************************************************************************************/
Real thin_shell_dec_radius(Real E_iso, Real n_ism, Real Gamma0) {
    return std::cbrt(3 * E_iso / (4 * con::pi * con::mp * con::c2 * n_ism * Gamma0 * Gamma0));
}

/********************************************************************************************************************
 * FUNCTION: thick_shell_dec_radius
 * DESCRIPTION: Computes the deceleration radius for the thick shell case using the formula:
 *                  R_dec = [3 E_iso engine_dura c / (4π n_ism mp c^2)]^(1/4)
 ********************************************************************************************************************/
Real thick_shell_dec_radius(Real E_iso, Real n_ism, Real engine_dura) {
    return std::sqrt(std::sqrt(3 * E_iso * engine_dura / n_ism * con::c / (4 * con::pi * con::mp * con::c2)));
}

/********************************************************************************************************************
 * FUNCTION: shell_spreading_radius
 * DESCRIPTION: Computes the radius at which shell spreading becomes significant.
 *              The formula is: R_spread = Gamma0^2 * c * engine_dura.
 ********************************************************************************************************************/
Real shell_spreading_radius(Real Gamma0, Real engine_dura) { return Gamma0 * Gamma0 * con::c * engine_dura; }

/********************************************************************************************************************
 * FUNCTION: RS_transition_radius
 * DESCRIPTION: Computes the radius at which the reverse shock transitions, based on the Sedov length,
 *              engine duration, and initial Lorentz factor.
 *              The formula is: R_RS = (SedovLength^(1.5)) / (sqrt(c * engine_dura) * Gamma0^2)
 ********************************************************************************************************************/
Real RS_transition_radius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura) {
    return std::pow(sedov_length(E_iso, n_ism), 1.5) / std::sqrt(con::c * engine_dura) / Gamma0 / Gamma0;
}

/********************************************************************************************************************
 * FUNCTION: shell_thickness_param
 * DESCRIPTION: Computes the dimensionless parameter (ξ) that characterizes the shell geometry.
 *              This parameter helps determine whether the shell behaves as thick or thin.
 *              The formula is: ξ = sqrt(Sedov_length / shell_width) * Gamma0^(-4/3)
 ********************************************************************************************************************/
Real shell_thickness_param(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura) {
    Real Sedov_l = sedov_length(E_iso, n_ism);
    Real shell_width = con::c * engine_dura;
    return std::sqrt(Sedov_l / shell_width) * std::pow(Gamma0, -4. / 3);
}

/********************************************************************************************************************
 * FUNCTION: calc_engine_duration
 * DESCRIPTION: Calculates the engine duration needed to achieve a specific shell thickness parameter (ξ)
 *              for given energy, density, and Lorentz factor values.
 *              The formula is: T_engine = Sedov_l / (ξ^2 * Gamma0^(8/3) * c)
 ********************************************************************************************************************/
Real calc_engine_duration(Real E_iso, Real n_ism, Real Gamma0, Real xi) {
    Real Sedov_l = sedov_length(E_iso, n_ism);
    return Sedov_l / (xi * xi * std::pow(Gamma0, 8. / 3) * con::c);
}
