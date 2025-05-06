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
 * @brief Computes the deceleration radius of the shock.
 * @details For a given isotropic energy E_iso, ISM density n_ism, initial Lorentz factor Gamma0,
 *          and engine duration, the deceleration radius is the maximum of the thin shell and thick shell
 *          deceleration radii.
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param Gamma0 Initial Lorentz factor
 * @param engine_dura Engine duration
 * @return The deceleration radius
 ********************************************************************************************************************/
Real dec_radius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura) {
    return std::max(thin_shell_dec_radius(E_iso, n_ism, Gamma0), thick_shell_dec_radius(E_iso, n_ism, engine_dura));
}

/********************************************************************************************************************
 * @brief Computes the deceleration radius for the thin shell case.
 * @details Uses the formula: R_dec = [3E_iso / (4π n_ism mp c^2 Gamma0^2)]^(1/3)
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param Gamma0 Initial Lorentz factor
 * @return The thin shell deceleration radius
 ********************************************************************************************************************/
Real thin_shell_dec_radius(Real E_iso, Real n_ism, Real Gamma0) {
    return std::cbrt(3 * E_iso / (4 * con::pi * con::mp * con::c2 * n_ism * Gamma0 * Gamma0));
}

/********************************************************************************************************************
 * @brief Computes the deceleration radius for the thick shell case.
 * @details Uses the formula: R_dec = [3 E_iso engine_dura c / (4π n_ism mp c^2)]^(1/4)
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param engine_dura Engine duration
 * @return The thick shell deceleration radius
 ********************************************************************************************************************/
Real thick_shell_dec_radius(Real E_iso, Real n_ism, Real engine_dura) {
    return std::sqrt(std::sqrt(3 * E_iso * engine_dura / n_ism * con::c / (4 * con::pi * con::mp * con::c2)));
}

/********************************************************************************************************************
 * @brief Computes the radius at which shell spreading becomes significant.
 * @details Uses the formula: R_spread = Gamma0^2 * c * engine_dura
 * @param Gamma0 Initial Lorentz factor
 * @param engine_dura Engine duration
 * @return The shell spreading radius
 ********************************************************************************************************************/
Real shell_spreading_radius(Real Gamma0, Real engine_dura) { return Gamma0 * Gamma0 * con::c * engine_dura; }

/********************************************************************************************************************
 * @brief Computes the radius at which the reverse shock transitions.
 * @details Based on the Sedov length, engine duration, and initial Lorentz factor.
 *          Uses the formula: R_RS = (SedovLength^(1.5)) / (sqrt(c * engine_dura) * Gamma0^2)
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param Gamma0 Initial Lorentz factor
 * @param engine_dura Engine duration
 * @return The reverse shock transition radius
 ********************************************************************************************************************/
Real RS_transition_radius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura) {
    return std::pow(sedov_length(E_iso, n_ism), 1.5) / std::sqrt(con::c * engine_dura) / Gamma0 / Gamma0;
}

/********************************************************************************************************************
 * @brief Computes the dimensionless parameter (ξ) that characterizes the shell geometry.
 * @details This parameter helps determine whether the shell behaves as thick or thin.
 *          Uses the formula: ξ = sqrt(Sedov_length / shell_width) * Gamma0^(-4/3)
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param Gamma0 Initial Lorentz factor
 * @param engine_dura Engine duration
 * @return The shell thickness parameter ξ
 ********************************************************************************************************************/
Real shell_thickness_param(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura) {
    Real Sedov_l = sedov_length(E_iso, n_ism);
    Real shell_width = con::c * engine_dura;
    return std::sqrt(Sedov_l / shell_width) * std::pow(Gamma0, -4. / 3);
}

/********************************************************************************************************************
 * @brief Calculates the engine duration needed to achieve a specific shell thickness parameter.
 * @details Uses the formula: T_engine = Sedov_l / (ξ^2 * Gamma0^(8/3) * c)
 * @param E_iso Isotropic energy
 * @param n_ism ISM density
 * @param Gamma0 Initial Lorentz factor
 * @param xi Target shell thickness parameter
 * @return The required engine duration
 ********************************************************************************************************************/
Real calc_engine_duration(Real E_iso, Real n_ism, Real Gamma0, Real xi) {
    Real Sedov_l = sedov_length(E_iso, n_ism);
    return Sedov_l / (xi * xi * std::pow(Gamma0, 8. / 3) * con::c);
}
