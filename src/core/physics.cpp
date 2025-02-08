//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "physics.h"

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