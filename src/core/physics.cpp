//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "physics.h"

#include <boost/numeric/odeint.hpp>

#include "shock.h"
#include "utilities.h"

/********************************************************************************************************************
 * FUNCTION: zToLuminosityDistance
 * DESCRIPTION: Computes the luminosity distance corresponding to a given redshift z.
 *              The integration is performed using boost::numeric::odeint with a dense output stepper.
 *              The differential equation solved is:
 *                  dL/dz = 1 / sqrt(Ω_m (1+z)³ + Ω_Λ)
 *              and the final luminosity distance is given by:
 *                  D_L = (1+z) * L * c / H0
 ********************************************************************************************************************/
Real zToLuminosityDistance(Real z) {
    using namespace boost::numeric::odeint;
    Real atol = 0;
    Real rtol = 1e-6;
    auto stepper = make_dense_output(atol, rtol, runge_kutta_dopri5<Real>());

    // Define the differential equation: dL/dz = 1 / sqrt(Ω_m (1+z)^3 + Ω_Λ)
    auto eqn = [&](Real const& y, Real& dydz, Real z0) {
        dydz = 1 / std::sqrt(con::Omega_m * (1 + z0) * (1 + z0) * (1 + z0) + con::Omega_L);
    };
    // Initialize the stepper with initial condition L(0)=0 and step size 1e-6
    stepper.initialize(0, 0, 1e-6);

    // Integrate the ODE until z is reached
    for (; stepper.current_time() <= z;) {
        stepper.do_step(eqn);
    }
    Real L = 0;
    // Calculate the state at z
    stepper.calc_state(z, L);
    // Return the luminosity distance using the relation D_L = (1+z) * L * c / H0
    return (1 + z) * L * con::c / con::H0;
}

/********************************************************************************************************************
 * FUNCTION: luminosityDistanceToz
 * DESCRIPTION: Computes the redshift z corresponding to a given luminosity distance L.
 *              The integration is performed using boost::numeric::odeint with a dense output stepper,
 *              solving the same differential equation as zToLuminosityDistance.
 *              The integration continues until (1+z) * L_current * c / H0 exceeds L, at which point z is returned.
 ********************************************************************************************************************/
Real luminosityDistanceToz(Real L) {
    using namespace boost::numeric::odeint;
    Real atol = 0;
    Real rtol = 1e-6;
    auto stepper = make_dense_output(atol, rtol, runge_kutta_dopri5<Real>());

    // Define the same differential equation: dL/dz = 1 / sqrt(Ω_m (1+z)^3 + Ω_Λ)
    auto eqn = [&](Real const& y, Real& dydz, Real z0) {
        dydz = 1 / std::sqrt(con::Omega_m * (1 + z0) * (1 + z0) * (1 + z0) + con::Omega_L);
    };
    stepper.initialize(0, 0, 1e-6);
    // Integrate indefinitely until the condition is met
    for (;;) {
        Real z = stepper.current_time();
        Real L_current = stepper.current_state();
        // If the computed luminosity distance exceeds the given L, return the current redshift z
        if ((1 + z) * L_current * con::c / con::H0 >= L) {
            return z;
        }
        stepper.do_step(eqn);
    }
}

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
    return std::pow(3 * E_iso / (4 * con::pi * n_ism * con::mp * con::c2 * Gamma0 * Gamma0), 1.0 / 3);
}

/********************************************************************************************************************
 * FUNCTION: thickShellDecRadius
 * DESCRIPTION: Computes the deceleration radius for the thick shell case using the formula:
 *                  R_dec = [3 E_iso engine_dura c / (4π n_ism mp c^2)]^(1/4)
 ********************************************************************************************************************/
Real thickShellDecRadius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura) {
    return std::pow(3 * E_iso * engine_dura * con::c / (4 * con::pi * n_ism * con::mp * con::c2), 0.25);
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