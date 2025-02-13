//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "jet.h"

#include <cmath>
#include <iostream>

#include "macros.h"
#include "physics.h"
#include "utilities.h"

/********************************************************************************************************************
 * FUNCTION: LiangGhirlanda2010_
 * DESCRIPTION: Computes the Lorentz factor directly from an energy value 'e' using the Liang & Ghirlanda (2010)
 *              prescription. The calculation is:
 *                  u = (e / e_max)^idx * gamma_max
 *                  return sqrt(1 + u^2)
 ********************************************************************************************************************/
Real LiangGhirlanda2010_(Real e, Real e_max, Real gamma_max, Real idx) {
    Real u = fastPow(e / e_max, idx) * gamma_max;
    return std::sqrt(1 + u * u);
}

/********************************************************************************************************************
 * METHOD: TophatJet::dEdOmega
 * DESCRIPTION: Returns the energy per unit solid angle for a TophatJet. For angles less than the core angle
 *              theta_c_, the energy per unit solid angle is dEdOmega_; otherwise, it is 0.
 ********************************************************************************************************************/
Real TophatJet::dEdOmega(Real phi, Real theta, Real t) const { return theta < theta_c_ ? dEdOmega_ : 0; }

/********************************************************************************************************************
 * METHOD: TophatJet::Gamma0
 * DESCRIPTION: Returns the initial Lorentz factor for a TophatJet. For angles less than the core angle theta_c_,
 *              the Lorentz factor is Gamma0_; otherwise, it is 0.
 ********************************************************************************************************************/
Real TophatJet::Gamma0(Real phi, Real theta, Real t) const { return theta < theta_c_ ? Gamma0_ : 0; }

/********************************************************************************************************************
 * METHOD: GaussianJet::dEdOmega
 * DESCRIPTION: Returns the energy per unit solid angle for a GaussianJet. The energy is attenuated
 *              according to a Gaussian profile based on theta relative to the core angle theta_c_.
 ********************************************************************************************************************/
Real GaussianJet::dEdOmega(Real phi, Real theta, Real t) const {
    return dEdOmega_ * fastExp(-theta * theta / two_theta_c_sq);
}

/********************************************************************************************************************
 * METHOD: GaussianJet::Gamma0
 * DESCRIPTION: Computes the initial Lorentz factor for a GaussianJet using the LiangGhirlanda2010_ helper.
 *              The function uses the computed dEdOmega (from the Gaussian profile), along with the jet's
 *              dEdOmega_ and Gamma0_ parameters and index idx_.
 ********************************************************************************************************************/
Real GaussianJet::Gamma0(Real phi, Real theta, Real t) const {
    return LiangGhirlanda2010_(dEdOmega(phi, theta, t), dEdOmega_, Gamma0_, idx_);
}

/********************************************************************************************************************
 * METHOD: PowerLawJet::dEdOmega
 * DESCRIPTION: Returns the energy per unit solid angle for a PowerLawJet. For angles beyond the core (theta >=
 *theta_c_), the energy falls off as a power law with index k_.
 ********************************************************************************************************************/
Real PowerLawJet::dEdOmega(Real phi, Real theta, Real t) const {
    return dEdOmega_ * fastPow(1 + theta / theta_c_, -k_);
}

/********************************************************************************************************************
 * METHOD: PowerLawJet::Gamma0
 * DESCRIPTION: Computes the initial Lorentz factor for a PowerLawJet using the LiangGhirlanda2010_ helper.
 *              The function uses the computed dEdOmega (from the power-law profile), along with the jet's
 *              dEdOmega_ and Gamma0_ parameters and index idx_.
 ********************************************************************************************************************/
Real PowerLawJet::Gamma0(Real phi, Real theta, Real t) const {
    return LiangGhirlanda2010_(dEdOmega(phi, theta, t), dEdOmega_, Gamma0_, idx_);
}