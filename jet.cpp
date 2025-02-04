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
 * FUNCTION: LiangGhirlanda2010
 * DESCRIPTION: Returns a lambda function that computes a Lorentz factor using the Liang & Ghirlanda (2010)
 *              prescription. The lambda takes an energy function (energy_func) and parameters e_max, gamma_max,
 *              and idx, and returns a function of (phi, theta, t). The returned function calculates:
 *                  e = energy_func(phi, theta, 0)
 *                  u = (e / e_max)^idx * gamma_max
 *                  return sqrt(1 + u^2)
 ********************************************************************************************************************/
auto LiangGhirlanda2010(auto energy_func, double e_max, double gamma_max, double idx) {
    return [=](double phi, double theta, double t = 0) -> double {
        double e = energy_func(phi, theta, 0);
        double u = std::pow(e / e_max, idx) * gamma_max;
        return std::sqrt(1 + u * u);
    };
}

/********************************************************************************************************************
 * FUNCTION: LiangGhirlanda2010_
 * DESCRIPTION: Computes the Lorentz factor directly from an energy value 'e' using the Liang & Ghirlanda (2010)
 *              prescription. The calculation is:
 *                  u = (e / e_max)^idx * gamma_max
 *                  return sqrt(1 + u^2)
 ********************************************************************************************************************/
double LiangGhirlanda2010_(double e, double e_max, double gamma_max, double idx) {
    double u = std::pow(e / e_max, idx) * gamma_max;
    return std::sqrt(1 + u * u);
}

/********************************************************************************************************************
 * METHOD: TophatJet::dEdOmega
 * DESCRIPTION: Returns the energy per unit solid angle for a TophatJet. For angles less than the core angle
 *              theta_c_, the energy per unit solid angle is dEdOmega_; otherwise, it is 0.
 ********************************************************************************************************************/
double TophatJet::dEdOmega(double phi, double theta, double t) const { return theta < theta_c_ ? dEdOmega_ : 0; }

/********************************************************************************************************************
 * METHOD: TophatJet::Gamma0
 * DESCRIPTION: Returns the initial Lorentz factor for a TophatJet. For angles less than the core angle theta_c_,
 *              the Lorentz factor is Gamma0_; otherwise, it is 0.
 ********************************************************************************************************************/
double TophatJet::Gamma0(double phi, double theta, double t) const { return theta < theta_c_ ? Gamma0_ : 0; }

/********************************************************************************************************************
 * METHOD: GaussianJet::dEdOmega
 * DESCRIPTION: Returns the energy per unit solid angle for a GaussianJet. The energy is attenuated
 *              according to a Gaussian profile based on theta relative to the core angle theta_c_.
 ********************************************************************************************************************/
double GaussianJet::dEdOmega(double phi, double theta, double t) const {
    return dEdOmega_ * fastExp(-theta * theta / (2 * theta_c_ * theta_c_));
}

/********************************************************************************************************************
 * METHOD: GaussianJet::Gamma0
 * DESCRIPTION: Computes the initial Lorentz factor for a GaussianJet using the LiangGhirlanda2010_ helper.
 *              The function uses the computed dEdOmega (from the Gaussian profile), along with the jet's
 *              dEdOmega_ and Gamma0_ parameters and index idx_.
 ********************************************************************************************************************/
double GaussianJet::Gamma0(double phi, double theta, double t) const {
    return LiangGhirlanda2010_(dEdOmega(phi, theta, t), dEdOmega_, Gamma0_, idx_);
}

/********************************************************************************************************************
 * METHOD: PowerLawJet::dEdOmega
 * DESCRIPTION: Returns the energy per unit solid angle for a PowerLawJet. For angles beyond the core (theta >=
 *theta_c_), the energy falls off as a power law with index k_.
 ********************************************************************************************************************/
double PowerLawJet::dEdOmega(double phi, double theta, double t) const {
    return dEdOmega_ * fastPow(theta / theta_c_, -k_);
}

/********************************************************************************************************************
 * METHOD: PowerLawJet::Gamma0
 * DESCRIPTION: Computes the initial Lorentz factor for a PowerLawJet using the LiangGhirlanda2010_ helper.
 *              The function uses the computed dEdOmega (from the power-law profile), along with the jet's
 *              dEdOmega_ and Gamma0_ parameters and index idx_.
 ********************************************************************************************************************/
double PowerLawJet::Gamma0(double phi, double theta, double t) const {
    return LiangGhirlanda2010_(dEdOmega(phi, theta, t), dEdOmega_, Gamma0_, idx_);
}