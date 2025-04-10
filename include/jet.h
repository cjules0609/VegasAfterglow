//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include "macros.h"
#include "mesh.h"
#include "utilities.h"

/********************************************************************************************************************
 * CLASS: Ejecta
 * DESCRIPTION: Represents ejecta properties for a simulation. It uses ternary functions (TernaryFunc) to
 *              describe various quantities as functions of phi, theta, and time.
 ********************************************************************************************************************/
class Ejecta {
   public:
    BinaryFunc dE0dOmega{func::zero0};   // Initial energy per unit solid angle
    BinaryFunc sigma0{func::zero0};      // Initial magnetization parameter
    BinaryFunc Gamma0{func::one0};       // Lorentz factor profile in the ejecta (default: one)
    TernaryFunc dEdtdOmega{func::zero};  // Energy injection rate per solid angle (default: zero)
    TernaryFunc dMdtdOmega{func::zero};  // Mass injection rate per unit solid angle (default: zero)

    Real T0{1 * con::sec};  // Duration of the ejecta
    bool spreading{false};  // Flag indicating if the ejecta spreads
    // (Additional member functions could be defined as needed.)
};

/********************************************************************************************************************
 * NAMESPACE: math
 * DESCRIPTION: Provides a collection of mathematical helper functions for combining functions,
 *              constructing various injection profiles (e.g., tophat, Gaussian, power-law), and computing integrals.
 ********************************************************************************************************************/
namespace math {
    // Combines a spatial function and a temporal function into one function of (phi, theta, t).
    template <typename F1, typename F2>
    inline auto combine(F1 f_spatial, F2 f_temporal) {
        return [=](Real phi, Real theta, Real t) -> Real { return f_spatial(phi, theta) * f_temporal(t); };
    }

    template <typename F1>
    inline auto t_indep(F1 f_spatial) {
        return [=](Real phi, Real theta, Real t) -> Real { return f_spatial(phi, theta); };
    }

    // Returns a constant (isotropic) function with a fixed height.
    inline auto isotropic(Real height) {
        return [=](Real phi, Real theta) -> Real { return height; };
    }

    // Returns a tophat function: it is constant (height) within the core angle theta_c and zero outside.
    inline auto tophat(Real theta_c, Real hight) {
        return [=](Real phi, Real theta) -> Real { return theta < theta_c ? hight : 0; };
    }

    // Returns a Gaussian profile function for jet properties.
    inline auto gaussian(Real theta_c, Real height) {
        return [=](Real phi, Real theta) -> Real { return height * fastExp(-theta * theta / (2 * theta_c * theta_c)); };
    }

    // Returns a power-law profile function for jet properties.
    inline auto powerLaw(Real theta_c, Real height, Real k) {
        return [=](Real phi, Real theta) -> Real { return height * fastPow(1 + theta / theta_c, -k); };
    }

    // Constant injection: returns 1 regardless of time.
    inline auto constInjection() {
        return [=](Real t) -> Real { return 1; };
    }

    // Step injection: returns 1 if t is greater than t0, else 0.
    inline auto stepInjection(Real t0) {
        return [=](Real t) -> Real { return t > t0 ? 1 : 0; };
    }

    // Square injection: returns 1 if t is between t0 and t1, else 0.
    inline auto squareInjection(Real t0, Real t1) {
        return [=](Real t) -> Real { return t > t0 && t < t1 ? 1 : 0; };
    }

    // Power-law injection: returns a decaying function with power-law index q.
    inline auto powerLawInjection(Real t0, Real q) {
        return [=](Real t) -> Real { return fastPow(1 + t / t0, -q); };
    }
}  // namespace math

/********************************************************************************************************************
 * FUNCTION: LiangGhirlanda2010
 * DESCRIPTION: Returns a lambda function that computes a Lorentz factor using the Liang & Ghirlanda (2010)
 *              prescription. The lambda takes an energy function (energy_func) and parameters e_max, gamma_max,
 *              and idx, and returns a function of (phi, theta, t). The returned function calculates:
 *                  e = energy_func(phi, theta, 0)
 *                  u = (e / e_max)^idx * gamma_max
 *                  return sqrt(1 + u^2)
 ********************************************************************************************************************/
template <typename F>
auto LiangGhirlanda2010(F energy_func, Real e_max, Real gamma_max, Real idx) {
    return [=](Real phi, Real theta, Real t = 0) -> Real {
        Real e = energy_func(phi, theta);
        Real u = fastPow(e / e_max, idx) * gamma_max;
        return std::sqrt(1 + u * u);
    };
}
