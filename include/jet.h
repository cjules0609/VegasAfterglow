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
 *              describe various quantities as functions of phi, theta, and time. This class encapsulates
 *              all the properties of the material ejected in a gamma-ray burst, including its energy,
 *              magnetization, and Lorentz factor profiles.
 ********************************************************************************************************************/
class Ejecta {
   public:
    // Initial energy per unit solid angle as a function of (phi, theta)
    BinaryFunc dE0dOmega{func::zero0};

    // Initial magnetization parameter as a function of (phi, theta)
    BinaryFunc sigma0{func::zero0};

    // Lorentz factor profile in the ejecta as a function of (phi, theta)
    // Default is uniform (one) across all angles
    BinaryFunc Gamma0{func::one0};

    // Energy injection rate per solid angle as a function of (phi, theta, t)
    // Default is no energy injection (zero)
    TernaryFunc dEdtdOmega{func::zero};

    // Mass injection rate per unit solid angle as a function of (phi, theta, t)
    // Default is no mass injection (zero)
    TernaryFunc dMdtdOmega{func::zero};

    // Duration of the ejecta in seconds
    Real T0{1 * con::sec};

    // Flag indicating if the ejecta spreads laterally during evolution
    bool spreading{false};
};

/********************************************************************************************************************
 * NAMESPACE: math
 * DESCRIPTION: Provides a collection of mathematical helper functions for combining functions,
 *              constructing various injection profiles (e.g., tophat, Gaussian, power-law), and computing integrals.
 *              These functions are used to create complex jet profiles by combining spatial and temporal dependencies.
 ********************************************************************************************************************/
namespace math {
    // Combines a spatial function (f_spatial) and a temporal function (f_temporal)
    // into one function of (phi, theta, t). The result is the product of both functions.
    template <typename F1, typename F2>
    inline auto combine(F1 f_spatial, F2 f_temporal) {
        return [=](Real phi, Real theta, Real t) -> Real { return f_spatial(phi, theta) * f_temporal(t); };
    }

    // Creates a time-independent function from a spatial function
    // by ignoring the time parameter t
    template <typename F1>
    inline auto t_indep(F1 f_spatial) {
        return [=](Real phi, Real theta, Real t) -> Real { return f_spatial(phi, theta); };
    }

    // Returns a constant (isotropic) function with a fixed height.
    // This creates a uniform distribution across all angles.
    inline auto isotropic(Real height) {
        return [=](Real phi, Real theta) -> Real { return height; };
    }

    // Returns a tophat function: it is constant (height) within the core angle theta_c
    // and zero outside. This creates a uniform jet with sharp edges.
    inline auto tophat(Real theta_c, Real hight) {
        return [=](Real phi, Real theta) -> Real { return theta < theta_c ? hight : 0; };
    }

    // Returns a Gaussian profile function for jet properties.
    // The profile peaks at theta = 0 and falls off exponentially with angle.
    inline auto gaussian(Real theta_c, Real height) {
        return [=](Real phi, Real theta) -> Real { return height * fastExp(-theta * theta / (2 * theta_c * theta_c)); };
    }

    // Returns a power-law profile function for jet properties.
    // The profile follows a power-law decay with angle, controlled by the index k.
    inline auto powerLaw(Real theta_c, Real height, Real k) {
        return [=](Real phi, Real theta) -> Real { return height / (1 + fastPow(theta / theta_c, k)); };
    }

    // Creates a constant injection profile: returns 1 regardless of time.
    // This represents continuous energy/mass injection.
    inline auto constInjection() {
        return [=](Real t) -> Real { return 1; };
    }

    // Creates a step injection profile: returns 1 if t > t0, else 0.
    // This represents a sudden turn-on of injection at time t0.
    inline auto stepInjection(Real t0) {
        return [=](Real t) -> Real { return t > t0 ? 1 : 0; };
    }

    // Creates a square injection profile: returns 1 if t is between t0 and t1, else 0.
    // This represents a finite duration of injection between t0 and t1.
    inline auto squareInjection(Real t0, Real t1) {
        return [=](Real t) -> Real { return t > t0 && t < t1 ? 1 : 0; };
    }

    // Creates a power-law injection profile: returns a decaying function with power-law index q.
    // The injection rate decays as (1 + t/t0)^(-q).
    inline auto powerLawInjection(Real t0, Real q) {
        return [=](Real t) -> Real { return fastPow(1 + t / t0, -q); };
    }
}  // namespace math

/********************************************************************************************************************
 * FUNCTION: LiangGhirlanda2010
 * DESCRIPTION: Returns a lambda function that computes a Lorentz factor using the Liang & Ghirlanda (2010)
 *              prescription. This is a commonly used model for the Lorentz factor distribution in GRB jets.
 *              The lambda takes an energy function (energy_func) and parameters e_max, gamma_max, and idx,
 *              and returns a function of (phi, theta, t). The returned function calculates:
 *                  e = energy_func(phi, theta, 0)
 *                  u = (e / e_max)^idx * gamma_max
 *                  return sqrt(1 + u^2)
 *              This creates a power-law relationship between the energy and Lorentz factor of the jet.
 ********************************************************************************************************************/
template <typename F>
auto LiangGhirlanda2010(F energy_func, Real e_max, Real gamma_max, Real idx) {
    return [=](Real phi, Real theta, Real t = 0) -> Real {
        // Get the energy at the given angle
        Real e = energy_func(phi, theta);

        // Calculate the velocity parameter u using power-law scaling
        Real u = fastPow(e / e_max, idx) * gamma_max;

        // Convert to Lorentz factor
        return std::sqrt(1 + u * u);
    };
}
