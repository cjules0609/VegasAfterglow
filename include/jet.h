//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <cmath>

#include "utilities.h"
/********************************************************************************************************************
 * CLASS: Ejecta
 * DESCRIPTION: Represents generic ejecta properties for a simulation. It uses ternary functions (TernaryFunc) to
 *              accept user-defined ejecta that describes various quantities as functions of phi, theta, and time. This
 *              class encapsulates all the properties of the material ejected in a gamma-ray burst, including its
 *              energy, magnetization, and Lorentz factor profiles.
 ********************************************************************************************************************/
class Ejecta {
   public:
    // Initial energy per unit solid angle as a function of (phi, theta)
    BinaryFunc dE0dOmega{func::zero_2d};

    // Lorentz factor profile in the ejecta as a function of (phi, theta)
    // Default is uniform (one) across all angles
    BinaryFunc Gamma0{func::one_2d};

    // Initial magnetization parameter as a function of (phi, theta)
    BinaryFunc sigma0{func::zero_2d};

    // Energy injection rate per solid angle as a function of (phi, theta, t)
    // Default is no energy injection (zero)
    TernaryFunc dEdtdOmega{func::zero_3d};

    // Mass injection rate per unit solid angle as a function of (phi, theta, t)
    // Default is no mass injection (zero)
    TernaryFunc dMdtdOmega{func::zero_3d};

    // Duration of the ejecta in seconds
    Real T0{1 * con::sec};

    // Flag indicating if the ejecta spreads laterally during evolution
    bool spreading{false};
};

class TophatJet {
   public:
    constexpr TophatJet(Real theta_c, Real E_iso, Real Gamma0) noexcept
        : theta_c_(theta_c), eps_iso_(E_iso / (4 * con::pi)), Gamma0_(Gamma0) {}

    constexpr inline Real dE0dOmega(Real phi, Real theta) const noexcept { return theta < theta_c_ ? eps_iso_ : 0; }

    constexpr inline Real Gamma0(Real phi, Real theta) const noexcept { return theta < theta_c_ ? Gamma0_ : 1; }

    Real T0{1 * con::sec};
    bool spreading{false};

   private:
    Real const theta_c_{0};
    Real const eps_iso_{0};
    Real const Gamma0_{1};
};

class GaussianJet {
   public:
    constexpr GaussianJet(Real theta_c, Real E_iso, Real Gamma0) noexcept
        : norm_(-1 / (2 * theta_c * theta_c)), eps_iso_(E_iso / (4 * con::pi)), Gamma0_(Gamma0) {}

    constexpr inline Real dE0dOmega(Real phi, Real theta) const noexcept {
        return eps_iso_ * fastExp(theta * theta * norm_);
    }

    constexpr inline Real Gamma0(Real phi, Real theta) const noexcept {
        return Gamma0_ * fastExp(theta * theta * norm_);
    }

    Real T0{1 * con::sec};
    bool spreading{false};

   private:
    Real const norm_{0};
    Real const eps_iso_{0};
    Real const Gamma0_{1};
};

class PowerLawJet {
   public:
    constexpr PowerLawJet(Real theta_c, Real E_iso, Real Gamma0, Real k) noexcept
        : theta_c_(theta_c), eps_iso_(E_iso / (4 * con::pi)), Gamma0_(Gamma0), k_(k) {}

    constexpr inline Real dE0dOmega(Real phi, Real theta) const noexcept {
        return eps_iso_ * fastPow(theta / theta_c_, -k_);
    }

    constexpr inline Real Gamma0(Real phi, Real theta) const noexcept {
        return Gamma0_ * fastPow(theta / theta_c_, -k_);
    }

    Real T0{1 * con::sec};
    bool spreading{false};

   private:
    Real const theta_c_{0};
    Real const eps_iso_{0};
    Real const Gamma0_{1};
    Real const k_{2};
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
    constexpr inline auto combine(F1 f_spatial, F2 f_temporal) noexcept {
        return [=](Real phi, Real theta, Real t) constexpr noexcept { return f_spatial(phi, theta) * f_temporal(t); };
    }

    // Creates a time-independent function from a spatial function
    // by ignoring the time parameter t
    template <typename F1>
    constexpr inline auto t_indep(F1 f_spatial) noexcept {
        return [=](Real phi, Real theta, Real t) constexpr noexcept { return f_spatial(phi, theta); };
    }

    // Returns a constant (isotropic) function with a fixed height.
    // This creates a uniform distribution across all angles.
    constexpr inline auto isotropic(Real height) noexcept {
        return [=](Real phi, Real theta) constexpr noexcept { return height; };
    }

    // Returns a tophat function: it is constant (height) within the core angle theta_c
    // and zero outside. This creates a uniform jet with sharp edges.
    constexpr inline auto tophat(Real theta_c, Real hight) noexcept {
        return [=](Real phi, Real theta) constexpr noexcept { return theta < theta_c ? hight : 0; };
    }

    // Returns a Gaussian profile function for jet properties.
    // The profile peaks at theta = 0 and falls off exponentially with angle.
    constexpr inline auto gaussian(Real theta_c, Real height) noexcept {
        Real spread = 2 * theta_c * theta_c;
        return [=](Real phi, Real theta) constexpr noexcept { return height * fastExp(-theta * theta / spread); };
    }

    // Returns a power-law profile function for jet properties.
    // The profile follows a power-law decay with angle, controlled by the index k.
    constexpr inline auto powerLaw(Real theta_c, Real height, Real k) noexcept {
        return [=](Real phi, Real theta) constexpr noexcept { return height / (1 + fastPow(theta / theta_c, k)); };
    }

    // Creates a constant injection profile: returns 1 regardless of time.
    // This represents continuous energy/mass injection.
    constexpr inline auto constInjection() noexcept {
        return [=](Real t) constexpr noexcept { return 1.; };
    }

    // Creates a step injection profile: returns 1 if t > t0, else 0.
    // This represents a sudden turn-on of injection at time t0.
    constexpr inline auto stepInjection(Real t0) noexcept {
        return [=](Real t) constexpr noexcept { return t > t0 ? 1. : 0.; };
    }

    // Creates a square injection profile: returns 1 if t is between t0 and t1, else 0.
    // This represents a finite duration of injection between t0 and t1.
    constexpr inline auto squareInjection(Real t0, Real t1) noexcept {
        return [=](Real t) constexpr noexcept { return t > t0 && t < t1 ? 1. : 0.; };
    }

    // Creates a power-law injection profile: returns a decaying function with power-law index q.
    // The injection rate decays as (1 + t/t0)^(-q).
    constexpr inline auto powerLawInjection(Real t0, Real q) noexcept {
        return [=](Real t) constexpr noexcept { return fastPow(1 + t / t0, -q); };
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
    return [=](Real phi, Real theta, Real t = 0) noexcept {
        // Get the energy at the given angle
        Real e = energy_func(phi, theta);

        // Calculate the velocity parameter u using power-law scaling
        Real u = fastPow(e / e_max, idx) * gamma_max;

        // Convert to Lorentz factor
        return std::sqrt(1 + u * u);
    };
}
