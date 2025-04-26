//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include "macros.h"
#include "utilities.h"
/********************************************************************************************************************
 * CLASS: Medium
 * DESCRIPTION: Represents the generic medium or any user-defined surrounding medium that the GRB jet interacts with.
 *              The class provides methods to compute the density (rho) as a function of position (phi, theta, r).
 ********************************************************************************************************************/
class Medium {
   public:
    // Density function that returns the mass density at a given position (phi, theta, r)
    // The function is initialized to zero by default
    TernaryFunc rho{func::zero_3d};  // density(phi, theta, r)
};

class ISM {
   public:
    constexpr ISM(Real n_ism) noexcept : rho_(n_ism * con::mp) {}

    constexpr inline Real rho(Real phi, Real theta, Real r) const noexcept { return rho_; }

    constexpr inline Real mass(Real phi, Real theta, Real r) const noexcept { return rho_ * r * r * r / 3; }

   private:
    Real rho_{0};
};

class Wind {
   public:
    constexpr Wind(Real A_star) noexcept : A(A_star * 5e11 * con::g / con::cm) {}

    constexpr inline Real rho(Real phi, Real theta, Real r) const noexcept { return A / (r * r); }

    constexpr inline Real mass(Real phi, Real theta, Real r) const noexcept { return A * r; }

   private:
    Real A{0};
};

/********************************************************************************************************************
 * NAMESPACE: evn
 * DESCRIPTION: Provides functions to create different types of ambient medium profiles.
 *              These functions return lambda functions that compute the density at any given position.
 ********************************************************************************************************************/
namespace evn {
    // Creates a uniform interstellar medium (ISM) profile
    inline auto ISM(Real n_ism) {
        Real rho = n_ism * con::mp;

        return [rho](Real phi, Real theta, Real r) noexcept {
            return rho;  // Convert number density to mass density
        };
    }

    // Creates a stellar wind medium profile
    inline auto wind(Real A_star) {
        // Convert A_star to proper units: A_star * 5e11 g/cm
        Real A = A_star * 5e11 * con::g / con::cm;

        // Return a function that computes density = A/r^2
        // This represents a steady-state stellar wind where density falls off as 1/r^2
        return [A](Real phi, Real theta, Real r) noexcept { return A / (r * r); };
    }
}  // namespace evn
