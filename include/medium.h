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
    TernaryFunc rho{func::zero_3d};   // density(phi, theta, r)
    TernaryFunc mass{func::zero_3d};  // mass(phi, theta, r)
};

/********************************************************************************************************************
 * CLASS: ISM
 * DESCRIPTION: Implements a uniform interstellar medium (ISM) with constant density.
 *              Provides methods to compute density and integrated mass at any position.
 *              The ISM is characterized by the particle number density n_ism.
 ********************************************************************************************************************/
class ISM {
   public:
    // Constructor: Initialize with particle number density in cm^-3
    constexpr ISM(Real n_ism) noexcept : rho_(n_ism * con::mp) {}

    // Return density at given position (constant everywhere)
    constexpr inline Real rho(Real /*phi*/, Real /*theta*/, Real /*r*/) const noexcept { return rho_; }

    // Return integrated mass within radius r (proportional to r^3)
    constexpr inline Real mass(Real /*phi*/, Real /*theta*/, Real r) const noexcept { return rho_ * r * r * r / 3; }

   private:
    Real rho_{0};  // Mass density (particle number density × proton mass)
};

/********************************************************************************************************************
 * CLASS: Wind
 * DESCRIPTION: Implements a stellar wind medium with density proportional to 1/r².
 *              Provides methods to compute density and integrated mass at any position.
 *              The wind is characterized by the wind parameter A_star.
 ********************************************************************************************************************/
class Wind {
   public:
    // Constructor: Initialize with wind parameter A_star (in standard units)
    constexpr Wind(Real A_star) noexcept : A(A_star * 5e11 * unit::g / unit::cm) {}

    // Return density at given position (proportional to 1/r²)
    constexpr inline Real rho(Real /*phi*/, Real /*theta*/, Real r) const noexcept { return A / (r * r); }

    // Return integrated mass within radius r (proportional to r)
    constexpr inline Real mass(Real /*phi*/, Real /*theta*/, Real r) const noexcept { return A * r; }

   private:
    Real A{0};  // Wind density parameter in physical units
};

/********************************************************************************************************************
 * NAMESPACE: evn
 * DESCRIPTION: Provides functions to create different types of ambient medium profiles.
 *              These functions return lambda functions that compute the density at any given position.
 ********************************************************************************************************************/
namespace evn {
    // Creates a uniform interstellar medium (ISM) profile
    // Parameter n_ism: Number density of particles in cm^-3
    inline auto ISM(Real n_ism) {
        Real rho = n_ism * con::mp;

        return std::make_pair([rho](Real /*phi*/, Real /*theta*/, Real /*r*/) noexcept { return rho; },
                              [rho](Real /*phi*/, Real /*theta*/, Real r) noexcept { return rho * r * r * r / 3; });
    };

    // Creates a stellar wind medium profile
    // Parameter A_star: Wind parameter in standard units
    inline auto wind(Real A_star) {
        // Convert A_star to proper units: A_star * 5e11 g/cm
        Real A = A_star * 5e11 * unit::g / unit::cm;

        // Return a function that computes density = A/r^2
        // This represents a steady-state stellar wind where density falls off as 1/r^2
        return std::make_pair([A](Real /*phi*/, Real /*theta*/, Real r) noexcept { return A / (r * r); },
                              [A](Real /*phi*/, Real /*theta*/, Real r) noexcept { return A * r; });
    }
}  // namespace evn
