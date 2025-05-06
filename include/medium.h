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
 * @class Medium
 * @brief Represents the generic medium or any user-defined surrounding medium that the GRB jet interacts with.
 * @details The class provides methods to compute the density (rho) as a function of position (phi, theta, r).
 ********************************************************************************************************************/
class Medium {
   public:
    /// Density function that returns the mass density at a given position (phi, theta, r)
    /// The function is initialized to zero by default
    TernaryFunc rho{func::zero_3d};  // density(phi, theta, r)

    /// Mass function that returns the integrated mass within a radius r at position (phi, theta)
    TernaryFunc mass{func::zero_3d};  // mass(phi, theta, r)
};

/********************************************************************************************************************
 * @class ISM
 * @brief Implements a uniform interstellar medium (ISM) with constant density.
 * @details Provides methods to compute density and integrated mass at any position.
 *          The ISM is characterized by the particle number density n_ism.
 ********************************************************************************************************************/
class ISM {
   public:
    /********************************************************************************************************************
     * @brief Constructor: Initialize with particle number density in cm^-3
     * @param n_ism Particle number density in cm^-3
     ********************************************************************************************************************/
    ISM(Real n_ism) noexcept : rho_(n_ism * con::mp) {}

    /********************************************************************************************************************
     * @brief Return density at given position (constant everywhere)
     * @param phi Azimuthal angle (unused)
     * @param theta Polar angle (unused)
     * @param r Radial distance (unused)
     * @return Constant density value
     ********************************************************************************************************************/
    inline Real rho(Real /*phi*/, Real /*theta*/, Real /*r*/) const noexcept { return rho_; }

    /********************************************************************************************************************
     * @brief Return integrated mass within radius r (proportional to r^3)
     * @param phi Azimuthal angle (unused)
     * @param theta Polar angle (unused)
     * @param r Radial distance
     * @return Mass enclosed within radius r (= rho * r^3 / 3)
     ********************************************************************************************************************/
    inline Real mass(Real /*phi*/, Real /*theta*/, Real r) const noexcept { return rho_ * r * r * r / 3; }

   private:
    Real rho_{0};  ///< Mass density (particle number density × proton mass)
};

/********************************************************************************************************************
 * @class Wind
 * @brief Implements a stellar wind medium with density proportional to 1/r².
 * @details Provides methods to compute density and integrated mass at any position.
 *          The wind is characterized by the wind parameter A_star.
 ********************************************************************************************************************/
class Wind {
   public:
    /********************************************************************************************************************
     * @brief Constructor: Initialize with wind parameter A_star (in standard units)
     * @param A_star Wind density parameter in standard units
     ********************************************************************************************************************/
    Wind(Real A_star) noexcept : A(A_star * 5e11 * unit::g / unit::cm) {}

    /********************************************************************************************************************
     * @brief Return density at given position (proportional to 1/r²)
     * @param phi Azimuthal angle (unused)
     * @param theta Polar angle (unused)
     * @param r Radial distance
     * @return Density value at radius r (= A/r²)
     ********************************************************************************************************************/
    inline Real rho(Real /*phi*/, Real /*theta*/, Real r) const noexcept { return A / (r * r); }

    /********************************************************************************************************************
     * @brief Return integrated mass within radius r (proportional to r)
     * @param phi Azimuthal angle (unused)
     * @param theta Polar angle (unused)
     * @param r Radial distance
     * @return Mass enclosed within radius r (= A * r)
     ********************************************************************************************************************/
    inline Real mass(Real /*phi*/, Real /*theta*/, Real r) const noexcept { return A * r; }

   private:
    Real A{0};  ///< Wind density parameter in physical units
};

/********************************************************************************************************************
 * @namespace evn
 * @brief Provides functions to create different types of ambient medium profiles.
 * @details These functions return lambda functions that compute the density at any given position.
 ********************************************************************************************************************/
namespace evn {
    /********************************************************************************************************************
     * @brief Creates a uniform interstellar medium (ISM) profile
     * @param n_ism Number density of particles in cm^-3
     * @return Pair of functions for density and mass calculation
     ********************************************************************************************************************/
    inline auto ISM(Real n_ism) {
        Real rho = n_ism * con::mp;

        return std::make_pair([rho](Real /*phi*/, Real /*theta*/, Real /*r*/) noexcept { return rho; },
                              [rho](Real /*phi*/, Real /*theta*/, Real r) noexcept { return rho * r * r * r / 3; });
    };

    /********************************************************************************************************************
     * @brief Creates a stellar wind medium profile
     * @param A_star Wind parameter in standard units
     * @return Pair of functions for density and mass calculation
     * @details Converts A_star to proper units (A_star * 5e11 g/cm) and returns functions that compute
     *          density = A/r² and mass = A*r, representing a steady-state stellar wind where density
     *          falls off as 1/r²
     ********************************************************************************************************************/
    inline auto wind(Real A_star) {
        // Convert A_star to proper units: A_star * 5e11 g/cm
        Real A = A_star * 5e11 * unit::g / unit::cm;

        // Return a function that computes density = A/r^2
        // This represents a steady-state stellar wind where density falls off as 1/r^2
        return std::make_pair([A](Real /*phi*/, Real /*theta*/, Real r) noexcept { return A / (r * r); },
                              [A](Real /*phi*/, Real /*theta*/, Real r) noexcept { return A * r; });
    }
}  // namespace evn
