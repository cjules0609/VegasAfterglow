//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <iostream>

#include "macros.h"
#include "mesh.h"
#include "utilities.h"

/********************************************************************************************************************
 * CLASS: Medium
 * DESCRIPTION: Represents the interstellar medium (ISM) or any surrounding medium that the GRB jet interacts with.
 *              The class provides methods to compute the density (rho) as a function of position (phi, theta, r).
 *              This is crucial for modeling the shock dynamics and energy dissipation in the afterglow.
 ********************************************************************************************************************/
class Medium {
   public:
    // Density function that returns the mass density at a given position (phi, theta, r)
    // The function is initialized to zero by default
    TernaryFunc rho{func::zero};  // density(phi, theta, r)
};

/********************************************************************************************************************
 * NAMESPACE: evn
 * DESCRIPTION: Provides functions to create different types of ambient medium profiles.
 *              These functions return lambda functions that compute the density at any given position.
 *              Currently supports two types of media: uniform ISM and stellar wind.
 ********************************************************************************************************************/
namespace evn {
    // Creates a uniform interstellar medium (ISM) profile
    inline auto ISM(Real n_ism) {
        return [n_ism](Real phi, Real theta, Real r) {
            return n_ism * con::mp;  // Convert number density to mass density
        };
    }

    // Creates a stellar wind medium profile
    inline auto wind(Real A_star) {
        // Convert A_star to proper units: A_star * 5e11 g/cm
        Real A = A_star * 5e11 * con::g / con::cm;

        // Return a function that computes density = A/r^2
        // This represents a steady-state stellar wind where density falls off as 1/r^2
        return [A](Real phi, Real theta, Real r) { return A / (r * r); };
    }
}  // namespace evn
