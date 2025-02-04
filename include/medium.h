//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _MEDIUM_
#define _MEDIUM_

#include <iostream>

#include "macros.h"
#include "mesh.h"
#include "utilities.h"

/********************************************************************************************************************
 * CLASS: Medium
 * DESCRIPTION: Represents the interstellar medium (ISM) or any surrounding medium.
 *              The medium is characterized by a density profile with index k, a characteristic number density n_c,
 *              and a characteristic radius r_c. The class provides methods to compute the mass density
 *              (rho) and the enclosed mass up to radius r.
 ********************************************************************************************************************/
class Medium {
   public:
    // Constructor: Initializes a Medium object with density profile index k, characteristic density n_c, and radius
    // r_c.
    Medium(Real k, Real n_c, Real r_c);
    // Constructor: Initializes a Medium object with only characteristic density n_c. (k and r_c may be set to default
    // values.)
    Medium(Real n_c);

    // Returns the mass density  at radius phi, theta r.
    Real rho(Real r) const;
    // Returns the total mass enclosed within radius r.
    Real mass(Real r) const;

   private:
    Real const k;    // Density profile index (if 0, constant density)
    Real const n_c;  // Characteristic number density (in cm^-3)
    Real const r_c;  // Characteristic radius (in cm)
};

/********************************************************************************************************************
 * INLINE METHOD: Medium::rho
 * DESCRIPTION: Computes the mass density at a given radius r.
 *              - For a constant density profile (k == 0), returns n_c * con::mp.
 *              - Otherwise, scales the density as (r / r_c)^(-k) times n_c * con::mp.
 ********************************************************************************************************************/
inline Real Medium::rho(Real r) const {
    if (k == 0) {
        return n_c * con::mp;
    } else {
        return n_c * fastPow(r / r_c, -k) * con::mp;
    }
}

/********************************************************************************************************************
 * INLINE METHOD: Medium::mass
 * DESCRIPTION: Computes the enclosed mass per solid angle within radius r.
 *              - For a constant density (k == 0): Uses the standard spherical mass formula.
 *              - For k == 3: Returns a mass proportional to the logarithm of r.
 *              - For other k values: Uses the general power-law mass formula.
 ********************************************************************************************************************/
inline Real Medium::mass(Real r) const {
    if (k == 0) {
        // Mass for constant density: (4π/3) * r^3 * n_c * mp.
        return (4 * con::pi * con::mp / 3) * r * r * r * n_c;
    } else if (k == 3) {
        // For k = 3, the enclosed mass grows logarithmically with radius.
        return (4 * con::pi * con::mp) * n_c * r_c * r_c * r_c * fastLog(r / con::cm);
    } else {
        // For other values of k, use the general formula: 4π/(3-k) * r^3 * n_c * (r/r_c)^(-k) * mp.
        return (4 * con::pi * con::mp) / (3 - k) * r * r * r * n_c * fastPow(r / r_c, -k);
    }
}

/********************************************************************************************************************
 * FUNCTION: createISM
 * DESCRIPTION: Factory function to create a Medium object representing the interstellar medium (ISM)
 *              with the given number density n_ism. (Additional parameters may be set to default values internally.)
 ********************************************************************************************************************/
Medium createISM(Real n_ism);

#endif