//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "medium.h"

#include <cmath>

#include "macros.h"

/********************************************************************************************************************
 * CONSTRUCTOR: Medium::Medium(Real k, Real n_c, Real r_c)
 * DESCRIPTION: Initializes a Medium object with a specified density profile index (k), a characteristic number
 *              density (n_c), and a characteristic radius (r_c).
 ********************************************************************************************************************/
Medium::Medium(Real k, Real n_c, Real r_c) : k(k), n_c(n_c), r_c(r_c) {};

/********************************************************************************************************************
 * CONSTRUCTOR: Medium::Medium(Real n_c)
 * DESCRIPTION: Initializes a Medium object with constant density (k = 0). The characteristic density is set
 *              to n_c, and the characteristic radius is set to con::cm (a constant representing 1 cm).
 ********************************************************************************************************************/
Medium::Medium(Real n_c) : k(0), n_c(n_c), r_c(con::cm) {};

/********************************************************************************************************************
 * FUNCTION: createISM
 * DESCRIPTION: Factory function that creates a Medium object representing the interstellar medium (ISM)
 *              using the provided number density (n_ism). It uses the constructor for constant density.
 ********************************************************************************************************************/
Medium createISM(Real n_ism) { return Medium(n_ism); }
