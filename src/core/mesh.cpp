//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "mesh.h"

#include <cmath>
#include <iostream>

#include "macros.h"

/********************************************************************************************************************
 * FUNCTION: uniform_cos
 * DESCRIPTION: Creates and returns an Array of 'num' values uniformly spaced in cosine between angles 'start' and
 *'end'.
 ********************************************************************************************************************/
Array uniform_cos(Real start, Real end, size_t num) {
    // First generate a linearly spaced array in cosine space.
    Array result = xt::linspace(std::cos(start), std::cos(end), num);
    for (size_t i = 0; i < num; i++) {
        // Convert back to angle by taking the arccosine.
        result[i] = std::acos(result[i]);
    }
    return result;
}

/********************************************************************************************************************
 * FUNCTION: uniform_sin
 * DESCRIPTION: Creates and returns an Array of 'num' values uniformly spaced in sine between angles 'start' and
 *'end'.
 ********************************************************************************************************************/
Array uniform_sin(Real start, Real end, size_t num) {
    // First generate a linearly spaced array in sine space.
    Array result = xt::linspace(std::sin(start), std::sin(end), num);
    for (size_t i = 0; i < num; i++) {
        // Convert back to angle by taking the arcsine.
        result[i] = std::asin(result[i]);
    }
    return result;
}

/********************************************************************************************************************
 * FUNCTION: is_linear_scale
 * DESCRIPTION: Checks if the values in the given Array are approximately linearly spaced within the specified
 *tolerance.
 ********************************************************************************************************************/
bool is_linear_scale(Array const& arr, Real tolerance) {
    if (arr.size() < 2) return false;  // At least two elements are needed.

    Real diff = arr[1] - arr[0];
    for (size_t i = 2; i < arr.size(); ++i) {
        if (std::fabs((arr[i] - arr[i - 1] - diff) / diff) > tolerance) {
            return false;
        }
    }
    return true;
}

/********************************************************************************************************************
 * FUNCTION: is_log_scale
 * DESCRIPTION: Checks if the values in the given Array are approximately logarithmically spaced (constant ratio)
 *              within the specified tolerance.
 ********************************************************************************************************************/
bool is_log_scale(Array const& arr, Real tolerance) {
    if (arr.size() < 2) return false;  // At least two elements are needed.

    Real ratio = arr[1] / arr[0];
    for (size_t i = 2; i < arr.size(); ++i) {
        if (std::fabs((arr[i] / arr[i - 1] - ratio) / ratio) > tolerance) {
            return false;
        }
    }
    return true;
}

/********************************************************************************************************************
 * FUNCTION: boundary_to_center
 * DESCRIPTION: Converts a boundary array to center values by averaging adjacent boundaries.
 ********************************************************************************************************************/
Array boundary_to_center(Array const& boundary) {
    Array center({boundary.size() - 1}, 0);
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = 0.5 * (boundary[i] + boundary[i + 1]);
    }
    return center;
}

/********************************************************************************************************************
 * FUNCTION: boundary_to_center_log
 * DESCRIPTION: Converts a boundary array to center values in logarithmic space using the geometric mean.
 ********************************************************************************************************************/
Array boundary_to_center_log(Array const& boundary) {
    Array center({boundary.size() - 1}, 0);
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = std::sqrt(boundary[i] * boundary[i + 1]);
    }
    return center;
}