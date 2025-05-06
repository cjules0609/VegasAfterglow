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
 * @brief Creates and returns an Array of 'num' values uniformly spaced in cosine between angles 'start' and 'end'.
 * @details First generates a linearly spaced array in cosine space, then converts back to angles by taking
 *          the arccosine.
 * @param start The starting angle in radians
 * @param end The ending angle in radians
 * @param num The number of values to generate
 * @return Array of angles with uniform cosine spacing
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
 * @brief Creates and returns an Array of 'num' values uniformly spaced in sine between angles 'start' and 'end'.
 * @details First generates a linearly spaced array in sine space, then converts back to angles by taking
 *          the arcsine.
 * @param start The starting angle in radians
 * @param end The ending angle in radians
 * @param num The number of values to generate
 * @return Array of angles with uniform sine spacing
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
 * @brief Checks if the values in the given Array are approximately linearly spaced within the specified tolerance.
 * @param arr The array to check
 * @param tolerance Maximum allowed deviation from linearity
 * @return True if array elements are linearly spaced, false otherwise
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
 * @brief Checks if the values in the given Array are approximately logarithmically spaced within the specified
 *tolerance.
 * @details Tests if the ratio between consecutive elements is approximately constant.
 * @param arr The array to check
 * @param tolerance Maximum allowed deviation from constant ratio
 * @return True if array elements are logarithmically spaced, false otherwise
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
 * @brief Converts a boundary array to center values by averaging adjacent boundaries.
 * @param boundary Array of boundary values
 * @return Array of center values computed as the average of adjacent boundaries
 ********************************************************************************************************************/
Array boundary_to_center(Array const& boundary) {
    Array center({boundary.size() - 1}, 0);
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = 0.5 * (boundary[i] + boundary[i + 1]);
    }
    return center;
}

/********************************************************************************************************************
 * @brief Converts a boundary array to center values in logarithmic space using the geometric mean.
 * @param boundary Array of boundary values
 * @return Array of center values computed as the geometric mean of adjacent boundaries
 ********************************************************************************************************************/
Array boundary_to_center_log(Array const& boundary) {
    Array center({boundary.size() - 1}, 0);
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = std::sqrt(boundary[i] * boundary[i + 1]);
    }
    return center;
}