//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <cmath>

#include "macros.h"
#include "xtensor/containers/xadapt.hpp"
#include "xtensor/containers/xtensor.hpp"
#include "xtensor/core/xmath.hpp"
#include "xtensor/core/xnoalias.hpp"
#include "xtensor/views/xview.hpp"

// Type aliases for commonly used tensor types
using Array = xt::xtensor<Real, 1>;       // 1D array for storing 1D data (e.g., time points)
using MeshGrid = xt::xtensor<Real, 2>;    // 2D grid for storing 2D data (e.g., spatial coordinates)
using MeshGrid3d = xt::xtensor<Real, 3>;  // 3D grid for storing 3D data (e.g., full spatial-temporal data)
using MaskGrid = xt::xtensor<int, 3>;     // 3D boolean grid for masking operations

/********************************************************************************************************************
 * FUNCTION PROTOTYPES: Array and Grid Utilities
 * DESCRIPTION: Declares a set of functions for generating and processing arrays and grids. These functions include:
 *              - Converting boundary arrays to center arrays (linear and logarithmic),
 *              - Creating linearly and logarithmically spaced arrays,
 *              - Creating arrays with uniform spacing in cosine,
 *              - Generating arrays of zeros and ones,
 *              - Finding the minimum and maximum of grids,
 *              - Checking if an array is linearly or logarithmically scaled,
 *              - Creating 2D and 3D grids.
 ********************************************************************************************************************/

// Check if an array is linearly scaled within a given tolerance
bool isLinearScale(Array const& arr, Real tolerance = 1e-6);

// Check if an array is logarithmically scaled within a given tolerance
bool isLogScale(Array const& arr, Real tolerance = 1e-6);

/********************************************************************************************************************
 * CLASS: Coord
 * DESCRIPTION: Represents a coordinate system with arrays for phi, theta, and t. This class is used to define
 *              the computational grid for GRB simulations. It stores the angular coordinates (phi, theta) and
 *              time (t) arrays, along with derived quantities needed for numerical calculations.
 ********************************************************************************************************************/
class Coord {
   public:
    // Constructor taking phi, theta, and t arrays
    Coord(Array const& phi, Array const& theta, Array const& t);

    // Default constructor
    Coord() {};

    Array phi;    // Array of azimuthal angles (phi) in radians
    Array theta;  // Array of polar angles (theta) in radians
    Array t;      // Array of time points in seconds

    // Returns the dimensions of the coordinate arrays as a tuple (n_phi, n_theta, n_t)
    auto shape() const { return std::make_tuple(phi.size(), theta.size(), t.size()); }
};

/********************************************************************************************************************
 * TEMPLATE FUNCTION: logspace
 * DESCRIPTION: Creates an array with logarithmically spaced values between start and end.
 *              This is useful for creating time grids where we need more resolution at early times.
 * @param start: Starting value (log10)
 * @param end: Ending value (log10)
 * @param result: Output array to store the results
 ********************************************************************************************************************/
template <typename Arr>
void logspace(Real start, Real end, Arr& result) {
    size_t num = result.size();
    Real log_start = start;
    Real log_end = end;

    // Calculate step size in log space
    Real step = (log_end - log_start) / (num - 1);
    if (num == 1) {
        step = 0;  // Handle single-point array case
    }

    // Fill array with logarithmically spaced values
    for (size_t i = 0; i < num; i++) {
        result[i] = std::pow(10, log_start + i * step);
    }
}

/********************************************************************************************************************
 * TEMPLATE FUNCTION: boundaryToCenter (linear)
 * DESCRIPTION: Converts boundary values to center values using linear interpolation.
 *              This is used to compute cell-centered values from cell-boundary values.
 * @param boundary: Input array of boundary values
 * @param center: Output array of center values
 ********************************************************************************************************************/
template <typename Arr1, typename Arr2>
void boundaryToCenter(Arr1 const& boundary, Arr2& center) {
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = 0.5 * (boundary[i] + boundary[i + 1]);
    }
}

/********************************************************************************************************************
 * TEMPLATE FUNCTION: boundaryToCenterLog
 * DESCRIPTION: Converts boundary values to center values using geometric mean (logarithmic interpolation).
 *              This is used when working with logarithmically scaled quantities.
 * @param boundary: Input array of boundary values
 * @param center: Output array of center values
 ********************************************************************************************************************/
template <typename Arr1, typename Arr2>
void boundaryToCenterLog(Arr1 const& boundary, Arr2& center) {
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = std::sqrt(boundary[i] * boundary[i + 1]);
    }
}

// Non-template versions of the boundary-to-center conversion functions
Array boundaryToCenter(Array const& boundary);

Array boundaryToCenterLog(Array const& boundary);