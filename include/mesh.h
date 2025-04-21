//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <cmath>
#include <functional>
#include <vector>

#include "macros.h"
#include "xtensor/containers/xadapt.hpp"
#include "xtensor/containers/xtensor.hpp"
#include "xtensor/core/xmath.hpp"
#include "xtensor/views/xview.hpp"

using Array = xt::xtensor<Real, 1>;
using MeshGrid = xt::xtensor<Real, 2>;
using MeshGrid3d = xt::xtensor<Real, 3>;
using MaskGrid = xt::xtensor<bool, 3>;
/********************************************************************************************************************
 * FUNCTION TYPE DEFINITIONS
 * DESCRIPTION: Defines convenient aliases for unary, binary, and ternary functions operating on Reals.
 ********************************************************************************************************************/
using UnaryFunc = std::function<Real(Real)>;
using BinaryFunc = std::function<Real(Real, Real)>;
using TernaryFunc = std::function<Real(Real, Real, Real)>;

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

bool isLinearScale(Array const& arr, Real tolerance = 1e-6);
bool isLogScale(Array const& arr, Real tolerance = 1e-6);

/********************************************************************************************************************
 * CLASS: Coord
 * DESCRIPTION: Represents a coordinate system with arrays for phi, theta, and t. Also stores derived quantities
 *              such as differential phi (dphi) and differential cosine (dcos) for theta. The class also holds
 *              the minimum and maximum observation times.
 ********************************************************************************************************************/
class Coord {
   public:
    // Constructor taking phi, theta, and t arrays. Default constructor is deleted.
    Coord(Array const& phi, Array const& theta, Array const& t);
    Coord() {};
    // Coord() = delete;

    Array phi;    // Array of phi values
    Array theta;  // Array of theta values
    Array t;      // Array of radius values

    // Returns the dimensions of the coordinate arrays.
    auto shape() const { return std::make_tuple(phi.size(), theta.size(), t.size()); }
};

template <typename Arr>
void logspace(Real start, Real end, Arr& result) {
    size_t num = result.size();
    Real log_start = start;
    Real log_end = end;

    Real step = (log_end - log_start) / (num - 1);
    if (num == 1) {
        step = 0;
    }  // Handle empty array case.
    for (size_t i = 0; i < num; i++) {
        result[i] = std::pow(10, log_start + i * step);
    }
}
/********************************************************************************************************************
 * TEMPLATE FUNCTION: boundaryToCenter (linear)
 * DESCRIPTION: Given a boundary array, computes the center values by averaging adjacent boundaries.
 ********************************************************************************************************************/
template <typename Arr1, typename Arr2>
void boundaryToCenter(Arr1 const& boundary, Arr2& center) {
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = 0.5 * (boundary[i] + boundary[i + 1]);
    }
}

/********************************************************************************************************************
 * TEMPLATE FUNCTION: boundaryToCenterLog
 * DESCRIPTION: Given a boundary array, computes the center values in logarithmic space using the geometric mean.
 ********************************************************************************************************************/
template <typename Arr1, typename Arr2>
void boundaryToCenterLog(Arr1 const& boundary, Arr2& center) {
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = std::sqrt(boundary[i] * boundary[i + 1]);
    }
}

Array boundaryToCenter(Array const& boundary);
Array boundaryToCenterLog(Array const& boundary);