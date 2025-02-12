//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _MESHES_
#define _MESHES_

#include <boost/multi_array.hpp>
#include <cmath>
#include <functional>
#include <vector>

#include "macros.h"
/********************************************************************************************************************
 * CONDITIONAL TYPE DEFINITIONS
 * DESCRIPTION: If STATIC_ARRAY is defined, fixed sizes (PHI_SIZE, THETA_SIZE, R_SIZE) are used and specific
 *              array types (Array2D, Array3D) are employed. Otherwise, generic boost::multi_array types are used.
 ********************************************************************************************************************/

using Array = boost::multi_array<Real, 1>;
using MeshGrid = boost::multi_array<Real, 2>;
using MeshGrid3d = boost::multi_array<Real, 3>;

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
Array boundaryToCenter(Array const& boundary);
Array boundaryToCenterLog(Array const& boundary);
Array linspace(Real start, Real end, size_t num);
Array logspace(Real start, Real end, size_t num);
Array uniform_cos(Real start, Real end, size_t num);
Array uniform_sin(Real start, Real end, size_t num);
Array zeros(size_t num);
Array ones(size_t num);
Real min(MeshGrid const& grid);
Real max(MeshGrid const& grid);
bool isLinearScale(Array const& arr, Real tolerance = 1e-6);
bool isLogScale(Array const& arr, Real tolerance = 1e-6);
MeshGrid createGrid(size_t theta_size, size_t t_size, Real val = 0);
MeshGrid createGridLike(MeshGrid const& grid, Real val = 0);
MeshGrid3d create3DGrid(size_t phi_size, size_t theta_size, size_t t_size, Real val = 0);
MeshGrid3d create3DGridLike(MeshGrid3d const& grid, Real val = 0);

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
    Coord() = delete;

    Array phi;    // Array of phi values
    Array theta;  // Array of theta values
    Array t;      // Array of radius values

    Array dphi;  // Differential phi values
    Array dcos;  // Differential cosine of theta values

    // Returns the dimensions of the coordinate arrays.
    auto shape() const { return std::make_tuple(phi.size(), theta.size(), t.size()); }
};

/********************************************************************************************************************
 * TEMPLATE FUNCTIONS: min and max for MeshGrid types
 * DESCRIPTION: Provides variadic template functions to compute the minimum and maximum values across one or more grids.
 ********************************************************************************************************************/
template <typename... MeshGrid>
Real min(MeshGrid const&... grids) {
    return std::min(std::min(grids...));
}

template <typename... MeshGrid>
Real max(MeshGrid const&... grids) {
    return std::max(std::max(grids...));
}

/********************************************************************************************************************
 * TEMPLATE FUNCTION: linspace
 * DESCRIPTION: Fills the provided array 'result' with 'num' linearly spaced values between 'start' and 'end'.
 ********************************************************************************************************************/
template <typename Arr>
void linspace(Real start, Real end, Arr& result) {
    size_t num = result.size();
    // Handle empty array case.
    Real step = (end - start) / (num - 1);
    if (num == 1) {
        step = 0;
    }
    for (size_t i = 0; i < num; i++) {
        result[i] = start + i * step;
    }
}

/********************************************************************************************************************
 * TEMPLATE FUNCTION: logspace
 * DESCRIPTION: Fills the provided array 'result' with 'num' logarithmically spaced values between 'start' and 'end'.
 ********************************************************************************************************************/
template <typename Arr>
void logspace(Real start, Real end, Arr& result) {
    size_t num = result.size();
    Real log_start = std::log(start);
    Real log_end = std::log(end);

    Real step = (log_end - log_start) / (num - 1);
    if (num == 1) {
        step = 0;
    }  // Handle empty array case.
    for (size_t i = 0; i < num; i++) {
        result[i] = std::exp(log_start + i * step);
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
#endif