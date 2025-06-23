//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <cmath>

#include "macros.h"
#include "physics.h"
#include "xtensor/containers/xadapt.hpp"
#include "xtensor/containers/xfixed.hpp"
#include "xtensor/containers/xtensor.hpp"
#include "xtensor/core/xmath.hpp"
#include "xtensor/core/xnoalias.hpp"
#include "xtensor/io/xnpy.hpp"
#include "xtensor/views/xview.hpp"
/**
 * <!-- ************************************************************************************** -->
 * @defgroup Mesh_Utilities Array and Grid Utilities
 * @brief Functions for generating and processing arrays and grids.
 * @details Declares a set of functions for generating and processing arrays and grids. These functions include:
 *          - Converting boundary arrays to center arrays (linear and logarithmic),
 *          - Creating linearly and logarithmically spaced arrays,
 *          - Creating arrays with uniform spacing in cosine,
 *          - Generating arrays of zeros and ones,
 *          - Finding the minimum and maximum of grids,
 *          - Checking if an array is linearly or logarithmically scaled,
 *          - Creating 2D and 3D grids.
 * <!-- ************************************************************************************** -->
 */

/// Type alias for 1D arrays (e.g., time points)
using Array = xt::xtensor<Real, 1>;
/// Type alias for 2D grids (e.g., spatial coordinates)
using MeshGrid = xt::xtensor<Real, 2>;
/// Type alias for 3D grids (e.g., full spatial-temporal data)
using MeshGrid3d = xt::xtensor<Real, 3>;
/// Type alias for 3D boolean grids for masking operations
using MaskGrid = xt::xtensor<int, 3>;

/**
 * <!-- ************************************************************************************** -->
 * @class Coord
 * @brief Represents a coordinate system with arrays for phi, theta, and t.
 * @details This class is used to define the computational grid for GRB simulations.
 *          It stores the angular coordinates (phi, theta) and time (t) arrays,
 *          along with derived quantities needed for numerical calculations.
 * <!-- ************************************************************************************** -->
 */
class Coord {
   public:
    /// Default constructor
    Coord() = default;

    Array phi;           ///< Array of azimuthal angles (phi) in radians
    Array theta;         ///< Array of polar angles (theta) in radians
    MeshGrid3d t;        ///< Array of engine time points
    Real theta_view{0};  ///< Viewing angle
    Real phi_view{0};    ///< Viewing angle

    /**
     * <!-- ************************************************************************************** -->
     * @brief Returns the dimensions of the coordinate arrays
     * @return Tuple containing (n_phi, n_theta, n_t)
     * <!-- ************************************************************************************** -->
     */
    auto shape() const { return std::make_tuple(phi.size(), theta.size(), t.shape()[2]); }
};

/**
 * <!-- ************************************************************************************** -->
 * @brief Checks if an array is linearly scaled within a given tolerance
 * @param arr The array to check
 * @param tolerance Maximum allowed deviation from linearity (default: 1e-6)
 * @return True if the array is linearly scaled, false otherwise
 * <!-- ************************************************************************************** -->
 */
bool is_linear_scale(Array const& arr, Real tolerance = 1e-6);

/**
 * <!-- ************************************************************************************** -->
 * @brief Checks if an array is logarithmically scaled within a given tolerance
 * @param arr The array to check
 * @param tolerance Maximum allowed deviation from logarithmic scaling (default: 1e-6)
 * @return True if the array is logarithmically scaled, false otherwise
 * <!-- ************************************************************************************** -->
 */
bool is_log_scale(Array const& arr, Real tolerance = 1e-6);

/**
 * <!-- ************************************************************************************** -->
 * @brief Creates an array with logarithmically spaced values
 * @tparam Arr Type of the output array
 * @param start Starting value (log10)
 * @param end Ending value (log10)
 * @param result Output array to store the results
 * @details This is useful for creating time grids where we need more resolution at early times.
 * <!-- ************************************************************************************** -->
 */
template <typename Arr>
void logspace(Real start, Real end, Arr& result);

/**
 * <!-- ************************************************************************************** -->
 * @brief Converts boundary values to center values using linear interpolation
 * @tparam Arr1 Type of the input array
 * @tparam Arr2 Type of the output array
 * @param boundary Input array of boundary values
 * @param center Output array of center values
 * @details This is used to compute cell-centered values from cell-boundary values.
 * <!-- ************************************************************************************** -->
 */
template <typename Arr1, typename Arr2>
void boundary_to_center(Arr1 const& boundary, Arr2& center);

/**
 * <!-- ************************************************************************************** -->
 * @brief Converts boundary values to center values using geometric mean (logarithmic interpolation)
 * @tparam Arr1 Type of the input array
 * @tparam Arr2 Type of the output array
 * @param boundary Input array of boundary values
 * @param center Output array of center values
 * @details This is used when working with logarithmically scaled quantities.
 * <!-- ************************************************************************************** -->
 */
template <typename Arr1, typename Arr2>
void boundary_to_center_log(Arr1 const& boundary, Arr2& center);

/**
 * <!-- ************************************************************************************** -->
 * @brief Converts boundary values to center values using linear interpolation
 * @param boundary Input array of boundary values
 * @return Array of center values
 * <!-- ************************************************************************************** -->
 */
Array boundary_to_center(Array const& boundary);

/**
 * <!-- ************************************************************************************** -->
 * @brief Converts boundary values to center values using geometric mean
 * @param boundary Input array of boundary values
 * @return Array of center values
 * <!-- ************************************************************************************** -->
 */
Array boundary_to_center_log(Array const& boundary);

/**
 * <!-- ************************************************************************************** -->
 * @brief Constructs a coordinate grid (Coord) for shock evolution
 * @tparam Ejecta Type of the jet/ejecta class
 * @param jet The jet/ejecta object
 * @param t_obs Array of observation times
 * @param theta_cut Maximum theta value to include
 * @param theta_view Viewing angle
 * @param z Redshift
 * @param phi_resol Number of grid points per degree in phi (default: 0.1)
 * @param theta_resol Number of grid points per degree in theta (default: 1)
 * @param t_resol Number of grid points per decade in time (default: 5)
 * @param is_axisymmetric Whether the jet is axisymmetric (default: true)
 * @param phi_view Viewing angle (default: 0)
 * @return A Coord object with the constructed grid
 * @details The grid is based on the observation times (t_obs), maximum theta value (theta_cut), and
 *          specified numbers of grid points in phi, theta, and t. The radial grid is logarithmically
 *          spaced between t_min and t_max, and the theta grid is generated linearly.
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta>
Coord auto_grid(Ejecta const& jet, Array const& t_obs, Real theta_cut, Real theta_view, Real z, Real phi_resol = 0.25,
                Real theta_resol = 1, Real t_resol = 3, bool is_axisymmetric = true, Real phi_view = 0);

//========================================================================================================
//                                  template function implementation
//========================================================================================================
template <typename Arr1, typename Arr2>
void boundary_to_center(Arr1 const& boundary, Arr2& center) {
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = 0.5 * (boundary[i] + boundary[i + 1]);
    }
}

template <typename Arr1, typename Arr2>
void boundary_to_center_log(Arr1 const& boundary, Arr2& center) {
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = std::sqrt(boundary[i] * boundary[i + 1]);
    }
}

template <typename Ejecta>
Coord auto_grid(Ejecta const& jet, Array const& t_obs, Real theta_cut, Real theta_view, Real z, Real phi_resol,
                Real theta_resol, Real t_resol, bool is_axisymmetric, Real phi_view) {
    // constexpr size_t min_grid_size = 24;
    Coord coord;
    coord.theta_view = theta_view;

    size_t phi_num = std::max<size_t>(static_cast<size_t>(360 * phi_resol), 1);
    coord.phi = xt::linspace(0., 2 * con::pi, phi_num);  // Generate phi grid linearly spaced.

    Real jet_edge = find_jet_edge(jet, con::Gamma_cut);  // Determine the jet edge angle.
    Real theta_min = 1e-6;
    Real theta_max = std::min(jet_edge, theta_cut);
    size_t theta_num = std::max<size_t>(static_cast<size_t>((theta_max - theta_min) * 180 / con::pi * theta_resol), 32);
    coord.theta = xt::linspace(theta_min, theta_max, theta_num);  //

    Real t_max = *std::max_element(t_obs.begin(), t_obs.end());  // Maximum observation time.
    Real t_min = *std::min_element(t_obs.begin(), t_obs.end());  // Minimum observation time.
    size_t t_num = std::max<size_t>(static_cast<size_t>(std::log10(t_max / t_min) * t_resol), 24);

    size_t phi_size_needed = is_axisymmetric ? 1 : phi_num;
    coord.t = xt::zeros<Real>({phi_size_needed, theta_num, t_num});
    for (size_t i = 0; i < phi_size_needed; ++i) {
        for (size_t j = 0; j < theta_num; ++j) {
            Real b = gamma_to_beta(jet.Gamma0(coord.phi(i), coord.theta(j)));
            // Real theta_max = coord.theta(j) + theta_view;
            Real theta_max = coord.theta.back() + theta_view;
            Real t_start = 0.99 * t_min * (1 - b) / (1 - std::cos(theta_max) * b) / (1 + z);
            Real t_end = 1.01 * t_max / (1 + z);
            xt::view(coord.t, i, j, xt::all()) = xt::logspace(std::log10(t_start), std::log10(t_end), t_num);
        }
    }

    return coord;  // Construct coordinate object.
}
