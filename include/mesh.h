//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <cmath>

#include "boost/numeric/odeint.hpp"
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
 * @param phi_ppd Points per degree in phi direction
 * @param theta_ppd Points per degree in theta direction
 * @param t_ppd Points per decade in time direction
 * @param is_axisymmetric Whether the jet is axisymmetric (default: true)
 * @param phi_view Viewing angle (default: 0)
 * @return A Coord object with the constructed grid
 * @details The grid is based on the observation times (t_obs), maximum theta value (theta_cut), and
 *          specified numbers of grid points in phi, theta, and t. The radial grid is logarithmically
 *          spaced between t_min and t_max, and the theta grid is generated linearly.
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta>
Coord auto_grid(Ejecta const& jet, Array const& t_obs, Real theta_cut, Real theta_view, Real z, Real phi_ppd = 0.3,
                Real theta_ppd = 1, Real t_ppd = 5, bool is_axisymmetric = true, Real phi_view = 0);

/**
 * <!-- ************************************************************************************** -->
 * @brief Determines the edge of the jet based on a given gamma cut-off using binary search
 * @tparam Ejecta Type of the jet/ejecta class
 * @param jet The jet/ejecta object
 * @param gamma_cut Lorentz factor cutoff value
 * @param phi_ppd Azimuthal resolution (points per degree)
 * @param theta_ppd Polar resolution (points per degree)
 * @param is_axisymmetric Flag for axisymmetric jets
 * @return Angle (in radians) at which the jet's Lorentz factor drops to gamma_cut
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta>
Real find_jet_edge(Ejecta const& jet, Real gamma_cut, Real phi_ppd, Real theta_ppd, bool is_axisymmetric);

/**
 * <!-- ************************************************************************************** -->
 * @brief Determines the edge of the jet where the spreading is strongest
 * @tparam Ejecta Type of the jet/ejecta class
 * @tparam Medium Type of the ambient medium
 * @param jet The jet/ejecta object
 * @param medium The ambient medium object
 * @param phi Azimuthal angle
 * @param theta_min Minimum polar angle to consider
 * @param theta_max Maximum polar angle to consider
 * @param t0 Initial time
 * @return Angle (in radians) where the spreading is strongest
 * @details The spreading strength is measured by the derivative of the pressure with respect to theta,
 *          which is proportional to d((Gamma-1)Gamma rho)/dtheta
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium>
Real jet_spreading_edge(Ejecta const& jet, Medium const& medium, Real phi, Real theta_min, Real theta_max, Real t0);
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

template <typename Arr>
void logspace_boundary_center(Real lg2_min, Real lg2_max, size_t size, Arr& center, Array& bin_width) {
    center = Arr::from_shape({size});
    bin_width = Array::from_shape({size});
    if (size == 0) return;

    const Real dlg2 = (lg2_max - lg2_min) / static_cast<Real>(size);

    const Real r = std::exp2(dlg2);
    const Real s = std::sqrt(r);
    const Real w = r - 1.;

    Real left = std::exp2(lg2_min);

    for (std::size_t i = 0; i < size; ++i) {
        center(i) = left * s;
        bin_width(i) = left * w;
        left *= r;
    }
}

template <typename Ejecta>
Real find_jet_edge(Ejecta const& jet, Real gamma_cut, Real phi_resol, Real theta_resol, bool is_axisymmetric) {
    // binary search for the edge of the jet
    if (jet.Gamma0(0, con::pi / 2) >= gamma_cut) {
        return con::pi / 2;  // If the Lorentz factor at pi/2 is above the cut, the jet extends to pi/2.
    }
    Real low = 0;
    Real hi = con::pi / 2;
    Real eps = 1e-9;
    for (; hi - low > eps;) {
        Real mid = 0.5 * (low + hi);
        if (jet.Gamma0(0, mid) > gamma_cut) {
            low = mid;
        } else {
            hi = mid;
        }
    }

    // grid based search for the edge of the jet
    size_t phi_num = std::max<size_t>(static_cast<size_t>(360. * phi_resol), 1);
    phi_num = is_axisymmetric ? 1 : phi_num;
    size_t theta_num = std::max<size_t>(static_cast<size_t>(90. * theta_resol), 32);
    auto phi = xt::linspace(0., 2 * con::pi, phi_num);
    auto theta = xt::linspace(0., con::pi / 2, theta_num);

    Real theta_edge = con::pi / 2;

    for (size_t i = 0; i < phi_num; ++i) {
        for (int j = theta_num - 1; j >= 0; --j) {
            if (jet.Gamma0(phi[i], theta[j]) >= gamma_cut) {
                theta_edge = theta[j];
                break;  // Found the edge for this phi, no need to check lower theta values
            } else {
                theta_edge = theta[j];
            }
        }
    }

    return std::max(theta_edge, low);
}

template <typename Ejecta, typename Medium>
Real jet_spreading_edge(Ejecta const& jet, Medium const& medium, Real phi, Real theta_min, Real theta_max, Real t0) {
    Real step = (theta_max - theta_min) / 256;
    Real theta_s = theta_min;
    Real dp_min = 0;

    for (Real theta = theta_min; theta <= theta_max; theta += step) {
        // Real G = jet.Gamma0(phi, theta);
        // Real beta0 = gamma_to_beta(G);
        // Real r0 = beta0 * con::c * t0 / (1 - beta0);
        // Real rho = medium.rho(phi, theta, 0);
        Real th_lo = std::max(theta - step, theta_min);
        Real th_hi = std::min(theta + step, theta_max);
        Real dG = (jet.Gamma0(phi, th_hi) - jet.Gamma0(phi, th_lo)) / (th_hi - th_lo);
        // Real drho = (medium.rho(phi, th_hi, r0) - medium.rho(phi, th_lo, r0)) / (th_hi - th_lo);
        Real dp = dG;  //(2 * G - 1) * rho * dG + (G - 1) * G * drho;

        if (dp < dp_min) {
            dp_min = dp;
            theta_s = theta;
        }
    }
    if (dp_min == 0) {
        theta_s = theta_max;
    }

    return theta_s;
}

template <typename Func>
Array inverse_CFD_sampling(Func&& pdf, Real min, Real max, size_t num) {
    using namespace boost::numeric::odeint;
    constexpr Real rtol = 1e-6;
    constexpr size_t sample_num = 200;
    Array x_i = xt::linspace(min, max, sample_num);
    Array CDF_i = xt::zeros<Real>({sample_num});

    auto stepper = make_dense_output(rtol, rtol, runge_kutta_dopri5<Real>());
    stepper.initialize(0, min, (max - min) / 1e3);

    for (size_t k = 1; stepper.current_time() <= max;) {
        stepper.do_step(pdf);
        while (k < x_i.size() && stepper.current_time() > x_i(k)) {
            stepper.calc_state(x_i(k), CDF_i(k));
            ++k;
        }
    }

    Array CDF_out = xt::linspace(CDF_i.front(), CDF_i.back(), num);
    Array x_out = xt::zeros<Real>({num});

    for (size_t k = 0; k < num; ++k) {
        for (size_t j = 0; j < sample_num; ++j) {
            if (CDF_out(k) <= CDF_i(j)) {
                if (j == 0) {
                    x_out(k) = x_i(j);
                } else {
                    Real slope = (x_i(j) - x_i(j - 1)) / (CDF_i(j) - CDF_i(j - 1));
                    x_out(k) = x_i(j - 1) + slope * (CDF_out(k) - CDF_i(j - 1));
                }
                break;
            }
        }
    }
    return x_out;
}

template <typename Ejecta>
Array adaptive_theta(Ejecta const& jet, Real theta_min, Real theta_max, size_t theta_num, Real theta_v) {
    auto eqn = [=, &jet](Real const& cdf, Real& pdf, Real theta) {
        Real Gamma = jet.Gamma0(0, theta);
        Real beta = std::sqrt(std::fabs(Gamma * Gamma - 1)) / Gamma;
        Real a = (1 - beta) / (1 - beta * std::cos(theta - theta_v));
        pdf = a * Gamma * std::sqrt(Gamma * (Gamma - 1)) * std::sin(theta);
    };

    return inverse_CFD_sampling(eqn, theta_min, theta_max, theta_num);
}

template <typename Ejecta>
Array adaptive_phi(Ejecta const& jet, size_t phi_num, Real theta_v, Real theta_max, bool is_axisymmetric) {
    if ((theta_v == 0 && is_axisymmetric) || theta_v > theta_max) {
        return xt::linspace(0., 2 * con::pi, phi_num);
    } else {
        Real cos2 = std::cos(theta_v);
        cos2 = cos2 * cos2;
        Real sin2 = 1 - cos2;

        auto eqn = [=, &jet](Real const& cdf, Real& pdf, Real phi) {
            Real Gamma = jet.Gamma0(phi, theta_v);
            Real beta = std::sqrt(std::fabs(Gamma * Gamma - 1)) / Gamma;
            Real a = (1 - beta) / (1 - beta * (std::cos(phi) * sin2 + cos2));
            pdf = a * Gamma * std::sqrt(Gamma * (Gamma - 1));
        };

        return inverse_CFD_sampling(eqn, 0, 2 * con::pi, phi_num);
    }
}

template <typename Ejecta>
Coord auto_grid(Ejecta const& jet, Array const& t_obs, Real theta_cut, Real theta_view, Real z, Real phi_resol,
                Real theta_resol, Real t_resol, bool is_axisymmetric, Real phi_view) {
    Coord coord;
    coord.theta_view = theta_view;

    Real jet_edge = find_jet_edge(jet, con::Gamma_cut, phi_resol, theta_resol, is_axisymmetric);
    Real theta_min = 1e-6;
    Real theta_max = std::min(jet_edge, theta_cut);
    size_t theta_num = std::max<size_t>(static_cast<size_t>((theta_max - theta_min) * 180 / con::pi * theta_resol), 64);
    // coord.theta = xt::linspace(theta_min, theta_max, theta_num);  //
    coord.theta = adaptive_theta(jet, theta_min, theta_max, theta_num, theta_view);

    size_t phi_num = std::max<size_t>(static_cast<size_t>(360 * phi_resol), 1);
    // coord.phi = xt::linspace(0., 2 * con::pi, phi_num);
    coord.phi = adaptive_phi(jet, phi_num, theta_view, theta_max, is_axisymmetric);

    Real t_max = *std::max_element(t_obs.begin(), t_obs.end());
    Real t_min = *std::min_element(t_obs.begin(), t_obs.end());
    size_t t_num = std::max<size_t>(static_cast<size_t>(std::log10(t_max / t_min) * t_resol), 24);

    size_t phi_size_needed = is_axisymmetric ? 1 : phi_num;
    coord.t = xt::zeros<Real>({phi_size_needed, theta_num, t_num});
    for (size_t i = 0; i < phi_size_needed; ++i) {
        for (size_t j = 0; j < theta_num; ++j) {
            Real b = gamma_to_beta(jet.Gamma0(coord.phi(i), coord.theta(j)));
            // Real theta_max = coord.theta(j) + theta_view;
            Real theta_max = coord.theta.back() + theta_view;
            Real t_start =
                std::max<Real>(0.99 * t_min * (1 - b) / (1 - std::cos(theta_max) * b) / (1 + z), 1e-2 * unit::sec);
            Real t_end = 1.01 * t_max / (1 + z);
            xt::view(coord.t, i, j, xt::all()) = xt::logspace(std::log10(t_start), std::log10(t_end), t_num);
        }
    }

    return coord;  // Construct coordinate object.
}