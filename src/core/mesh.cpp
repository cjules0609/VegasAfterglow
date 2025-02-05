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
 * FUNCTION: linspace
 * DESCRIPTION: Creates and returns an Array with 'num' linearly spaced values between 'start' and 'end'.
 ********************************************************************************************************************/
Array linspace(Real start, Real end, size_t num) {
    // Create an Array with the specified size using boost::multi_array extents.
    Array result(boost::extents[num]);

    Real step = (end - start) / (num - 1);
    if (num == 1) {
        step = 0;
    }
    for (size_t i = 0; i < num; i++) {
        // Compute each value by linear interpolation between start and end.
        result[i] = start + i * step;
    }
    return result;
}

/********************************************************************************************************************
 * FUNCTION: logspace
 * DESCRIPTION: Creates and returns an Array with 'num' logarithmically spaced values between 'start' and 'end'.
 ********************************************************************************************************************/
Array logspace(Real start, Real end, size_t num) {
    Array result(boost::extents[num]);
    Real log_start = std::log(start);
    Real log_end = std::log(end);
    Real step = (log_end - log_start) / (num - 1);
    if (num == 1) step = 0;
    for (size_t i = 0; i < num; i++) {
        // Exponentiate the linearly spaced logarithmic values to obtain logarithmically spaced values.
        result[i] = std::exp(log_start + i * step);
    }
    return result;
}

/********************************************************************************************************************
 * FUNCTION: uniform_cos
 * DESCRIPTION: Creates and returns an Array of 'num' values uniformly spaced in cosine between angles 'start' and
 *'end'.
 ********************************************************************************************************************/
Array uniform_cos(Real start, Real end, size_t num) {
    // First generate a linearly spaced array in cosine space.
    Array result = linspace(std::cos(start), std::cos(end), num);
    for (size_t i = 0; i < num; i++) {
        // Convert back to angle by taking the arccosine.
        result[i] = std::acos(result[i]);
    }
    return result;
}

/********************************************************************************************************************
 * FUNCTION: zeros
 * DESCRIPTION: Creates and returns an Array of size 'num' initialized to zero.
 ********************************************************************************************************************/
Array zeros(size_t num) {
    Array result(boost::extents[num]);
    std::fill(result.data(), result.data() + result.num_elements(), 0);
    return result;
}

/********************************************************************************************************************
 * FUNCTION: ones
 * DESCRIPTION: Creates and returns an Array of size 'num' filled with ones.
 ********************************************************************************************************************/
Array ones(size_t num) {
    Array array(boost::extents[num]);
    std::fill(array.data(), array.data() + array.num_elements(), 1);
    return array;
}

/********************************************************************************************************************
 * FUNCTION: createGrid
 * DESCRIPTION: Creates and returns a 2D MeshGrid with dimensions (theta_size x r_size) filled with the value 'val'.
 ********************************************************************************************************************/
MeshGrid createGrid(size_t theta_size, size_t r_size, Real val) {
    MeshGrid grid(boost::extents[theta_size][r_size]);
    std::fill(grid.data(), grid.data() + grid.num_elements(), val);
    return grid;
}

/********************************************************************************************************************
 * FUNCTION: isLinearScale
 * DESCRIPTION: Checks if the values in the given Array are approximately linearly spaced within the specified
 *tolerance.
 ********************************************************************************************************************/
bool isLinearScale(Array const& arr, Real tolerance) {
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
 * FUNCTION: isLogScale
 * DESCRIPTION: Checks if the values in the given Array are approximately logarithmically spaced (constant ratio)
 *              within the specified tolerance.
 ********************************************************************************************************************/
bool isLogScale(Array const& arr, Real tolerance) {
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
 * FUNCTION: createGridLike
 * DESCRIPTION: Creates a 2D MeshGrid with the same shape as the given grid 'grid_old', filled with the value 'val'.
 ********************************************************************************************************************/
MeshGrid createGridLike(MeshGrid const& grid_old, Real val) {
    const size_t* shape = grid_old.shape();
    MeshGrid grid(boost::extents[shape[0]][shape[1]]);
    std::fill(grid.data(), grid.data() + grid.num_elements(), val);
    return grid;
}

/********************************************************************************************************************
 * FUNCTION: create3DGrid
 * DESCRIPTION: Creates and returns a 3D MeshGrid (MeshGrid3d) with dimensions (phi_size x theta_size x r_size)
 *              filled with the value 'val'.
 ********************************************************************************************************************/
MeshGrid3d create3DGrid(size_t phi_size, size_t theta_size, size_t r_size, Real val) {
    MeshGrid3d grid(boost::extents[phi_size][theta_size][r_size]);
    std::fill(grid.data(), grid.data() + grid.num_elements(), val);
    return grid;
}

/********************************************************************************************************************
 * FUNCTION: create3DGridLike
 * DESCRIPTION: Creates a 3D MeshGrid with the same shape as the provided grid 'grid_old', filled with the value 'val'.
 ********************************************************************************************************************/
MeshGrid3d create3DGridLike(MeshGrid3d const& grid_old, Real val) {
    const size_t* shape = grid_old.shape();
    MeshGrid3d grid(boost::extents[shape[0]][shape[1]][shape[2]]);
    std::fill(grid.data(), grid.data() + grid.num_elements(), val);
    return grid;
}

/********************************************************************************************************************
 * CONSTRUCTOR: Coord::Coord
 * DESCRIPTION: Constructs a Coord object with the provided phi, theta, and r arrays. It also computes the
 *              differential arrays dphi and dcos for phi and theta, respectively.
 ********************************************************************************************************************/
Coord::Coord(Array const& phi, Array const& theta, Array const& r)
    : phi(phi), theta(theta), r(r), dphi(zeros(phi.size())), dcos(zeros(theta.size())) {
    size_t phi_size = phi.size();
    if (phi_size > 1) {
        // Compute dphi for the first element.
        dphi[0] = std::fabs(phi[1] - phi[0]) / 2;
        // Compute dphi for intermediate elements.
        for (size_t i = 1; i < phi_size - 1; i++) {
            dphi[i] = std::fabs(phi[i + 1] - phi[i - 1]) / 2;
        }
        // Compute dphi for the last element.
        dphi[phi_size - 1] = std::fabs(phi[phi_size - 1] - phi[phi_size - 2]) / 2;
    } else if (phi_size == 1) {
        dphi[0] = 2 * con::pi;
    }

    size_t theta_size = theta.size();
    if (theta_size > 1) {
        // Compute dcos for the first element.
        Real theta_hi = (theta[1] + theta[0]) / 2;
        dcos[0] = std::fabs(std::cos(theta_hi) - std::cos(theta[0]));

        // Compute dcos for intermediate theta values.
        for (size_t j = 1; j < theta_size - 1; j++) {
            Real theta_lo = (theta[j] + theta[j - 1]) / 2;
            Real theta_hi = (theta[j] + theta[j + 1]) / 2;
            dcos[j] = std::fabs(std::cos(theta_hi) - std::cos(theta_lo));
        }
        // Compute dcos for the last element.
        Real theta_lo = (theta[theta_size - 1] + theta[theta_size - 2]) / 2;
        dcos[theta_size - 1] = std::fabs(std::cos(theta[theta_size - 1]) - std::cos(theta_lo));
    } else if (theta_size == 1) {
        dcos[0] = 1;
    }
}

/********************************************************************************************************************
 * FUNCTION: min (MeshGrid)
 * DESCRIPTION: Returns the minimum value found in the provided 2D MeshGrid.
 ********************************************************************************************************************/
Real min(MeshGrid const& grid) {
    Real min = grid[0][0];
    for (size_t i = 0; i < grid.size(); ++i) {
        for (size_t j = 0; j < grid[i].size(); ++j) {
            if (grid[i][j] < min) {
                min = grid[i][j];
            }
        }
    }
    return min;
}

/********************************************************************************************************************
 * FUNCTION: max (MeshGrid)
 * DESCRIPTION: Returns the maximum value found in the provided 2D MeshGrid.
 ********************************************************************************************************************/
Real max(MeshGrid const& grid) {
    Real max = grid[0][0];
    for (size_t i = 0; i < grid.size(); ++i) {
        for (size_t j = 0; j < grid[i].size(); ++j) {
            if (grid[i][j] > max) {
                max = grid[i][j];
            }
        }
    }
    return max;
}

/********************************************************************************************************************
 * FUNCTION: boundaryToCenter (linear)
 * DESCRIPTION: Converts a boundary array to center values by averaging adjacent boundaries.
 ********************************************************************************************************************/
Array boundaryToCenter(Array const& boundary) {
    Array center(boost::extents[boundary.size() - 1]);
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = 0.5 * (boundary[i] + boundary[i + 1]);
    }
    return center;
}

/********************************************************************************************************************
 * FUNCTION: boundaryToCenterLog
 * DESCRIPTION: Converts a boundary array to center values in logarithmic space using the geometric mean.
 ********************************************************************************************************************/
Array boundaryToCenterLog(Array const& boundary) {
    Array center(boost::extents[boundary.size() - 1]);
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = std::sqrt(boundary[i] * boundary[i + 1]);
    }
    return center;
}