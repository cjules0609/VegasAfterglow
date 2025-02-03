#include "mesh.h"

#include <cmath>

#include "macros.h"
Array linspace(double start, double end, size_t num) {
    Array result(boost::extents[num]);
    double step = (end - start) / (num - 1);
    for (size_t i = 0; i < num; i++) {
        result[i] = start + i * step;
    }
    return result;
}

Array logspace(double start, double end, size_t num) {
    Array result(boost::extents[num]);
    double log_start = std::log(start);
    double log_end = std::log(end);
    double step = (log_end - log_start) / (num - 1);
    for (size_t i = 0; i < num; i++) {
        result[i] = std::exp(log_start + i * step);
    }
    return result;
}

Array uniform_cos(double start, double end, size_t num) {
    Array result = linspace(std::cos(start), std::cos(end), num);
    for (size_t i = 0; i < num; i++) {
        result[i] = std::acos(result[i]);
    }
    return result;
}

Array zeros(size_t num) { return Array(boost::extents[num]); }

Array ones(size_t num) {
    Array array(boost::extents[num]);
    std::fill(array.data(), array.data() + array.num_elements(), 1);
    return array;
}

MeshGrid createGrid(size_t theta_size, size_t r_size, double val) {
    MeshGrid grid(boost::extents[theta_size][r_size]);
    std::fill(grid.data(), grid.data() + grid.num_elements(), val);
    return grid;
}

bool isLinearScale(Array const& arr, double tolerance) {
    if (arr.size() < 2) return false;  // Need at least two elements

    double diff = arr[1] - arr[0];
    for (size_t i = 2; i < arr.size(); ++i) {
        if (std::fabs((arr[i] - arr[i - 1] - diff) / diff) > tolerance) {
            return false;
        }
    }
    return true;
}

bool isLogScale(Array const& arr, double tolerance) {
    if (arr.size() < 2) return false;  // Need at least two elements

    double ratio = arr[1] / arr[0];
    for (size_t i = 2; i < arr.size(); ++i) {
        if (std::fabs((arr[i] / arr[i - 1] - ratio) / ratio) > tolerance) {
            return false;
        }
    }
    return true;
}

Coord::Coord(Array const& phi, Array const& theta, Array const& r)
    : phi(phi), theta(theta), r(r), dphi(zeros(phi.size())), dcos(zeros(theta.size())) {
    size_t phi_size = phi.size();
    if (phi_size > 1) {
        dphi[0] = std::fabs(phi[1] - phi[0]) / 2;
        for (size_t i = 1; i < phi_size - 1; i++) {
            dphi[i] = std::fabs(phi[i + 1] - phi[i - 1]) / 2;
        }
        dphi[phi_size - 1] = std::fabs(phi[phi_size - 1] - phi[phi_size - 2]) / 2;
    }

    size_t theta_size = theta.size();
    if (theta_size > 1) {
        double theta_hi = (theta[1] + theta[0]) / 2;
        dcos[0] = std::fabs(std::cos(theta_hi) - std::cos(theta[0]));

        for (size_t j = 1; j < theta_size - 1; j++) {
            double theta_lo = (theta[j] + theta[j - 1]) / 2;
            double theta_hi = (theta[j] + theta[j + 1]) / 2;
            dcos[j] = std::fabs(std::cos(theta_hi) - std::cos(theta_lo));
        }
        double theta_lo = (theta[theta_size - 1] + theta[theta_size - 2]) / 2;
        dcos[theta_size - 1] = std::fabs(std::cos(theta[theta_size - 1]) - std::cos(theta_lo));
    }
}

double min(MeshGrid const& grid) {
    double min = grid[0][0];
    for (size_t i = 0; i < grid.size(); ++i) {
        for (size_t j = 0; j < grid[i].size(); ++j) {
            if (grid[i][j] < min) {
                min = grid[i][j];
            }
        }
    }
    return min;
}

double max(MeshGrid const& grid) {
    double max = grid[0][0];
    for (size_t i = 0; i < grid.size(); ++i) {
        for (size_t j = 0; j < grid[i].size(); ++j) {
            if (grid[i][j] > max) {
                max = grid[i][j];
            }
        }
    }
    return max;
}

MeshGrid createGridLike(MeshGrid const& grid_old, double val) {
    const size_t* shape = grid_old.shape();
    MeshGrid grid(boost::extents[shape[0]][shape[1]]);
    std::fill(grid.data(), grid.data() + grid.num_elements(), val);
    return grid;
}

MeshGrid3d create3DGrid(size_t phi_size, size_t theta_size, size_t r_size, double val) {
    MeshGrid3d grid(boost::extents[phi_size][theta_size][r_size]);
    std::fill(grid.data(), grid.data() + grid.num_elements(), val);
    return grid;
}

MeshGrid3d create3DGridLike(MeshGrid3d const& grid_old, double val) {
    const size_t* shape = grid_old.shape();
    MeshGrid3d grid(boost::extents[shape[0]][shape[1]][shape[2]]);
    std::fill(grid.data(), grid.data() + grid.num_elements(), val);
    return grid;
}

Array boundaryToCenter(Array const& boundary) {
    Array center(boost::extents[boundary.size() - 1]);
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = 0.5 * (boundary[i] + boundary[i + 1]);
    }
    return center;
}

Array boundaryToCenterLog(Array const& boundary) {
    Array center(boost::extents[boundary.size() - 1]);
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = std::sqrt(boundary[i] * boundary[i + 1]);
    }
    return center;
}
