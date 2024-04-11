#include "mesh.h"

#include <cmath>

Array linspace(double start, double end, size_t num) {
    Array result(num);
    double step = (end - start) / (num - 1);
    for (size_t i = 0; i < num; i++) {
        result[i] = start + i * step;
    }
    return result;
}

Array logspace(double start, double end, size_t num) {
    Array result(num);
    double log_start = std::log(start);
    double log_end = std::log(end);
    double step = (log_end - log_start) / (num - 1);
    for (size_t i = 0; i < num; i++) {
        result[i] = std::exp(log_start + i * step);
    }
    return result;
}

Array zeros(size_t num) { return Array(num, 0); }

Array ones(size_t num) { return Array(num, 1); }

MeshGrid create_grid(size_t theta_size, size_t r_size, double val) { return MeshGrid(theta_size, Array(r_size, val)); }

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

MeshGrid create_grid_like(MeshGrid const& grid, double val) {
    return MeshGrid(grid.size(), Array(grid[0].size(), val));
}

MeshGrid3d create_3d_grid(size_t phi_size, size_t theta_size, size_t r_size, double val) {
    return MeshGrid3d(phi_size, MeshGrid(theta_size, Array(r_size, val)));
}

MeshGrid3d create_3d_grid_like(MeshGrid3d const& grid, double val) {
    return MeshGrid3d(grid.size(), MeshGrid(grid[0].size(), Array(grid[0][0].size(), val)));
}

Array boundary2center(Array const& boundary) {
    Array center(boundary.size() - 1);
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = 0.5 * (boundary[i] + boundary[i + 1]);
    }
    return center;
}

Array boundary2centerlog(Array const& boundary) {
    Array center(boundary.size() - 1);
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = std::sqrt(boundary[i] * boundary[i + 1]);
    }
    return center;
}