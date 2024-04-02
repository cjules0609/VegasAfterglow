#include "mesh.h"

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
    double log_start = std::log10(start);
    double log_end = std::log10(end);
    double step = (log_end - log_start) / (num - 1);
    for (size_t i = 0; i < num; i++) {
        result[i] = std::pow(10, log_start + i * step);
    }
    return result;
}

Array zeros(size_t num) {
    Array result(num);
    std::fill(result.begin(), result.end(), 0);
    return result;
}

Array ones(size_t num) {
    Array result(num);
    std::fill(result.begin(), result.end(), 1);
    return result;
}

MeshGrid createGrid(size_t theta_size, size_t r_size, double val = 0) {
    return MeshGrid(theta_size, Array(r_size, val));
}

MeshGrid3d createGrid3d(size_t phi_size, size_t theta_size, size_t r_size, double val = 0) {
    return MeshGrid3d(phi_size, MeshGrid(theta_size, Array(r_size, val)));
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