#ifndef _MESHES_
#define _MESHES_

#define NDEBUG
#include <boost/multi_array.hpp>
#include <cmath>
#include <functional>
#include <vector>

using Array = boost::multi_array<double, 1>;
using MeshGrid = boost::multi_array<double, 2>;
using MeshGrid3d = boost::multi_array<double, 3>;
using UnaryFunc = std::function<double(double)>;
using BinaryFunc = std::function<double(double, double)>;
using TernaryFunc = std::function<double(double, double, double)>;

Array boundaryToCenter(Array const& boundary);
Array boundaryToCenterLog(Array const& boundary);
Array linspace(double start, double end, size_t num);
Array logspace(double start, double end, size_t num);
Array uniform_cos(double start, double end, size_t num);
Array zeros(size_t num);
Array ones(size_t num);
double min(MeshGrid const& grid);
double max(MeshGrid const& grid);
bool isLinearScale(Array const& arr, double tolerance = 1e-6);
bool isLogScale(Array const& arr, double tolerance = 1e-6);
MeshGrid createGrid(size_t theta_size, size_t r_size, double val = 0);
MeshGrid createGridLike(MeshGrid const& grid, double val = 0);
MeshGrid3d create3DGrid(size_t phi_size, size_t theta_size, size_t r_size, double val = 0);
MeshGrid3d create3DGridLike(MeshGrid3d const& grid, double val = 0);

class Coord {
   public:
    Coord(Array const& phi, Array const& theta, Array const& r);
    Coord() = delete;

    Array phi;
    Array theta;
    Array r;

    Array dphi;
    Array dcos;

    double t_min{0};
    double t_max{std::numeric_limits<double>::infinity()};

    auto shape() const { return std::make_tuple(phi.size(), theta.size(), r.size()); }
};

template <typename... MeshGrid>
double min(MeshGrid const&... grids) {
    return std::min(std::min(grids...));
}

template <typename... MeshGrid>
double max(MeshGrid const&... grids) {
    return std::max(std::max(grids...));
}

template <typename Arr>
void linspace(double start, double end, Arr& result) {
    size_t num = result.size();
    double step = (end - start) / (num - 1);
    for (size_t i = 0; i < num; i++) {
        result[i] = start + i * step;
    }
}

template <typename Arr>
void logspace(double start, double end, Arr& result) {
    size_t num = result.size();
    double log_start = std::log(start);
    double log_end = std::log(end);
    double step = (log_end - log_start) / (num - 1);
    for (size_t i = 0; i < num; i++) {
        result[i] = std::exp(log_start + i * step);
    }
}

template <typename Arr1, typename Arr2>
void boundaryToCenter(Arr1 const& boundary, Arr2& center) {
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = 0.5 * (boundary[i] + boundary[i + 1]);
    }
}

template <typename Arr1, typename Arr2>
void boundaryToCenterLog(Arr1 const& boundary, Arr2& center) {
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = std::sqrt(boundary[i] * boundary[i + 1]);
    }
}
#endif