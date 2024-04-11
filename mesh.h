#ifndef _MESHES_
#define _MESHES_

#include <cmath>
#include <functional>
#include <vector>

using MeshGrid = std::vector<std::vector<double>>;

using MeshGrid3d = std::vector<std::vector<std::vector<double>>>;

using Array = std::vector<double>;

using Profile = std::function<double(double)>;

using Profile2d = std::function<double(double, double)>;

Array linspace(double start, double end, size_t num);

Array logspace(double start, double end, size_t num);

Array zeros(size_t num);

Array ones(size_t num);

MeshGrid create_grid(size_t theta_size, size_t r_size, double val = 0);

double min(MeshGrid const& grid);

double max(MeshGrid const& grid);

MeshGrid create_grid_like(MeshGrid const& grid, double val = 0);

MeshGrid3d create_3d_grid(size_t phi_size, size_t theta_size, size_t r_size, double val = 0);

MeshGrid3d create_3d_grid_like(MeshGrid3d const& grid, double val = 0);

Array boundary2center(Array const& boundary);

Array boundary2centerlog(Array const& boundary);

class Coord {
   public:
    Coord(Array const& r_b, Array const& theta_b, Array const& phi_b);

    Coord(double r_min, double r_max, double theta_max, size_t r_num, size_t theta_num, size_t phi_num);

    Coord() = delete;
    Array r_b;
    Array theta_b;
    Array phi_b;
    Array r;
    Array theta;
    Array phi;
};

template <typename... MeshGrid>
double min(MeshGrid const&... grids) {
    return std::min(min(grids...));
}

template <typename... MeshGrid>
double max(MeshGrid const&... grids) {
    return std::max(max(grids...));
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
void boundary2center(Arr1 const& boundary, Arr2& center) {
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = 0.5 * (boundary[i] + boundary[i + 1]);
    }
}
template <typename Arr1, typename Arr2>
void boundary2centerlog(Arr1 const& boundary, Arr2& center) {
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = std::sqrt(boundary[i] * boundary[i + 1]);
    }
}
#endif