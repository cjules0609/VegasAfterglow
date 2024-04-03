#ifndef _MESHES_
#define _MESHES_

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

template <typename... MeshGrid>
double min(MeshGrid const&... grids) {
    return std::min(min(grids...));
}

template <typename... MeshGrid>
double max(MeshGrid const&... grids) {
    return std::max(max(grids...));
}

MeshGrid create_grid_like(MeshGrid const& grid, double val = 0);

MeshGrid3d create_3d_grid(size_t phi_size, size_t theta_size, size_t r_size, double val = 0);

MeshGrid3d create_3d_grid_like(MeshGrid3d const& grid, double val = 0);

Array boundary2center(Array const& boundary);

Array boundary2centerlog(Array const& boundary);

class Coord {
   public:
    Coord(Array r_b, Array theta_b, Array phi_b) : r_b(r_b), theta_b(theta_b), phi_b(phi_b) {
        r = boundary2centerlog(r_b);
        theta = boundary2center(theta_b);
        phi = boundary2center(phi_b);
    }
    Coord() = delete;
    Array r_b;
    Array theta_b;
    Array phi_b;
    Array r;
    Array theta;
    Array phi;
};

#endif