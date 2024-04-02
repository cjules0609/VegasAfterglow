#ifndef _SYNCHROTRON_
#define _SYNCHROTRON_
#include <vector>

#include "medium.h"
#include "mesh.h"
struct SynRad {
    double I_nu_peak;
    double nu_m;
    double nu_c;
    double nu_a;
    double nu_M;

    double I_nu(double nu, double pel) const;

   private:
    inline double I_nu_(double nu, double pel) const;
};
using SynRadArray = std::vector<SynRad>;
using SynRadMesh = std::vector<std::vector<SynRad>>;

SynRadMesh createSynRadGrid(size_t theta_size, size_t r_size, SynRad val = {0, 0, 0, 0, 0});

SynRadMesh calc_syn_radiation(Coord const& coord, MeshGrid const& Gamma, MeshGrid const& t_com, Medium const& medium);
#endif