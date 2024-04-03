#ifndef _SYNCHROTRON_
#define _SYNCHROTRON_
#include <vector>

#include "medium.h"
#include "mesh.h"
struct SynRad {
    double I_nu_peak{0};
    double nu_E_peak{0};
    double nu_m{0};
    double nu_c{0};
    double nu_a{0};
    double nu_M{0};
    double pel{2.3};

    double I_nu(double nu) const;

   private:
    inline double I_nu_(double nu) const;
};
using SynRadArray = std::vector<SynRad>;
using SynRadMesh = std::vector<std::vector<SynRad>>;

SynRadMesh createSynRadGrid(size_t theta_size, size_t r_size, SynRad val = {0, 0, 0, 0, 0, 0, 2.3});

SynRadMesh calc_syn_radiation(Coord const& coord, MeshGrid const& Gamma, MeshGrid const& t_com, Medium const& medium);

double calc_syn_gamma_c(double Gamma, double t_com, double Bprime, double Y_tilt);
double calc_syn_gamma_M(double Bprime, double zeta, double Y_tilt);
double syn_nu(double gamma, double B);
double syn_gamma(double nu, double B);
double syn_nu_E_peak(SynRad const& rad);
MeshGrid get_B_field(Coord const& coord, MeshGrid const& Gamma, Medium const& medium);
#endif