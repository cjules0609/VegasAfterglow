#ifndef _SYNCHROTRON_
#define _SYNCHROTRON_
#include <vector>

#include "medium.h"
#include "mesh.h"
#include "shock.h"
struct SynElectron {
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
using SynElectronArray = std::vector<SynElectron>;
using SynElectronMesh = std::vector<std::vector<SynElectron>>;

SynElectronMesh create_syn_electron_grid(size_t theta_size, size_t r_size, SynElectron val = {0, 0, 0, 0, 0, 0, 2.3});

SynElectronMesh gen_syn_electrons(Coord const& coord, Shock const& shock, Medium const& medium);

double syn_gamma_c(double Gamma, double t_com, double Bprime, double Y_tilt);
double syn_gamma_M(double Bprime, double zeta, double Y_tilt);
double syn_nu(double gamma, double B);
double syn_gamma(double nu, double B);
double syn_nu_E_peak(SynElectron const& rad);

#endif