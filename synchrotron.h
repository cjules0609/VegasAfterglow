#ifndef _SYNCHROTRON_
#define _SYNCHROTRON_
#include <vector>

#include "medium.h"
#include "mesh.h"
#include "shock.h"
struct SynElectrons {
    // all in comoving frame
    double gamma_m{0};
    double gamma_c{0};
    double gamma_a{0};
    double gamma_M{0};
    double p{2.3};
    double N_tot{0};
    double n_tot{0};
    double gamma_N_peak;
    size_t regime{0};
    double N(double gamma) const;
    double n(double gamma) const;

   private:
    inline double gamma_spectrum_(double gamma) const;
};

struct SynPhotons {
    // all in comoving frame
    double L_nu_peak{0};
    double nu_E_peak{0};
    double nu_m{0};
    double nu_c{0};
    double nu_a{0};
    double nu_M{0};
    double p{2.3};
    size_t regime{0};
    double L_nu(double nu) const;

   private:
    inline double spectrum_(double nu) const;
};

using SynPhotonsArray = std::vector<SynPhotons>;
using SynPhotonsMesh = std::vector<std::vector<SynPhotons>>;

using SynElectronsArray = std::vector<SynElectrons>;
using SynElectronsMesh = std::vector<std::vector<SynElectrons>>;

SynPhotonsMesh create_syn_photons_grid(size_t theta_size, size_t r_size);

SynElectronsMesh create_syn_electrons_grid(size_t theta_size, size_t r_size);

SynElectronsMesh gen_syn_electrons(double p, Coord const& coord, Shock const& shock);

SynElectronsMesh gen_syn_electrons(double p, Coord const& coord, Shock const& shock, MeshGrid const& Y_tilt);

SynPhotonsMesh gen_syn_photons(SynElectronsMesh const& electrons, Coord const& coord, Shock const& shock);

double syn_gamma_c(double t_com, double B, double Y_tilt);
double syn_gamma_N_peak(double gamma_a, double gamma_m, double gamma_c);
double syn_gamma_M(double B, double zeta, double Y_tilt);
double syn_gamma_a(double Gamma, double B, double I_syn_peak, double gamma_m, double gamma_c, double gamma_M);
// double syn_j_nu_peak(double r, double Gamma, double B, double rho, double xi, double p);
double syn_nu(double gamma, double B);
#endif