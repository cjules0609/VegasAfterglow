#ifndef _SYNCHROTRON_
#define _SYNCHROTRON_
#include <vector>

#include "medium.h"
#include "mesh.h"
#include "shock.h"
struct SynElectrons {
    double gamma_m{0};
    double gamma_c{0};
    double gamma_a{0};
    double gamma_M{0};
    double p{2.3};
};

struct SynPhotons {
    double I_nu_peak{0};
    double nu_E_peak{0};
    double nu_m{0};
    double nu_c{0};
    double nu_a{0};
    double nu_M{0};
    double p{2.3};

    double I_nu(double nu) const;

   private:
    inline double I_nu_(double nu) const;
};

/*struct SynRad {
    SynElectronsMesh electrons;
    SynPhotonsMesh photons;
};*/

using SynPhotonsArray = std::vector<SynPhotons>;
using SynPhotonsMesh = std::vector<std::vector<SynPhotons>>;

using SynElectronsArray = std::vector<SynElectrons>;
using SynElectronsMesh = std::vector<std::vector<SynElectrons>>;

SynPhotonsMesh create_syn_photons_grid(size_t theta_size, size_t r_size);

SynElectronsMesh create_syn_electrons_grid(size_t theta_size, size_t r_size);

SynElectronsMesh gen_syn_electrons(Coord const& coord, Shock const& shock, Medium const& medium);

SynPhotonsMesh gen_syn_photons(Coord const& coord, SynElectronsMesh const& electrons, Shock const& shock,
                               Medium const& medium);

double syn_gamma_c(double Gamma, double t_com, double B, double Y_tilt);
double syn_gamma_N_peak(double gamma_a, double gamma_m, double gamma_c);
double syn_gamma_M(double B, double zeta, double Y_tilt);
/*
double syn_nu(double gamma, double B);
double syn_gamma(double nu, double B);
double syn_nu_E_peak(SynPhotons const& rad);*/
#endif