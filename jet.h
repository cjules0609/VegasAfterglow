#ifndef _JET_
#define _JET_

#include "mesh.h"

class Jet {
   public:
    Profile2d dEdOmega;
    Profile Gamma0;
};

static Profile2d noInjection = Profile2d([](double theta, double t_lab) { return 0; });

Jet create_isotropic_jet(double E_iso, double Gamma0, Profile2d inject = noInjection);

Jet create_tophat_jet(double E_iso, double Gamma0, double theta_c, Profile2d inject = noInjection);

Jet create_power_law_jet(double E_iso, double Gamma0, double theta_m, double k, Profile2d inject = noInjection);

Jet create_gaussian_jet(double E_iso, double Gamma0, double theta_c, Profile2d inject = noInjection);

Profile2d create_iso_power_law_injection(double L0, double t0, double t_wait, double q);
#endif