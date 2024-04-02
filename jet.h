#ifndef _JET_
#define _JET_

#include "mesh.h"

class Jet {
   public:
    Profile2d dEdOmega;
    Profile Gamma0;
};

static Profile2d noInjection = Profile2d([](double theta, double t_lab) { return 0; });

Jet createIsotropicJet(double E_iso, double Gamma0, Profile2d inject = noInjection);

Jet createTopHatJet(double E_iso, double Gamma0, double theta_c, Profile2d inject = noInjection);

Jet createPowerLawJet(double E_iso, double Gamma0, double theta_m, double k, Profile2d inject = noInjection);

Jet createGaussianJet(double E_iso, double Gamma0, double theta_c, Profile2d inject = noInjection);

Profile2d createIsoPowerLawInjection(double L0, double t0, double t_wait, double q);
#endif