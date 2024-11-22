#ifndef _JET_
#define _JET_

#include "macros.h"
#include "mesh.h"

struct Injection {
    Profile2d dLdOmega;
    Profile2d dEdOmega;
};
class Jet {
   public:
    double duration{1 * con::sec};
    double theta_c0{0};
    bool spreading{false};
    Profile3d dEdOmega_spread;
    Profile2d dEdOmega;
    Profile dE0dOmega;
    Profile Gamma0_profile;
    Profile sigma_profile;
    Injection inj;
};

static Injection noInjection = {Profile2d([](double theta, double t_lab) { return 0; }),
                                Profile2d([](double theta, double t_lab) { return 0; })};
class IsoJet : public Jet {
   public:
    IsoJet(double E_iso, double Gamma0, double duration = 1 * con::sec, double sigma0 = 0,
           Injection inject = noInjection);
};

class TophatJet : public Jet {
   public:
    TophatJet(double theta_c0, double E_iso, double Gamma0, double duration = 1 * con::sec, double sigma0 = 0,
              Injection inject = noInjection);
};

class GaussianJet : public Jet {
   public:
    GaussianJet(double theta_c0, double E_iso, double Gamma0, double Gamma_idx = 1, double duration = 1 * con::sec,
                double sigma0 = 0, Injection inject = noInjection);
};

class PowerLawJet : public Jet {
   public:
    PowerLawJet(double theta_c0, double k, double E_iso, double Gamma0, double Gamma_idx = 1,
                double duration = 1 * con::sec, double sigma0 = 0, Injection inject = noInjection);
};

Injection create_iso_const_injection(double L0, double t0);
Injection create_iso_power_law_injection(double L0, double t0, double q);
#endif