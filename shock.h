#ifndef _FSDYNAMICS_
#define _FSDYNAMICS_

#include <tuple>

#include "jet.h"
#include "medium.h"
#include "mesh.h"
#include "physics.h"
#include "shock.h"
class Shock {
   public:
    Shock(Coord const& coord, double eps_e, double eps_B, double p, double xi = 1, double zeta = 1);
    MeshGrid t_com;
    MeshGrid dt_com;
    MeshGrid Gamma;      // bulk Lorentz factor
    MeshGrid e_th;       // internal energy density
    MeshGrid B;          // comoving magnetic field
    MeshGrid width_eff;  // comoving frame effective shock width that contains all shocked electrons
    MeshGrid n_p;        // n2
    double eps_e;
    double eps_B;
    double p;
    double xi;
    double zeta;
};

class ForwardShockEqn {
   public:
    ForwardShockEqn(Medium const& medium, Jet const& jet, double theta, double eps_e);

    void operator()(Array const& y, Array& dydr, double r);

    Medium const& medium;
    Jet const& jet;
    double theta;
    double eps_e;
    double gamma4;
    double spreading_factor;

   private:
    inline double dGammadr(double r, double Gamma, double u, double t_eng);

    inline double dUdr(double r, double Gamma, double u, double t_eng);

   private:
    double dM0dOmega;
};

class FRShockEqn {
   public:
    FRShockEqn(Medium const& medium, Jet const& jet, double theta);

    void operator()(Array const& y, Array& dydr, double r);

    Medium const& medium;
    Jet const& jet;
    double sigma;
    double theta;
    double gamma4;

   private:
    inline double dN3dr_per_Omega(double r, double n1, double n4, double gamma3);
};

// void solve_shocks(Coord const& coord, Jet const& jet, Medium const& medium, Shock& f_shock, Shock& r_shock);
// void solve_shocks(Coord const& coord, Jet const& jet, Medium const& medium, Shock& f_shock);
Shock gen_forward_shock(Coord const& coord, Jet const& jet, Medium const& medium, double eps_e, double eps_B, double p,
                        double xi = 1, double zeta = 1);

std::pair<Shock, Shock> gen_fr_shocks(Coord const& coord, Jet const& jet, Medium const& medium, double eps_e,
                                      double eps_B, double p, double xi = 1, double zeta = 1);
#endif