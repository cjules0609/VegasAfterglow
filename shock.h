#ifndef _FSDYNAMICS_
#define _FSDYNAMICS_

#include <tuple>

#include "jet.h"
#include "medium.h"
#include "mesh.h"
#include "shock.h"
#include "physics.h"
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

class BlastWaveEqn {
   public:
    BlastWaveEqn(Medium const& medium, Jet const& jet, double theta_lo, double theta_hi, double eps_e,
                 bool reverse_shock);

    void operator()(Array const& y, Array& dydr, double r);

    Medium medium;
    Jet jet;
    UDownStr u_down;
    double sigma;
    double dOmega;
    double theta;
    double eps_e;
    bool reverse_shock;

   private:
    inline double dGammadr(double r, double Gamma, double u, double t_eng);

    inline double dUdr(double r, double Gamma, double u, double t_eng);

    inline double dtdr_eng(double Gamma);

    inline double dtdr_com(double Gamma);  // co-moving time

    inline double dDdr_RS(double r, double Gamma, double D_com, double t_eng);  // lab frame reverse shock width

    inline double dDdr_FS(double r, double Gamma);  // lab frame shell width

   private:
    double theta_lo;
    double theta_hi;
};

void solve_shocks(Coord const& coord, Jet const& jet, Medium const& medium, Shock& f_shock, Shock& r_shock);
void solve_shocks(Coord const& coord, Jet const& jet, Medium const& medium, Shock& f_shock);
// std::pair<Shock, Shock> gen_shocks(Coord const& coord, Jet const& blast, Medium const& medium);
#endif