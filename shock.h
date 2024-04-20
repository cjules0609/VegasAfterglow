#ifndef _FSDYNAMICS_
#define _FSDYNAMICS_

#include <tuple>

#include "jet.h"
#include "medium.h"
#include "mesh.h"
#include "shock.h"

class Shock {
   public:
    Shock(Coord const& coord, double eps_e, double eps_B, double xi, double zeta);
    MeshGrid t_com;
    MeshGrid Gamma;
    MeshGrid B;
    MeshGrid D_com;  // co-moving shock width
    MeshGrid n_p;
    double eps_e;
    double eps_B;
    double xi;
    double zeta;
};

class BlastWaveEqn {
   public:
    BlastWaveEqn(Medium const& medium, Jet const& jet, double theta_lo, double theta_hi);

    void operator()(Array const& y, Array& dydr, double r);

    Medium medium;
    Jet jet;

   private:
    inline double dGammadr(double r, double Gamma, double u, double t_eng);

    inline double dUdr(double r, double Gamma, double u, double t_eng);

    inline double dtdr_eng(double Gamma);

    inline double dtdr_com(double Gamma);  // co-moving time

    inline double dxdr_com(double r, double Gamma, double D_com, double t_eng);  // co-moving reverse shock width

    inline double dDdr_com(double Gamma);  // co-moving shell width

   private:
    double theta_lo;
    double theta_hi;
    double theta;
    double dOmega;
};

std::pair<Shock, Shock> gen_shocks(Coord const& coord, Jet const& blast, Medium const& medium);
#endif