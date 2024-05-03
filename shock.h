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
    MeshGrid t_com_b;
    MeshGrid Gamma;  // bulk Lorentz factor
    MeshGrid e_th;   // internal energy density
    MeshGrid B;      // comoving magnetic field
    MeshGrid width;  // comoving frame shock width
    MeshGrid n_p;    // n2
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
    double dOmega;
    double theta;

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

std::pair<Shock, Shock> gen_shocks(Coord const& coord, Jet const& blast, Medium const& medium);
#endif