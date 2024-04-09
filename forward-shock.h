#ifndef _FSDYNAMICS_
#define _FSDYNAMICS_

#include <tuple>

#include "jet.h"
#include "medium.h"
#include "mesh.h"
#include "shock.h"

class ForwardShockEqn {
   public:
    ForwardShockEqn(Medium const& medium, Jet const& jet, double theta_lo, double theta_hi);

    void operator()(Array const& y, Array& dydr, double r);

    Medium medium;
    Jet blast;

   private:
    inline double dGammadr(double r, double Gamma, double u, double t_eng);

    inline double dUdr(double r, double Gamma, double u, double t_eng);

    inline double dtdr_eng(double Gamma);

    inline double dtdr_com(double Gamma);  // co-moving time

   private:
    double theta_lo;
    double theta_hi;
    double theta;
    double dOmega;
};

Shock gen_forward_shock(Coord const& coord, Jet const& blast, Medium const& medium);
#endif