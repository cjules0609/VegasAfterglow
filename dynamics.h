#ifndef _DYNAMICS_
#define _DYNAMICS_

#include <tuple>

#include "jet.h"
#include "medium.h"
class DynamicsEqn {
   public:
    DynamicsEqn(Medium medium, Jet blast, double theta_lo, double theta_hi);

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

//return t_co-moving, Gamma
std::tuple<MeshGrid, MeshGrid> solve_dynamics(Coord const& coord, Jet const& blast, Medium const& medium);
#endif