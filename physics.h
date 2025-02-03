#ifndef _RELATIVITY_H_
#define _RELATIVITY_H_
#include <cmath>

#include "jet.h"
#include "macros.h"
#include "medium.h"
#include "mesh.h"

double zToLuminosityDistance(double z);
double luminosityDistanceToz(double L);
double decRadius(double E_iso, double n_ism, double Gamma0, double engine_dura);
double thinShellDecRadius(double E_iso, double n_ism, double Gamma0);
double thickShellDecRadius(double E_iso, double n_ism, double Gamma0, double engine_dura);
double shellSpreadingRadius(double Gamma0, double engine_dura);
double RSTransitionRadius(double E_iso, double n_ism, double Gamma0, double engine_dura);
inline double gammaTobeta(double gamma) { return std::sqrt(1 - 1 / (gamma * gamma)); }
inline double adiabaticIndex(double gamma) { return (4 * gamma + 1) / (3 * gamma); }
inline double RShockCrossingRadius(double E_iso, double n_ism, double Gamma0, double engine_dura) {
    return thickShellDecRadius(E_iso, n_ism, Gamma0, engine_dura);
}
inline double SedovLength(double E_iso, double n_ism) {
    return std::cbrt(E_iso / (4 * con::pi / 3 * n_ism * con::mp * con::c2));
}
double jetEdge(TernaryFunc const& gamma, double gamma_cut);

std::tuple<double, double> findRadiusRange(Medium const& medium, Ejecta const& jet, Ejecta const& inj, double t_min,
                                           double t_max, double z = 0);

Coord adaptiveGrid(Medium const& medium, Ejecta const& jet, Ejecta const& inj, Array const& t_obs, double theta_max,
                   size_t phi_num = 32, size_t theta_num = 32, size_t r_num = 32);
#endif  // _RELATIVITY_H_
