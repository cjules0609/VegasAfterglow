#ifndef _RELATIVITY_H_
#define _RELATIVITY_H_
#include <cmath>

#include "macros.h"
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
#endif  // _RELATIVITY_H_
