#ifndef _RELATIVITY_H_
#define _RELATIVITY_H_
#include <cmath>

#include "macros.h"

inline double gamma_to_beta(double gamma) { return sqrt(1 - 1 / (gamma * gamma)); }

inline double adiabatic_index(double gamma) { return (4 * gamma + 1) / (3 * gamma); }

double z_to_luminosity_distance(double z);

double luminosity_distance_to_z(double L);

double dec_radius(double E_iso, double n_ism, double Gamma0, double engine_dura);

double thin_shell_dec_radius(double E_iso, double n_ism, double Gamma0);

double thick_shell_dec_radius(double E_iso, double n_ism, double Gamma0, double engine_dura);

double shell_spreading_radius(double Gamma0, double engine_dura);

double RS_transition_radius(double E_iso, double n_ism, double Gamma0, double engine_dura);

inline double RS_crossing_radius(double E_iso, double n_ism, double Gamma0, double engine_dura) {
    return thick_shell_dec_radius(E_iso, n_ism, Gamma0, engine_dura);
}

inline double Sedov_length(double E_iso, double n_ism) {
    return pow(E_iso / (4 * con::pi / 3 * n_ism * con::mp * con::c2), 1.0 / 3);
}
#endif  // _RELATIVITY_H_
