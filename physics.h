//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

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

template <typename Jet>
double jetEdge(Jet const& jet, double gamma_cut) {
    if (jet.Gamma0(0, con::pi / 2, 0) > gamma_cut) {
        return con::pi / 2;
    }
    double low = 0;
    double hi = con::pi / 2;
    double eps = 1e-6;
    for (; hi - low > eps;) {
        double mid = 0.5 * (low + hi);
        if (jet.Gamma0(0, mid, 0) > gamma_cut) {
            low = mid;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (low + hi);
}

template <typename Jet, typename Injector>
Coord adaptiveGrid(Medium const& medium, Jet const& jet, Injector const& inj, Array const& t_obs, double theta_max,
                   size_t phi_num = 32, size_t theta_num = 32, size_t r_num = 32) {
    double t_max = *std::max_element(t_obs.begin(), t_obs.end());
    double t_min = *std::min_element(t_obs.begin(), t_obs.end());
    auto [r_min, r_max] = findRadiusRange(medium, jet, inj, t_min, t_max);
    Array r = logspace(r_min, r_max, r_num);
    double jet_edge = jetEdge(jet, con::Gamma_cut);
    Array theta = uniform_cos(0, std::min(jet_edge, theta_max), theta_num);
    Array phi = linspace(0, 2 * con::pi, phi_num);
    Coord coord{phi, theta, r};
    coord.t_min = t_min;
    coord.t_max = t_max;
    return coord;
}
#endif  // _RELATIVITY_H_
