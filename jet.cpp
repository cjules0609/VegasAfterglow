
#include "jet.h"

#include <cmath>
#include <iostream>

#include "macros.h"
#include "physics.h"
#include "utilities.h"

auto LiangGhirlanda2010(auto energy_func, double e_max, double gamma_max, double idx) {
    return [=](double phi, double theta, double t = 0) -> double {
        double e = energy_func(phi, theta, 0);
        double u = std::pow(e / e_max, idx) * gamma_max;
        return std::sqrt(1 + u * u);
    };
}
Ejecta tophatJet(double theta_c, double E_iso, double Gamma0) {
    Ejecta jet;
    jet.dEdOmega = math::tophat(theta_c, E_iso / (4 * con::pi));
    jet.Gamma0 = math::tophat(theta_c, Gamma0);
    return jet;
}

Ejecta gaussianJet(double theta_c, double E_iso, double Gamma0, double idx) {
    Ejecta jet;
    jet.dEdOmega = math::gaussian(theta_c, E_iso / (4 * con::pi));
    jet.Gamma0 = LiangGhirlanda2010(jet.dEdOmega, E_iso / (4 * con::pi), Gamma0, idx);
    return jet;
}

Ejecta powerLawJet(double theta_c, double E_iso, double Gamma0, double k, double idx) {
    Ejecta jet;
    jet.dEdOmega = math::powerLaw(theta_c, E_iso / (4 * con::pi), k);
    jet.Gamma0 = LiangGhirlanda2010(jet.dEdOmega, E_iso / (4 * con::pi), Gamma0, idx);
    return jet;
}

std::tuple<double, double> findRadiusRange(double t_min, double t_max, double z, Ejecta const& jet) {
    auto theta = linspace(0, con::pi / 2, 50);
    double r_min = 0;
    double r_max = 0;

    /*double Gamma_min = jet.Gamma0_profile(0);

    double theta_max = 0;
    for (auto th : theta) {
        double G0 = jet.Gamma0_profile(th);
        if (G0 > con::Gamma_cut && G0 < Gamma_min) {
            Gamma_min = G0;
            theta_max = th;
        }
    }
    double beta_min = gammaTobeta(Gamma_min);

    r_min = t_min * con::c / ((1 + z) * (1 / beta_min - std::cos(theta_max)));*/

    r_min = t_min * con::c / (1 + z);

    double beta_max = gammaTobeta(jet.Gamma0(0, 0, 0));
    r_max = t_max * con::c / ((1 + z) * (1 / beta_max - 1)) / jet.Gamma0(0, 0, 0);

    return {r_min, r_max};
}