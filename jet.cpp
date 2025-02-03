
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
    return Ejecta(math::tophat(theta_c, E_iso / (4 * con::pi)), math::tophat(theta_c, Gamma0));
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
