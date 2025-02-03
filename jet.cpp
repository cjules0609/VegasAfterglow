//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

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

double LiangGhirlanda2010_(double e, double e_max, double gamma_max, double idx) {
    double u = std::pow(e / e_max, idx) * gamma_max;
    return std::sqrt(1 + u * u);
}

double TophatJet::dEdOmega(double phi, double theta, double t) const { return theta < theta_c_ ? dEdOmega_ : 0; }

double TophatJet::Gamma0(double phi, double theta, double t) const { return theta < theta_c_ ? Gamma0_ : 0; }

double GaussianJet::dEdOmega(double phi, double theta, double t) const {
    return dEdOmega_ * fastExp(-theta * theta / (2 * theta_c_ * theta_c_));
}

double GaussianJet::Gamma0(double phi, double theta, double t) const {
    return LiangGhirlanda2010_(dEdOmega(phi, theta, t), dEdOmega_, Gamma0_, idx_);
}

double PowerLawJet::dEdOmega(double phi, double theta, double t) const {
    return dEdOmega_ * fastPow(theta / theta_c_, -k_);
}

double PowerLawJet::Gamma0(double phi, double theta, double t) const {
    return LiangGhirlanda2010_(dEdOmega(phi, theta, t), dEdOmega_, Gamma0_, idx_);
}
