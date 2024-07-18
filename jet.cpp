
#include "jet.h"

#include <cmath>

#include "macros.h"
#include "utilities.h"

double E_iso2Gamma0(double E_iso, double gamma_max, double E) {
    double u = pow(E / E_iso, 0.25) * gamma_max;
    double gamma = sqrt(1 + u * u);
    return gamma;
}

Jet create_isotropic_jet(double E_iso, double Gamma0, double sigma0, double t_eng, Profile2d inject) {
    auto dEdOmega = [=](double theta, double t_lab) { return E_iso / (4 * con::pi) + inject(theta, t_lab); };

    auto Gamma = [=](double theta) { return Gamma0; };

    auto sigma = [=](double theta) { return sigma0; };

    return Jet{t_eng, dEdOmega, Gamma, sigma};
}

Jet create_tophat_jet(double E_iso, double Gamma0, double sigma0, double theta_c, double t_eng, Profile2d inject) {
    auto dEdOmega = [=](double theta, double t_lab) {
        return (theta < theta_c ? (E_iso / (4 * con::pi)) : 0) + inject(theta, t_lab);
    };

    auto Gamma = [=](double theta) { return theta < theta_c ? Gamma0 : 1; };

    auto sigma = [=](double theta) { return theta < theta_c ? sigma0 : 0; };

    return Jet{t_eng, dEdOmega, Gamma, sigma};
}

Jet create_power_law_jet(double E_iso_on_axis, double Gamma0_on_axis, double sigma0, double theta_m, double k,
                         double t_eng, Profile2d inject) {
    double e0 = E_iso_on_axis / (4 * con::pi);
    auto dEdOmega = [=](double theta, double t_lab) {
        return (theta < theta_m ? e0 : e0 * std::pow(theta / theta_m, -k)) + inject(theta, t_lab);
    };

    auto Gamma = [=](double theta) { return E_iso2Gamma0(e0, Gamma0_on_axis, dEdOmega(theta, 0)); };

    auto sigma = [=](double theta) { return sigma0; };

    return Jet{t_eng, dEdOmega, Gamma, sigma};
}

Jet create_gaussian_jet(double E_iso_on_axis, double Gamma0_on_axis, double sigma0, double theta_c, double t_eng,
                        Profile2d inject) {
    double e0 = E_iso_on_axis / (4 * con::pi);
    auto dEdOmega = [=](double theta, double t_lab) {
        return e0 * std::exp(-theta * theta / (2 * theta_c * theta_c)) + inject(theta, t_lab);
    };

    auto Gamma = [=](double theta) { return E_iso2Gamma0(e0, Gamma0_on_axis, dEdOmega(theta, 0)); };
    // auto Gamma = [=](double theta) { return Gamma0_on_axis; };

    auto sigma = [=](double theta) { return sigma0; };

    return Jet{t_eng, dEdOmega, Gamma, sigma};
}

Profile2d create_iso_power_law_injection(double L0, double t0, double t_wait, double q) {
    double t0q = std::pow(t0, q);
    double tw1q = std::pow(t0, 1 - q);
    double Omega = 4 * con::pi;
    if (std::fabs(q - 1) > 1e-6) {
        return [=](double theta, double t_lab) {
            return step_func(t_lab - t_wait) * L0 * t0q / (1 - q) * (std::pow(t_lab + t0 - t_wait, 1 - q) - tw1q) /
                   Omega;
        };
    } else {
        return [=](double theta, double t_lab) {
            return step_func(t_lab - t_wait) * L0 * t0 * std::log((t_lab + t0 - t_wait) / t_wait) / Omega;
        };
    }
}