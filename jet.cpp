
#include "jet.h"

#include <cmath>
#include <iostream>

#include "macros.h"
#include "utilities.h"

double E2Gamma0(double gamma_max, double e_iso, double e, double idex) {
    double u = pow(e / e_iso, idex) * gamma_max;
    double gamma = sqrt(1 + u * u);
    return gamma;
}

/*void Jet::jet_spread(double Gamma, double cs, double r, double dr) {
    if (this->theta_c0 < con::pi / 2) {
        this->theta_c0 += cs / con::c * dr / (r * Gamma * Gamma * theta_c0);
    }
}*/

IsoJet::IsoJet(double E_iso, double Gamma0, double duration, double sigma0, Injection inject) {
    this->duration = duration;

    this->theta_c0 = con::pi / 2;

    this->spreading = false;

    this->inj = inject;

    this->dEdOmega = [=](double theta, double t_lab) { return E_iso / (4 * con::pi) + inj.dEdOmega(theta, t_lab); };

    this->dE0dOmega = [=](double theta) { return E_iso / (4 * con::pi); };

    this->Gamma0_profile = [=](double theta) { return Gamma0; };

    this->sigma_profile = [=](double theta) { return sigma0; };

    this->dEdOmega_spread = [=](double theta, double theta_c, double t_lab) { return this->dEdOmega(theta, t_lab); };
}

TophatJet::TophatJet(double theta_c0, double E_iso, double Gamma0, double duration, double sigma0, Injection inject) {
    this->duration = duration;

    this->theta_c0 = theta_c0;

    this->spreading = false;

    double e_iso = E_iso / (4 * con::pi);

    this->dE0dOmega = [=](double theta) { return theta < theta_c0 ? e_iso : 0; };

    this->Gamma0_profile = [=](double theta) { return theta < theta_c0 ? Gamma0 : 1; };

    this->sigma_profile = [=](double theta) { return theta < theta_c0 ? sigma0 : 0; };

    this->inj.dLdOmega = [=](double theta, double t_lab) {
        return theta < theta_c0 ? inject.dLdOmega(theta, t_lab) : 0;
    };

    this->inj.dEdOmega = [=](double theta, double t_lab) {
        return theta < theta_c0 ? inject.dEdOmega(theta, t_lab) : 0;
    };

    this->dEdOmega = [=](double theta, double t_lab) { return dE0dOmega(theta) + inj.dEdOmega(theta, t_lab); };

    this->dEdOmega_spread = [=](double theta, double theta_c, double t_lab) {
        return theta < theta_c ? (1 - cos(theta_c0)) / (1 - cos(theta_c)) * (e_iso + inj.dEdOmega(theta, t_lab)) : 0;
    };
}

GaussianJet::GaussianJet(double theta_c0, double E_iso, double Gamma0, double Gamma_idx, double duration, double sigma0,
                         Injection inject) {
    this->duration = duration;

    this->theta_c0 = theta_c0;

    this->spreading = false;

    double e_iso = E_iso / (4 * con::pi);

    this->dE0dOmega = [=](double theta) { return e_iso * std::exp(-theta * theta / (2 * theta_c0 * theta_c0)); };

    this->Gamma0_profile = [=](double theta) { return E2Gamma0(Gamma0, e_iso, dE0dOmega(theta), Gamma_idx); };

    this->sigma_profile = [=](double theta) { return sigma0; };

    this->inj.dLdOmega = [=](double theta, double t_lab) {
        return theta < 10 * theta_c0 ? inject.dLdOmega(theta, t_lab) : 0;
    };

    this->inj.dEdOmega = [=](double theta, double t_lab) {
        return theta < 10 * theta_c0 ? inject.dEdOmega(theta, t_lab) : 0;
    };

    this->dEdOmega = [=](double theta, double t_lab) { return dE0dOmega(theta) + inj.dEdOmega(theta, t_lab); };

    this->dEdOmega_spread = [=](double theta, double theta_c, double t_lab) {
        return (1 - cos(theta_c0)) / (1 - cos(theta_c)) *
               (e_iso * std::exp(-theta * theta / (2 * theta_c * theta_c)) + inj.dEdOmega(theta, t_lab));
    };
}

PowerLawJet::PowerLawJet(double theta_c0, double k, double E_iso, double Gamma0, double Gamma_idx, double duration,
                         double sigma0, Injection inject) {
    this->duration = duration;

    this->theta_c0 = theta_c0;

    this->spreading = false;

    double e_iso = E_iso / (4 * con::pi);

    this->dE0dOmega = [=](double theta) { return (theta < theta_c0 ? e_iso : e_iso * std::pow(theta / theta_c0, -k)); };

    this->Gamma0_profile = [=](double theta) { return E2Gamma0(Gamma0, e_iso, dE0dOmega(theta), Gamma_idx); };

    this->sigma_profile = [=](double theta) { return sigma0; };

    this->inj.dLdOmega = [=](double theta, double t_lab) {
        return theta < theta_c0 ? inject.dLdOmega(theta, t_lab) : 0;
    };

    this->inj.dEdOmega = [=](double theta, double t_lab) {
        return theta < theta_c0 ? inject.dEdOmega(theta, t_lab) : 0;
    };

    this->dEdOmega = [=](double theta, double t_lab) { return dE0dOmega(theta) + inj.dEdOmega(theta, t_lab); };

    this->dEdOmega_spread = [=](double theta, double theta_c, double t_lab) {
        return (1 - cos(theta_c0)) / (1 - cos(theta_c)) *
               ((theta < theta_c ? e_iso : e_iso * std::pow(theta / theta_c, -k)) + inj.dEdOmega(theta, t_lab));
    };
}

CosJet::CosJet(double theta_c0, double E_iso, double Gamma0, double Gamma_idx, double duration, double sigma0,
               Injection inject) {
    this->duration = duration;

    this->theta_c0 = theta_c0;

    this->spreading = false;

    double e_iso = E_iso / (4 * con::pi);

    this->dE0dOmega = [=](double theta) {
        return e_iso * std::fabs(std::exp(-theta * 30 / con::pi) * std::cos(theta / theta_c0));
    };

    this->Gamma0_profile = [=](double theta) { return E2Gamma0(Gamma0, e_iso, dE0dOmega(theta), Gamma_idx); };

    this->sigma_profile = [=](double theta) { return sigma0; };

    this->inj.dLdOmega = [=](double theta, double t_lab) {
        return theta < theta_c0 ? inject.dLdOmega(theta, t_lab) : 0;
    };

    this->inj.dEdOmega = [=](double theta, double t_lab) {
        return theta < theta_c0 ? inject.dEdOmega(theta, t_lab) : 0;
    };

    this->dEdOmega = [=](double theta, double t_lab) { return dE0dOmega(theta) + inj.dEdOmega(theta, t_lab); };

    this->dEdOmega_spread = [=](double theta, double theta_c, double t_lab) {
        return (1 - cos(theta_c0)) / (1 - cos(theta_c)) *
               (e_iso * std::fabs(std::cos(theta / theta_c)) + inj.dEdOmega(theta, t_lab));
    };
}

Injection create_iso_const_injection(double L0, double t0) {
    auto dLdOmega = [=](double theta, double t_lab) {
        if (t_lab < t0) {
            return L0 / (4 * con::pi);
        } else {
            return 0.0;
        }
    };
    auto dEdOmega = [=](double theta, double t_lab) {
        if (t_lab < t0) {
            return L0 * t_lab / (4 * con::pi);
        } else {
            return L0 * t0 / (4 * con::pi);
        }
    };
    return Injection{dLdOmega, dEdOmega};
}

Injection create_iso_power_law_injection(double L0, double t0, double q) {
    auto dLdOmega = [=](double theta, double t_lab) { return L0 * std::pow(1 + t_lab / t0, -q) / (4 * con::pi); };
    auto dEdOmega = [=](double theta, double t_lab) {
        if (fabs(q - 1) > 1e-6) {
            return L0 * t0 / (1 - q) * (std::pow(1 + t_lab / t0, 1 - q) - 1) / (4 * con::pi);
        } else {
            return L0 * t0 * log(1 + t_lab / t0) / (4 * con::pi);
        }
    };
    return Injection{dLdOmega, dEdOmega};
}
