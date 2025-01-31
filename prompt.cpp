#include "prompt.h"

PromptPhotonsMesh genPromptPhotons(Coord const& coord, Jet const& jet, double R0, double nu_0, double alpha) {
    PromptPhotonsMesh ph = PromptPhotonsMesh(coord.theta.size(), PromptPhotonsArray(coord.r.size()));

    double Gamma_c = jet.Gamma0_profile(0);
    double beta_c = gammaTobeta(Gamma_c);
    for (size_t j = 0; j < coord.theta.size(); ++j) {
        double theta = coord.theta[j];
        double Gamma = jet.Gamma0_profile(theta);
        double beta = gammaTobeta(Gamma);
        double R = R0 / (beta_c) * (beta);
        for (size_t k = 0; k < coord.r.size(); ++k) {
            if (coord.r_b[k + 1] > R && coord.r_b[k] < R) {
                double dOmega = std::fabs(std::cos(coord.theta_b[j + 1]) - std::cos(coord.theta_b[j])) * 2 * con::pi;
                ph[j][k].E_nu_peak = jet.dEdOmega(theta, 0) / Gamma * dOmega;

            } else {
                ph[j][k].E_nu_peak = 0;
            }
            ph[j][k].nu_0 = nu_0;
            ph[j][k].alpha = alpha;
        }
    }
    return ph;
}

double PromptPhotons::I_nu(double nu) const { return E_nu_peak * std::pow(nu / nu_0, -alpha); }