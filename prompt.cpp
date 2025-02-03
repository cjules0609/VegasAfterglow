#include "prompt.h"

PromptPhotonsGrid createPromptPhotonsGrid(size_t phi_size, size_t theta_size, size_t r_size) {
    return PromptPhotonsGrid(boost::extents[phi_size][theta_size][r_size]);
}
PromptPhotonsGrid genPromptPhotons(Coord const& coord, Ejecta const& jet, double R0, double nu_0, double alpha) {
    auto [phi_size, theta_size, r_size] = coord.shape();

    PromptPhotonsGrid ph = createPromptPhotonsGrid(phi_size, theta_size, r_size);

    double Gamma_c = jet.Gamma0(0, 0, 0);
    double beta_c = gammaTobeta(Gamma_c);

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            double theta = coord.theta[j];
            double Gamma = jet.Gamma0(0, theta, 0);
            double beta = gammaTobeta(Gamma);
            double R = R0 / (beta_c) * (beta);
            for (size_t k = 0; k < r_size - 1; ++k) {
                if (coord.r[k + 1] > R && coord.r[k] < R) {
                    ph[i][j][k].E_nu_peak = jet.dEdOmega(0, theta, 0) / Gamma;

                } else {
                    ph[i][j][k].E_nu_peak = 0;
                }
                ph[i][j][k].nu_0 = nu_0;
                ph[i][j][k].alpha = alpha;
            }
        }
    }
    return ph;
}

double PromptPhotons::I_nu(double nu) const { return E_nu_peak * std::pow(nu / nu_0, -alpha); }