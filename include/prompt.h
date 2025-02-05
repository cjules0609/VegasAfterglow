//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _PROMPT_
#define _PROMPT_

#include "jet.h"
#include "mesh.h"
#include "physics.h"

struct PromptPhotons {
    Real E_nu_peak{0};
    Real nu_0{0};
    Real alpha{0};

    Real I_nu(Real nu) const;
};

using PromptPhotonsGrid = boost::multi_array<PromptPhotons, 3>;
PromptPhotonsGrid createPromptPhotonsGrid(size_t phi_size, size_t theta_size, size_t r_size);

template <typename Jet>
PromptPhotonsGrid genPromptPhotons(Coord const& coord, Jet const& jet, Real R0, Real nu_0, Real alpha) {
    auto [phi_size, theta_size, r_size] = coord.shape();

    PromptPhotonsGrid ph = createPromptPhotonsGrid(phi_size, theta_size, r_size);

    Real Gamma_c = jet.Gamma0(0, 0, 0);
    Real beta_c = gammaTobeta(Gamma_c);

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            Real theta = coord.theta[j];
            Real Gamma = jet.Gamma0(0, theta, 0);
            Real beta = gammaTobeta(Gamma);
            Real R = R0 / (beta_c) * (beta);
            for (size_t k = 0; k < r_size - 1; ++k) {
                if (coord.t[k + 1] > R && coord.t[k] < R) {
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

#endif