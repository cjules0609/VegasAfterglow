#ifndef _AFTERGLOW_
#define _AFTERGLOW_

#include "IO.h"
#include "forward-shock.h"
#include "inverse-compton.h"
#include "jet.h"
#include "macros.h"
#include "medium.h"
#include "mesh.h"
#include "observer.h"
#include "synchrotron.h"
#include "utilities.h"

template <typename... Photons>
MeshGrid co_moving_spectrum(size_t spectrum_resol, double nu_min, double nu_max, Photons const&... ph) {
    Array nu_bin = logspace(nu_min, nu_max, spectrum_resol + 1);
    Array nu_c = boundary2center(nu_bin);
    MeshGrid j_nu = create_grid(2, spectrum_resol, 0);
    for (size_t i = 0; i < j_nu[0].size(); ++i) {
        j_nu[0][i] = nu_c[i];
        j_nu[1][i] = (ph.j_nu(nu_c[i]) + ...);
    }
    return j_nu;
}

template <typename Photons1Array, typename... Photons2Array>
MeshGrid co_moving_spectrums(size_t data_points, double nu_min, double nu_max, Photons1Array const& ph0,
                             Photons2Array const&... ph) {
    Array nu_bin = logspace(nu_min, nu_max, data_points + 1);
    Array nu_c = boundary2center(nu_bin);
    MeshGrid j_nu = create_grid(1 + ph0.size(), data_points, 0);
    j_nu[0] = nu_c;
    for (size_t i = 0; i < ph0.size(); ++i) {
        for (size_t j = 0; j < j_nu[0].size(); ++j) {
            j_nu[i + 1][j] = ph0[i].j_nu(nu_c[j]);
            if constexpr (sizeof...(ph) > 0) {
                j_nu[i + 1][j] += (ph[i].j_nu(nu_c[j]) + ...);
            }
        }
    }
    return j_nu;
}

template <typename Electrons1Array, typename... Electrons2Array>
MeshGrid co_moving_n_spectrums(size_t data_points, double gamma_min, double gamma_max, Electrons1Array const& e0,
                               Electrons2Array const&... e) {
    Array gamma_bin = logspace(gamma_min, gamma_max, data_points + 1);
    Array gamma = boundary2center(gamma_bin);
    MeshGrid n_gamma = create_grid(1 + e0.size(), data_points, 0);
    n_gamma[0] = gamma;
    for (size_t i = 0; i < e0.size(); ++i) {
        for (size_t j = 0; j < n_gamma[0].size(); ++j) {
            n_gamma[i + 1][j] = e0[i].n(gamma[j]);
            if constexpr (sizeof...(e) > 0) {
                n_gamma[i + 1][j] += (e[i].n(gamma[j]) + ...);
            }
        }
    }
    return n_gamma;
}
#endif