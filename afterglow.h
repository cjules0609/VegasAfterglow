#ifndef _AFTERGLOW_
#define _AFTERGLOW_

#include "BS_thread_pool.hpp"
// inline BS::thread_pool th_pool;

#include "IO.h"
#include "inverse-compton.h"
#include "jet.h"
#include "macros.h"
#include "medium.h"
#include "mesh.h"
#include "observer.h"
#include "physics.h"
#include "prompt.h"
#include "shock.h"
#include "synchrotron.h"
#include "utilities.h"

/*
template <typename... Photons>
MeshGrid co_moving_spectrum(size_t spectrum_resol, double nu_min, double nu_max, Photons const&... ph) {
    Array nu_bin = logspace(nu_min, nu_max, spectrum_resol + 1);
    Array nu_c = boundaryToCenter(nu_bin);
    MeshGrid L_nu = createGrid(2, spectrum_resol, 0);
    for (size_t i = 0; i < L_nu[0].size(); ++i) {
        L_nu[0][i] = nu_c[i];
        L_nu[1][i] = (ph.L_nu(nu_c[i]) + ...);
    }
    return L_nu;
}

template <typename Photons1Array, typename... Photons2Array>
MeshGrid co_moving_spectrums(size_t spectrum_resol, double nu_min, double nu_max, Photons1Array const& ph0,
                             Photons2Array const&... ph) {
    Array nu_bin = logspace(nu_min, nu_max, spectrum_resol + 1);
    Array nu_c = boundaryToCenter(nu_bin);
    MeshGrid L_nu = createGrid(1 + ph0.size(), spectrum_resol, 0);
    L_nu[0] = nu_c;
    for (size_t i = 0; i < ph0.size(); ++i) {
        for (size_t j = 0; j < L_nu[0].size(); ++j) {
            L_nu[i + 1][j] = ph0[i].L_nu(nu_c[j]);
            if constexpr (sizeof...(ph) > 0) {
                L_nu[i + 1][j] += (ph[i].L_nu(nu_c[j]) + ...);
            }
        }
    }
    return L_nu;
}

template <typename... Electrons>
MeshGrid co_moving_e_spectrum(size_t spectrum_resol, double gamma_min, double gamma_max, Electrons const&... e) {
    Array gamma_bin = logspace(gamma_min, gamma_max, spectrum_resol + 1);
    Array gamma = boundaryToCenter(gamma_bin);
    MeshGrid n_gamma = createGrid(2, spectrum_resol, 0);
    for (size_t i = 0; i < n_gamma[0].size(); ++i) {
        n_gamma[0][i] = gamma[i];
        n_gamma[1][i] = (e.N(gamma[i]) + ...);
    }
    return n_gamma;
}

template <typename Electrons1Array, typename... Electrons2Array>
MeshGrid co_moving_e_spectrums(size_t spectrum_resol, double gamma_min, double gamma_max, Electrons1Array const& e0,
                               Electrons2Array const&... e) {
    Array gamma_bin = logspace(gamma_min, gamma_max, spectrum_resol + 1);
    Array gamma = boundaryToCenter(gamma_bin);
    MeshGrid n_gamma = createGrid(1 + e0.size(), spectrum_resol, 0);
    n_gamma[0] = gamma;
    for (size_t i = 0; i < e0.size(); ++i) {
        for (size_t j = 0; j < n_gamma[0].size(); ++j) {
            n_gamma[i + 1][j] = e0[i].N(gamma[j]);
            if constexpr (sizeof...(e) > 0) {
                n_gamma[i + 1][j] += (e[i].N(gamma[j]) + ...);
            }
        }
    }
    return n_gamma;
}*/
#endif