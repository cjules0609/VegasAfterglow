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

template <typename Photons>
MeshGrid spectrum(Photons const& ph, double nu_min, double nu_max, size_t data_points = 100) {
    Array nu_bin = logspace(nu_min, nu_max, data_points + 1);
    Array nu_c = boundary2center(nu_bin);
    MeshGrid L_nu = create_grid(2, data_points, 0);
    for (size_t i = 0; i < L_nu[0].size(); ++i) {
        L_nu[0][i] = nu_c[i];
        L_nu[1][i] = ph.j_nu(nu_c[i]);
    }
    return L_nu;
}

template <typename PhotonsArray>
MeshGrid full_spectrum(PhotonsArray const& ph, double nu_min, double nu_max, size_t data_points = 100) {
    Array nu_bin = logspace(nu_min, nu_max, data_points + 1);
    Array nu_c = boundary2center(nu_bin);
    MeshGrid L_nu = create_grid(1 + ph.size(), data_points, 0);
    L_nu[0] = nu_c;
    for (size_t i = 0; i < ph.size(); ++i) {
        for (size_t j = 0; j < L_nu[0].size(); ++j) {
            L_nu[i + 1][j] = ph[i].j_nu(nu_c[j]);
        }
    }
    return L_nu;
}
#endif