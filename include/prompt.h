//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

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
PromptPhotonsGrid createPromptPhotonsGrid(size_t phi_size, size_t theta_size, size_t t_size);

class CoastingShock {
   public:
    CoastingShock(size_t phi_size, size_t theta_size, size_t t_size);
    CoastingShock() = delete;

    MeshGrid3d r;        // radius
    MeshGrid3d theta;    // theta for jet spreading
    MeshGrid3d Gamma;    // relative lorentz factor between down stream and up stream
    MeshGrid3d epsilon;  // relative energy per solid angle

    auto shape() const { return std::make_tuple(phi_size, theta_size, t_size); }  // Returns grid dimensions

   private:
    size_t const phi_size{0};    // Number of grid points in phi direction
    size_t const theta_size{0};  // Number of grid points in theta direction
    size_t const t_size{0};      // Number of grid points in time direction
};

template <typename Ejecta>
CoastingShock genCoastingShock(Coord const& coord, Ejecta const& jet) {
    auto [phi_size, theta_size, t_size] = coord.shape();
    CoastingShock shock(1, theta_size, t_size);

    for (size_t j = 0; j < theta_size; ++j) {
        Real Gamma = jet.Gamma0(coord.phi[0], coord.theta[j]);
        Real beta = gammaTobeta(Gamma);
        Real epsilon = jet.dE0dOmega(coord.phi[0], coord.theta[j]);
        for (size_t k = 0; k < t_size; ++k) {
            shock.Gamma[0][j][k] = Gamma;
            shock.epsilon[0][j][k] = epsilon;
            shock.r[0][j][k] = (beta * con::c) / std::abs(1 - beta) * coord.t[k];
            shock.theta[0][j][k] = coord.theta[j];
        }
    }

    return shock;
}

PromptPhotonsGrid genPromptPhotons(CoastingShock const& shock, Real R0, Real nu_0, Real alpha);
