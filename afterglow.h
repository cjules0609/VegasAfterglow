//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

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
class VegasAfterglow {
   public:
    // Constructor: accepts the medium, jet, energy fractions, observation parameters, and observation times.
    VegasAfterglow(Coord const& coord, Medium const& medium, TophatJet const& jet, double eps_e, double eps_B,
                   double theta_obs, double lumi_dist, double z, Array const& t_obs)
        : medium_(medium),
          jet_(jet),
          eps_e_(eps_e),
          eps_B_(eps_B),
          theta_obs_(theta_obs),
          lumi_dist_(lumi_dist),
          z_(z),
          t_obs_(t_obs),
          coord_(coord),
          observer_(Observer(coord_)) {}

    VegasAfterglow() = delete;

    // Runs the simulation: generates shocks, electron and photon grids, applies electron cooling, and updates the
    // observer.
    void run() {
        // Generate the forward shock using a 3D interface.
        shock_ = genForwardShock3D(coord_, medium_, jet_, inject::none, eps_e_, eps_B_);
        // Generate synchrotron electron grid; assume xi = 1 for now and that jet_ provides parameter p.
        electrons_ = genSynElectrons(shock_, jet_.p, 1);
        // Generate synchrotron photon grid.
        photons_ = genSynPhotons(shock_, electrons_);
        // Apply electron cooling in the Thomson regime (could also use eCoolingKleinNishina).
        eCoolingThomson(electrons_, photons_, shock_);
        // Update the observer with the current shock parameters.
        observer_.observe(shock_, theta_obs_, lumi_dist_, z_);
    }

    // Outputs simulation results to files.
    void outputResults(std::string const& base_filename) {
        // Output the shock and coordinate grid using your existing output functions.
        output(shock_, base_filename + "_shock");
        output(coord_, base_filename + "_coord");
        // Example: compute and output specific flux for a chosen observed frequency (nu_obs).
        double nu_obs = 1e9;  // Example frequency [Hz]
        Array flux = observer_.specificFlux(t_obs_, nu_obs, photons_);
        // Use the unified output interface (for 1D arrays) to write the flux.
        // unifiedOutput(flux, base_filename + "_flux", 1);
    }

   private:
    Medium medium_;
    TophatJet jet_;
    double eps_e_;
    double eps_B_;
    double theta_obs_;
    double lumi_dist_;
    double z_;
    Array t_obs_;
    Coord coord_;
    Shock shock_;
    SynElectronGrid electrons_;
    SynPhotonGrid photons_;
    Observer observer_;  // Observer now constructed with the coordinate grid.
};*/
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