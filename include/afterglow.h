//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _AFTERGLOW_
#define _AFTERGLOW_

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
    VegasAfterglow(Coord const& coord, Medium const& medium, TophatJet const& jet, Real eps_e, Real eps_B,
                   Real theta_obs, Real lumi_dist, Real z, Array const& t_obs)
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
        Real nu_obs = 1e9;  // Example frequency [Hz]
        Array flux = observer_.specificFlux(t_obs_, nu_obs, photons_);
        // Use the unified output interface (for 1D arrays) to write the flux.
        // unifiedOutput(flux, base_filename + "_flux", 1);
    }

   private:
    Medium medium_;
    TophatJet jet_;
    Real eps_e_;
    Real eps_B_;
    Real theta_obs_;
    Real lumi_dist_;
    Real z_;
    Array t_obs_;
    Coord coord_;
    Shock shock_;
    SynElectronGrid electrons_;
    SynPhotonGrid photons_;
    Observer observer_;  // Observer now constructed with the coordinate grid.
};*/
/*
template <typename... Photons>
MeshGrid co_moving_spectrum(size_t spectrum_resol, Real nu_min, Real nu_max, Photons const&... ph) {
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
MeshGrid co_moving_spectrums(size_t spectrum_resol, Real nu_min, Real nu_max, Photons1Array const& ph0,
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
MeshGrid co_moving_e_spectrum(size_t spectrum_resol, Real gamma_min, Real gamma_max, Electrons const&... e) {
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
MeshGrid co_moving_e_spectrums(size_t spectrum_resol, Real gamma_min, Real gamma_max, Electrons1Array const& e0,
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

/********************************************************************************************************************
 * FUNCTION: find_r_max
 * DESCRIPTION: Finds the maximum radius (r_max) for the shock evolution by integrating the shock equations until
 *              the engine time exceeds t_max. Uses a dense output stepper.
 ********************************************************************************************************************/
template <typename ShockEqn>
Real find_r_max(ShockEqn& eqn, Real r_min, Real t_max) {
    using namespace boost::numeric::odeint;
    Real atol = 0, rtol = 1e-6, r0 = r_min;
    Real dr = r_min / 100;

    typename ShockEqn::State state;
    setForwardInit(eqn, state, r0);  // Initialize the state at r0

    auto stepper = make_dense_output(atol, rtol, runge_kutta_dopri5<typename ShockEqn::State>());
    stepper.initialize(state, r0, dr);
    // Integrate until the engine time in the state exceeds t_max
    for (; state[2] <= t_max;) {
        stepper.do_step(eqn);
        state = stepper.current_state();
    }
    return stepper.current_time() + stepper.current_time_step();
}

/********************************************************************************************************************
 * TEMPLATE HELP FUNCTIONS
 * DESCRIPTION: The following template functions are helper functions for the shock generation and evolution process.
 ********************************************************************************************************************/

// Determines the range of radii (r_min, r_max) for the shock evolution using medium and jet properties.
// The minimum radius is based on t_min and the maximum is found by solving the forward shock equation on-axis.
template <typename Jet, typename Injector>
std::tuple<Real, Real> findRadiusRange(Medium const& medium, Jet const& jet, Injector const& inj, Real t_min,
                                       Real t_max, Real z = 0) {
    Real gamma0 = jet.Gamma0(0, 0, 0);  // On-axis Gamma
    Real gamma_min = (gamma0 - 1) / 100 + 1;
    Real beta_min = std::sqrt(1 - 1 / (gamma_min * gamma_min));
    Real r_min = beta_min * con::c * t_min / (1 + beta_min);

    // Find the on-axis r_max by solving for theta = 0, phi = 0
    auto eqn = ForwardShockEqn<Jet, Injector>(medium, jet, inj, 0, 0, 0);
    Real r_max = find_r_max(eqn, r_min, t_max);
    return {r_min, r_max};
}
#endif