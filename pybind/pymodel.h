//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <iostream>
#include <optional>
#include <vector>

#include "afterglow.h"
#include "macros.h"
#include "pybind.h"

Ejecta PyTophatJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading = false, Real T0 = 1 * unit::sec);

Ejecta PyGaussianJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading = false, Real T0 = 1 * unit::sec);

Ejecta PyPowerLawJet(Real theta_c, Real E_iso, Real Gamma0, Real k, bool spreading = false, Real T0 = 1 * unit::sec);

Medium PyISM(Real n_ism);

Medium PyWind(Real A_star);

class PyObserver {
   public:
    PyObserver(Real lumi_dist, Real z, Real theta_obs, Real phi_obs = 0)
        : lumi_dist(lumi_dist * unit::cm), z(z), theta_obs(theta_obs), phi_obs(phi_obs) {}

    Real lumi_dist{1e28};
    Real z{0};
    Real theta_obs{0};
    Real phi_obs{0};
};

class PyRadiation {
   public:
    PyRadiation(Real eps_e, Real eps_B, Real p, Real xi_e, bool SSC = false, bool Klein_Nishina = true)
        : eps_e(eps_e), eps_B(eps_B), p(p), xi_e(xi_e), SSC(SSC), Klein_Nishina(Klein_Nishina) {}

    Real eps_e{1e-1};
    Real eps_B{1e-2};
    Real p{2.3};
    Real xi_e{1};
    bool SSC{false};
    bool Klein_Nishina{true};
};

class PyModel {
    using FluxDict = std::unordered_map<std::string, MeshGrid>;

   public:
    PyModel(Ejecta jet, Medium medium, PyObserver observer, PyRadiation fwd_rad,
            std::optional<PyRadiation> rvs_rad = std::nullopt)
        : jet(jet), medium(medium), obs_setup(observer), fwd_rad(fwd_rad), rvs_rad_opt(rvs_rad) {}

    FluxDict specific_flux(PyArray const& t, PyArray const& nu);

    FluxDict spectra(PyArray const& nu, PyArray const& t);

   private:
    FluxDict specific_flux_(Array const& t, Array const& nu);

    void specific_flux_for(Shock const& shock, Coord const& coord, Array const& t, Array const& nu, Observer& obs,
                           PyRadiation rad, FluxDict& flux_dict, std::string suffix);

    Ejecta jet;
    Medium medium;
    PyObserver obs_setup;
    PyRadiation fwd_rad;
    std::optional<PyRadiation> rvs_rad_opt;
    Real theta_w{con::pi / 2};
    size_t phi_num{64};
    size_t theta_num{64};
    size_t t_num{100};
    Real rtol{1e-5};
    bool axisymmetric{true};
};