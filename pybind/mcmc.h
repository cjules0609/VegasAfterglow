//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <iostream>
#include <vector>

#include "afterglow.h"
#include "macros.h"
#include "mesh.h"
#include "pybind.h"
#include "utilities.h"
struct MultiBandData {
    double estimate_chi2() const;

    void add_light_curve(double nu, PyArray const& t, PyArray const& Fv_obs, PyArray const& Fv_err);

    void add_spectrum(double t, PyArray const& nu, PyArray const& Fv_obs, PyArray const& Fv_err);

    void fill_data_arrays();

    Array times;
    Array frequencies;
    Array fluxes;
    Array errors;
    Array model_fluxes;

   private:
    std::vector<std::tuple<double, double, double, double>> tuple_data;
};

struct Params {
    double theta_v{0};

    double n_ism{1};
    double n0{con::inf};
    double A_star{0};

    double E_iso{1e52};
    double Gamma0{300};
    double theta_c{0.1};
    double k{2};
    double duration{1 * unit::sec};

    double E_iso_w{1e52};
    double Gamma0_w{300};
    double theta_w{con::pi / 2};

    double L0{0};
    double t0{1};
    double q{2};

    double p{2.3};
    double eps_e{0.1};
    double eps_B{0.01};
    double xi_e{1};

    double p_r{2.3};
    double eps_e_r{0.1};
    double eps_B_r{0.01};
    double xi_e_r{1};
};

struct ConfigParams {
    double lumi_dist{1e26};
    double z{0};
    std::string medium{"ism"};
    std::string jet{"tophat"};
    Real phi_resol{0.3};
    Real theta_resol{1};
    Real t_resol{10};
    double rtol{1e-5};
    bool rvs_shock{false};
    bool fwd_SSC{false};
    bool rvs_SSC{false};
    bool IC_cooling{false};
    bool KN{false};
    bool magnetar{false};
};

struct MultiBandModel {
    MultiBandModel() = delete;
    MultiBandModel(MultiBandData const& data);

    void configure(ConfigParams const& param);
    double estimate_chi2(Params const& param);
    PyGrid light_curves(Params const& param, PyArray const& t, PyArray const& nu);
    PyGrid spectra(Params const& param, PyArray const& nu, PyArray const& t);

   private:
    template <typename View>
    void build_system(Params const& param, Array const& t_eval, Array const& nu_eval, View& F_nu);
    MultiBandData obs_data;
    ConfigParams config;
};

template <typename View>
void MultiBandModel::build_system(Params const& param, Array const& t_eval, Array const& nu_eval, View& F_nu) {
    Real eps_iso = param.E_iso * unit::erg / (4 * con::pi);
    Real Gamma0 = param.Gamma0;
    Real theta_c = param.theta_c;
    Real theta_v = param.theta_v;
    Real theta_w = param.theta_w;
    Real eps_iso_w = param.E_iso_w * unit::erg / (4 * con::pi);
    Real Gamma0_w = param.Gamma0_w;
    RadParams rad;

    rad.p = param.p;
    rad.eps_e = param.eps_e;
    rad.eps_B = param.eps_B;
    rad.xi_e = param.xi_e;

    Real lumi_dist = config.lumi_dist * unit::cm;
    Real z = config.z;

    // create model
    Medium medium;
    if (config.medium == "ism") {
        Real n_ism = param.n_ism / unit::cm3;
        medium.rho = evn::ISM(n_ism);
    } else if (config.medium == "wind") {
        medium.rho = evn::wind(param.A_star, param.n_ism, param.n0);
    } else {
        std::cerr << "Error: Unknown medium type" << std::endl;
    }

    Ejecta jet;
    jet.T0 = param.duration * unit::sec;
    if (config.jet == "tophat") {
        jet.eps_k = math::tophat(theta_c, eps_iso);
        jet.Gamma0 = math::tophat(theta_c, Gamma0);
    } else if (config.jet == "gaussian") {
        jet.eps_k = math::gaussian(theta_c, eps_iso);
        jet.Gamma0 = math::gaussian_plus_one(theta_c, Gamma0 - 1);
    } else if (config.jet == "powerlaw") {
        jet.eps_k = math::powerlaw(theta_c, eps_iso, param.k);
        jet.Gamma0 = math::powerlaw_plus_one(theta_c, Gamma0 - 1, param.k);
    } else if (config.jet == "twocomponent") {
        jet.eps_k = math::two_component(theta_c, theta_w, eps_iso, eps_iso_w);
        jet.Gamma0 = math::two_component_plus_one(theta_c, theta_w, Gamma0 - 1, Gamma0_w - 1);
    } else if (config.jet == "steppowerlaw") {
        jet.eps_k = math::step_powerlaw(theta_c, eps_iso, eps_iso_w, param.k);
        jet.Gamma0 = math::step_powerlaw_plus_one(theta_c, Gamma0 - 1, Gamma0_w - 1, param.k);
    } else {
        std::cerr << "Error: Unknown jet type" << std::endl;
    }

    if (config.magnetar == true) {
        jet.deps_dt =
            math::magnetar_injection(param.t0 * unit::sec, param.q, param.L0 * unit::erg / unit::sec, theta_c);
    }

    Real t_resol = config.t_resol;
    Real theta_resol = config.theta_resol;
    Real phi_resol = config.phi_resol;

    auto coord = auto_grid(jet, t_eval, theta_w, theta_v, z, phi_resol, theta_resol, t_resol);

    if (config.rvs_shock == false) {
        auto shock = generate_fwd_shock(coord, medium, jet, rad, config.rtol);

        Observer obs;
        // obs.observe_at(t_eval, coord, shock, lumi_dist, z);
        obs.observe(coord, shock, lumi_dist, z);

        auto electrons = generate_syn_electrons(shock);
        auto photons = generate_syn_photons(shock, electrons);

        if (config.IC_cooling) {
            if (config.KN) {
                KN_cooling(electrons, photons, shock);
            } else {
                Thomson_cooling(electrons, photons, shock);
            }
        }

        F_nu = obs.specific_flux_series(t_eval, nu_eval, photons);

        if (config.fwd_SSC) {
            auto IC_photon = generate_IC_photons(electrons, photons, config.KN);
            auto IC_F_nu = obs.specific_flux_series(t_eval, nu_eval, IC_photon);
            F_nu += IC_F_nu;
        }

    } else {
        RadParams rad_rvs;

        rad_rvs.p = param.p_r;
        rad_rvs.eps_e = param.eps_e_r;
        rad_rvs.eps_B = param.eps_B_r;
        rad_rvs.xi_e = param.xi_e_r;

        auto [f_shock, r_shock] = generate_shock_pair(coord, medium, jet, rad, rad_rvs, config.rtol);

        Observer obs;
        // obs.observe_at(t_eval, coord, shock, lumi_dist, z);
        obs.observe(coord, f_shock, lumi_dist, z);

        auto f_electrons = generate_syn_electrons(f_shock);
        auto f_photons = generate_syn_photons(f_shock, f_electrons);

        auto r_electrons = generate_syn_electrons(r_shock);
        auto r_photons = generate_syn_photons(r_shock, r_electrons);

        if (config.IC_cooling) {
            if (config.KN) {
                KN_cooling(f_electrons, f_photons, f_shock);
                KN_cooling(r_electrons, r_photons, r_shock);
            } else {
                Thomson_cooling(f_electrons, f_photons, f_shock);
                Thomson_cooling(r_electrons, r_photons, r_shock);
            }
        }

        F_nu = obs.specific_flux_series(t_eval, nu_eval, f_photons);
        auto F_nu_rvs = obs.specific_flux_series(t_eval, nu_eval, r_photons);
        F_nu += F_nu_rvs;

        if (config.fwd_SSC) {
            auto IC_photon_f = generate_IC_photons(f_electrons, f_photons, config.KN);
            auto IC_F_nu_fwd = obs.specific_flux_series(t_eval, nu_eval, IC_photon_f);
            F_nu += IC_F_nu_fwd;
        }

        if (config.rvs_SSC) {
            auto IC_photon_r = generate_IC_photons(r_electrons, r_photons, config.KN);
            auto IC_F_nu_rvs = obs.specific_flux_series(t_eval, nu_eval, IC_photon_r);
            F_nu += IC_F_nu_rvs;
        }
    }
}