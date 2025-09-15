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

using ArrayDict = std::unordered_map<std::string, xt::xarray<Real>>;

struct FluxData {
    Array t;
    Array nu;
    Array Fv_obs;
    Array Fv_err;
    Array Fv_model;
    Array weights;

    double estimate_chi2() const;
};

struct MultiBandData {
    double estimate_chi2() const;

    void add_light_curve(double nu, PyArray const& t, PyArray const& Fv_obs, PyArray const& Fv_err,
                         std::optional<PyArray> weights = std::nullopt);

    void add_light_curve(double nu_min, double nu_max, size_t num_points, PyArray const& t, PyArray const& Fv_obs,
                         PyArray const& Fv_err, std::optional<PyArray> weights = std::nullopt);

    void add_spectrum(double t, PyArray const& nu, PyArray const& Fv_obs, PyArray const& Fv_err,
                      std::optional<PyArray> weights = std::nullopt);

    std::vector<size_t> logscale_screen(PyArray const& data, size_t num_order);

    void fill_data_arrays();

    size_t data_points_num() const;

    Array times;
    Array frequencies;
    Array fluxes;
    Array errors;
    Array weights;
    Array model_fluxes;
    double t_min{con::inf};
    double t_max{0};

    std::vector<FluxData> flux_data;

  private:
    std::vector<std::tuple<double, double, double, double, double>> tuple_data;
};

struct Params {
    double theta_v{0};

    double n_ism{0};
    double n0{con::inf};
    double A_star{0};

    double E_iso{1e52};
    double Gamma0{300};
    double theta_c{0.1};
    double k_e{2};
    double k_g{2};
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
    double rtol{1e-6};
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

    PyGrid specific_flux(Params const& param, PyArray const& t, PyArray const& nu);

    PyArray flux(Params const& param, PyArray const& t, double nu_min, double nu_max, size_t num_points);

  private:
    template <typename ICPhotonGrid>
    void generate_photons(Params const& param, double t_min, double t_max, Observer& obs, SynPhotonGrid& f_photons,
                          SynPhotonGrid& r_photons, ICPhotonGrid& f_IC_photons, ICPhotonGrid& r_IC_photons);

    Ejecta select_jet(Params const& param);

    Medium select_medium(Params const& param);

    MultiBandData obs_data;

    ConfigParams config;
};

template <typename ICPhotonGrid>
void MultiBandModel::generate_photons(Params const& param, double t_min, double t_max, Observer& obs,
                                      SynPhotonGrid& fwd_photons, SynPhotonGrid& rvs_photons,
                                      ICPhotonGrid& fwd_IC_photons, ICPhotonGrid& rvs_IC_photons) {
    Real theta_v = param.theta_v;
    Real theta_w = param.theta_w;
    RadParams rad;
    rad.p = param.p;
    rad.eps_e = param.eps_e;
    rad.eps_B = param.eps_B;
    rad.xi_e = param.xi_e;

    Real lumi_dist = config.lumi_dist * unit::cm;
    Real z = config.z;

    Medium medium = select_medium(param);
    Ejecta jet = select_jet(param);

    Real t_resol = config.t_resol;
    Real theta_resol = config.theta_resol;
    Real phi_resol = config.phi_resol;

    Array t_eval = xt::linspace<Real>(t_min, t_max, 5);

    auto coord = auto_grid(jet, t_eval, theta_w, theta_v, z, phi_resol, theta_resol, t_resol);

    if (config.rvs_shock == false) {
        auto shock = generate_fwd_shock(coord, medium, jet, rad, config.rtol);

        obs.observe(coord, shock, lumi_dist, z);

        auto electrons = generate_syn_electrons(shock);
        fwd_photons = generate_syn_photons(shock, electrons);

        if (config.IC_cooling) {
            if (config.KN) {
                KN_cooling(electrons, fwd_photons, shock);
            } else {
                Thomson_cooling(electrons, fwd_photons, shock);
            }
        }

        if (config.fwd_SSC) {
            fwd_IC_photons = generate_IC_photons(electrons, fwd_photons, config.KN);
        }

    } else {
        RadParams rad_rvs;

        rad_rvs.p = param.p_r;
        rad_rvs.eps_e = param.eps_e_r;
        rad_rvs.eps_B = param.eps_B_r;
        rad_rvs.xi_e = param.xi_e_r;

        auto [f_shock, r_shock] = generate_shock_pair(coord, medium, jet, rad, rad_rvs, config.rtol);

        obs.observe(coord, f_shock, lumi_dist, z);

        auto f_electrons = generate_syn_electrons(f_shock);
        fwd_photons = generate_syn_photons(f_shock, f_electrons);

        auto r_electrons = generate_syn_electrons(r_shock);
        rvs_photons = generate_syn_photons(r_shock, r_electrons);

        if (config.IC_cooling) {
            if (config.KN) {
                KN_cooling(f_electrons, fwd_photons, f_shock);
                KN_cooling(r_electrons, rvs_photons, r_shock);
            } else {
                Thomson_cooling(f_electrons, fwd_photons, f_shock);
                Thomson_cooling(r_electrons, rvs_photons, r_shock);
            }
        }

        if (config.fwd_SSC) {
            fwd_IC_photons = generate_IC_photons(f_electrons, fwd_photons, config.KN);
        }

        if (config.rvs_SSC) {
            rvs_IC_photons = generate_IC_photons(r_electrons, rvs_photons, config.KN);
        }
    }
}
