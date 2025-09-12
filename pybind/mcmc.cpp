//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "mcmc.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>

#include "pybind.h"

double FluxData::estimate_chi2() const {
    double chi_square = 0;
    for (size_t i = 0; i < t.size(); ++i) {
        double error = Fv_err(i);
        if (error == 0)
            continue;
        double diff = Fv_obs(i) - Fv_model(i);
        chi_square += weights(i) * (diff * diff) / (error * error);
    }
    return chi_square;
}

double MultiBandData::estimate_chi2() const {
    double chi_square = 0;
    for (size_t i = 0; i < times.size(); ++i) {
        double error = errors(i);
        if (error == 0)
            continue;
        double diff = fluxes(i) - model_fluxes(i);
        chi_square += weights(i) * (diff * diff) / (error * error);
    }
    for (auto& d : flux_data) {
        chi_square += d.estimate_chi2();
    }
    return chi_square;
}

Ejecta MultiBandModel::select_jet(Params const& param) {
    Real eps_iso = param.E_iso * unit::erg / (4 * con::pi);
    Real Gamma0 = param.Gamma0;
    Real theta_c = param.theta_c;
    Real theta_w = param.theta_w;
    Real eps_iso_w = param.E_iso_w * unit::erg / (4 * con::pi);
    Real Gamma0_w = param.Gamma0_w;
    Ejecta jet;
    jet.T0 = param.duration * unit::sec;
    if (config.jet == "tophat") {
        jet.eps_k = math::tophat(theta_c, eps_iso);
        jet.Gamma0 = math::tophat_plus_one(theta_c, Gamma0 - 1);
    } else if (config.jet == "gaussian") {
        jet.eps_k = math::gaussian(theta_c, eps_iso);
        jet.Gamma0 = math::gaussian_plus_one(theta_c, Gamma0 - 1);
    } else if (config.jet == "powerlaw") {
        jet.eps_k = math::powerlaw(theta_c, eps_iso, param.k_e);
        jet.Gamma0 = math::powerlaw_plus_one(theta_c, Gamma0 - 1, param.k_g);
    } else if (config.jet == "twocomponent") {
        jet.eps_k = math::two_component(theta_c, theta_w, eps_iso, eps_iso_w);
        jet.Gamma0 = math::two_component_plus_one(theta_c, theta_w, Gamma0 - 1, Gamma0_w - 1);
    } else if (config.jet == "steppowerlaw") {
        jet.eps_k = math::step_powerlaw(theta_c, eps_iso, eps_iso_w, param.k_e);
        jet.Gamma0 = math::step_powerlaw_plus_one(theta_c, Gamma0 - 1, Gamma0_w - 1, param.k_g);
    } else {
        assert(false && "Unknown jet type");
    }

    if (config.magnetar == true) {
        jet.deps_dt =
            math::magnetar_injection(param.t0 * unit::sec, param.q, param.L0 * unit::erg / unit::sec, theta_c);
    }
    return jet;
}

Medium MultiBandModel::select_medium(Params const& param) {
    Medium medium;
    if (config.medium == "ism") {
        medium.rho = evn::ISM(param.n_ism / unit::cm3);
    } else if (config.medium == "wind") {
        medium.rho = evn::wind(param.A_star, param.n_ism / unit::cm3, param.n0 / unit::cm3);
    } else {
        assert(false && "Unknown medium type");
    }
    return medium;
}

void MultiBandData::add_light_curve(double nu, PyArray const& t, PyArray const& Fv_obs, PyArray const& Fv_err,
                                    std::optional<PyArray> weights) {
    assert(t.size() == Fv_obs.size() && t.size() == Fv_err.size() && "light curve array inconsistent length!");

    Array w = xt::ones<Real>({t.size()});

    if (weights) {
        w = *weights;
        assert(t.size() == w.size() && "weights array inconsistent length!");
    }

    for (size_t i = 0; i < t.size(); i++) {
        tuple_data.push_back(std::make_tuple(t(i) * unit::sec, nu * unit::Hz, Fv_obs(i) * unit::flux_den_cgs,
                                             Fv_err(i) * unit::flux_den_cgs, w(i)));
    }
}

void MultiBandData::add_light_curve(double nu_min, double nu_max, size_t num_points, PyArray const& t,
                                    PyArray const& Fv_obs, PyArray const& Fv_err, std::optional<PyArray> weights) {
    assert(t.size() == Fv_obs.size() && t.size() == Fv_err.size() && "light curve array inconsistent length!");
    assert(is_ascending(t) && "Time array must be in ascending order!");
    assert(nu_min < nu_max && "nu_min must be less than nu_max!");

    Array w = xt::ones<Real>({t.size()});

    if (weights) {
        w = *weights;
        assert(t.size() == w.size() && "weights array inconsistent length!");

        size_t len = w.size();
        Real weight_sum = 0;
        for (size_t i = 0; i < len; ++i) {
            weight_sum += w(i);
        }
        w /= (weight_sum / len);
    }

    Array nu = xt::logspace(std::log10(nu_min * unit::Hz), std::log10(nu_max * unit::Hz), num_points);

    flux_data.emplace_back(
        FluxData{t * unit::sec, nu, Fv_obs * unit::flux_cgs, Fv_err * unit::flux_cgs, xt::zeros<Real>({t.size()}), w});
}

void MultiBandData::add_spectrum(double t, PyArray const& nu, PyArray const& Fv_obs, PyArray const& Fv_err,
                                 std::optional<PyArray> weights) {
    assert(nu.size() == Fv_obs.size() && nu.size() == Fv_err.size() && "spectrum array inconsistent length!");

    Array w = xt::ones<Real>({nu.size()});

    if (weights) {
        w = *weights;
        assert(nu.size() == w.size() && "weights array inconsistent length!");
    }

    for (size_t i = 0; i < nu.size(); i++) {
        tuple_data.push_back(std::make_tuple(t * unit::sec, nu(i) * unit::Hz, Fv_obs(i) * unit::flux_den_cgs,
                                             Fv_err(i) * unit::flux_den_cgs, w(i)));
    }
}

void MultiBandData::fill_data_arrays() {
    const size_t len = tuple_data.size();
    std::sort(tuple_data.begin(), tuple_data.end(),
              [](auto const& a, auto const& b) { return std::get<0>(a) < std::get<0>(b); });
    times = Array::from_shape({len});
    frequencies = Array::from_shape({len});
    fluxes = Array::from_shape({len});
    errors = Array::from_shape({len});
    model_fluxes = Array::from_shape({len});
    weights = Array::from_shape({len});

    Real weight_sum = 0;
    for (size_t i = 0; i < len; ++i) {
        times(i) = std::get<0>(tuple_data[i]);
        frequencies(i) = std::get<1>(tuple_data[i]);
        fluxes(i) = std::get<2>(tuple_data[i]);
        errors(i) = std::get<3>(tuple_data[i]);
        weights(i) = std::get<4>(tuple_data[i]);
        model_fluxes(i) = 0; // Placeholder for model fluxes
        weight_sum += weights(i);
    }
    weights /= (weight_sum / len);

    if (len > 0) {
        this->t_min = times.front();
        this->t_max = times.back();
    }

    for (auto& d : flux_data) {
        if (d.t.front() < t_min)
            t_min = d.t.front();
        if (d.t.back() > t_max)
            t_max = d.t.back();
    }
}

MultiBandModel::MultiBandModel(MultiBandData const& data) : obs_data(data) {
    obs_data.fill_data_arrays();

    assert((obs_data.times.size() > 0 || obs_data.fluxes.size() > 0) && "No observation time data provided!");
}

void MultiBandModel::configure(ConfigParams const& param) {
    this->config = param;
}

double MultiBandModel::estimate_chi2(Params const& param) {
    Observer obs;
    SynPhotonGrid f_photons;
    SynPhotonGrid r_photons;
    ICPhotonGrid<SynElectrons, SynPhotons> f_IC_photons;
    ICPhotonGrid<SynElectrons, SynPhotons> r_IC_photons;

    generate_photons(param, obs_data.t_min, obs_data.t_max, obs, f_photons, r_photons, f_IC_photons, r_IC_photons);

    obs_data.model_fluxes = obs.specific_flux_series(obs_data.times, obs_data.frequencies, f_photons);
    for (auto& d : obs_data.flux_data) {
        d.Fv_model = obs.flux(d.t, d.nu, f_photons);
    }

    if (r_photons.size() > 0) {
        obs_data.model_fluxes += obs.specific_flux_series(obs_data.times, obs_data.frequencies, r_photons);
        for (auto& d : obs_data.flux_data) {
            d.Fv_model += obs.flux(d.t, d.nu, r_photons);
        }
    }

    if (f_IC_photons.size() > 0) {
        obs_data.model_fluxes += obs.specific_flux_series(obs_data.times, obs_data.frequencies, f_IC_photons);
        for (auto& d : obs_data.flux_data) {
            d.Fv_model += obs.flux(d.t, d.nu, f_IC_photons);
        }
    }

    if (r_IC_photons.size() > 0) {
        obs_data.model_fluxes += obs.specific_flux_series(obs_data.times, obs_data.frequencies, r_IC_photons);
        for (auto& d : obs_data.flux_data) {
            d.Fv_model += obs.flux(d.t, d.nu, r_IC_photons);
        }
    }

    return obs_data.estimate_chi2();
}

auto MultiBandModel::specific_flux(Params const& param, PyArray const& t, PyArray const& nu) -> PyGrid {
    Array t_bins = t * unit::sec;
    Array nu_bins = nu * unit::Hz;
    MeshGrid F_nu = MeshGrid::from_shape({nu.size(), t.size()});

    Observer obs;
    SynPhotonGrid f_photons;
    SynPhotonGrid r_photons;
    ICPhotonGrid<SynElectrons, SynPhotons> f_IC_photons;
    ICPhotonGrid<SynElectrons, SynPhotons> r_IC_photons;

    generate_photons(param, t_bins.front(), t_bins.back(), obs, f_photons, r_photons, f_IC_photons, r_IC_photons);

    F_nu = obs.specific_flux(t_bins, nu_bins, f_photons);

    if (r_photons.size() > 0) {
        F_nu += obs.specific_flux(t_bins, nu_bins, r_photons);
    }

    if (f_IC_photons.size() > 0) {
        F_nu += obs.specific_flux(t_bins, nu_bins, f_IC_photons);
    }

    if (r_IC_photons.size() > 0) {
        F_nu += obs.specific_flux(t_bins, nu_bins, r_IC_photons);
    }

    // we bind this function for GIL free. As the return will create a pyobject, we need to get the GIL.
    pybind11::gil_scoped_acquire acquire;
    return F_nu / unit::flux_den_cgs;
}

auto MultiBandModel::flux(Params const& param, PyArray const& t, double nu_min, double nu_max, size_t num_points)
    -> PyArray {
    Array t_bins = t * unit::sec;
    Array nu_bins = xt::logspace(std::log10(nu_min * unit::Hz), std::log10(nu_max * unit::Hz), num_points);
    Array F_nu;

    Observer obs;
    SynPhotonGrid f_photons;
    SynPhotonGrid r_photons;
    ICPhotonGrid<SynElectrons, SynPhotons> f_IC_photons;
    ICPhotonGrid<SynElectrons, SynPhotons> r_IC_photons;

    generate_photons(param, t_bins.front(), t_bins.back(), obs, f_photons, r_photons, f_IC_photons, r_IC_photons);

    F_nu = obs.flux(t_bins, nu_bins, f_photons);

    if (r_photons.size() > 0) {
        F_nu += obs.flux(t_bins, nu_bins, r_photons);
    }

    if (f_IC_photons.size() > 0) {
        F_nu += obs.flux(t_bins, nu_bins, f_IC_photons);
    }

    if (r_IC_photons.size() > 0) {
        F_nu += obs.flux(t_bins, nu_bins, r_IC_photons);
    }

    // we bind this function for GIL free. As the return will create a pyobject, we need to get the GIL.
    pybind11::gil_scoped_acquire acquire;
    return F_nu / unit::flux_cgs;
}
