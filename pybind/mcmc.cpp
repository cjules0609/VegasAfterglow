//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "mcmc.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

#include "pybind.h"

double MultiBandData::estimate_chi2() const {
    double chi_square = 0;
    for (size_t i = 0; i < times.size(); ++i) {
        double error = errors(i);
        if (error == 0) continue;
        double diff = fluxes(i) - model_fluxes(i);
        chi_square += (diff * diff) / (error * error);
    }
    return chi_square;
}

void MultiBandData::add_light_curve(double nu, PyArray const& t, PyArray const& Fv_obs, PyArray const& Fv_err) {
    assert(t.size() == Fv_obs.size() && t.size() == Fv_err.size() && "light curve array inconsistent length!");

    for (size_t i = 0; i < t.size(); i++) {
        tuple_data.push_back(std::make_tuple(t(i) * unit::sec, nu * unit::Hz, Fv_obs(i) * unit::flux_den_cgs,
                                             Fv_err(i) * unit::flux_den_cgs));
    }
}

void MultiBandData::add_spectrum(double t, PyArray const& nu, PyArray const& Fv_obs, PyArray const& Fv_err) {
    assert(nu.size() == Fv_obs.size() && nu.size() == Fv_err.size() && "spectrum array inconsistent length!");

    for (size_t i = 0; i < nu.size(); i++) {
        tuple_data.push_back(std::make_tuple(t * unit::sec, nu(i) * unit::Hz, Fv_obs(i) * unit::flux_den_cgs,
                                             Fv_err(i) * unit::flux_den_cgs));
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

    for (size_t i = 0; i < len; ++i) {
        times(i) = std::get<0>(tuple_data[i]);
        frequencies(i) = std::get<1>(tuple_data[i]);
        fluxes(i) = std::get<2>(tuple_data[i]);
        errors(i) = std::get<3>(tuple_data[i]);
        model_fluxes(i) = 0;  // Placeholder for model fluxes
    }
}

MultiBandModel::MultiBandModel(MultiBandData const& data) : obs_data(data) {
    obs_data.fill_data_arrays();

    if (obs_data.times.size() == 0) {
        std::cerr << "Error: No observation time data provided!" << std::endl;
    }
}

void MultiBandModel::configure(ConfigParams const& param) { this->config = param; }

double MultiBandModel::estimate_chi2(Params const& param) {
    build_system(param, obs_data.times, obs_data.frequencies, obs_data.model_fluxes);

    return obs_data.estimate_chi2();
}

auto MultiBandModel::light_curves(Params const& param, PyArray const& t, PyArray const& nu) -> PyGrid {
    Array t_bins = t * unit::sec;

    MeshGrid F_nu = MeshGrid::from_shape({nu.size(), t.size()});

    for (size_t i = 0; i < nu.size(); ++i) {
        Array nus = Array({t.size()}, nu(i) * unit::Hz);
        auto view = xt::view(F_nu, i, xt::all());
        build_system(param, t_bins, nus, view);
    }

    // we bind this function for GIL free. As the return will create a pyobject, we need to get the GIL.
    pybind11::gil_scoped_acquire acquire;
    return F_nu / unit::flux_den_cgs;
}

auto MultiBandModel::spectra(Params const& param, PyArray const& nu, PyArray const& t) -> PyGrid {
    Array t_bins = t * unit::sec;

    MeshGrid F_nu = MeshGrid::from_shape({nu.size(), t.size()});

    for (size_t i = 0; i < nu.size(); ++i) {
        Array nus = Array({t.size()}, nu(i) * unit::Hz);
        auto view = xt::view(F_nu, i, xt::all());
        build_system(param, t_bins, nus, view);
    }

    // we bind this function for GIL free. As the return will create a pyobject, we need to get the GIL.
    pybind11::gil_scoped_acquire acquire;
    return xt::transpose(F_nu / unit::flux_den_cgs);
}