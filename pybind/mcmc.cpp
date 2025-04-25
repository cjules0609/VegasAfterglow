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

template <typename Array>
void sort_synchronized(Array& a, Array& b, Array& c) {
    std::size_t N = a.size();

    std::vector<std::size_t> idx(N);
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(), [&](std::size_t i, std::size_t j) { return a(i) < a(j); });

    Array a_sorted = xt::empty<double>({N});
    Array b_sorted = xt::empty<double>({N});
    Array c_sorted = xt::empty<double>({N});

    for (std::size_t i = 0; i < N; ++i) {
        a_sorted(i) = a(idx[i]);
        b_sorted(i) = b(idx[i]);
        c_sorted(i) = c(idx[i]);
    }

    a = std::move(a_sorted);
    b = std::move(b_sorted);
    c = std::move(c_sorted);
}

double LightCurveData::calcChiSquare() const {
    double chi_square = 0;
    for (size_t i = 0; i < t.size(); ++i) {
        if (Fv_err[i] == 0) continue;
        double diff = Fv_obs[i] - Fv_model[i];
        chi_square += (diff * diff) / (Fv_err[i] * Fv_err[i]);
    }
    return chi_square;
}

double SpectrumData::calcChiSquare() const {
    double chi_square = 0;
    for (size_t i = 0; i < nu.size(); ++i) {
        if (Fv_err[i] == 0) continue;
        double diff = Fv_obs[i] - Fv_model[i];
        chi_square += (diff * diff) / (Fv_err[i] * Fv_err[i]);
    }
    return chi_square;
}

double MultiBandData::calcChiSquare() const {
    double chi_square = 0;
    for (auto const& data : light_curve) {
        chi_square += data.calcChiSquare();
    }
    for (auto const& data : spectrum) {
        chi_square += data.calcChiSquare();
    }
    return chi_square;
}

void MultiBandData::addObsLightCurve(double nu, List const& t, List const& Fv_obs, List const& Fv_err) {
    assert(t.size() == Fv_obs.size() && t.size() == Fv_err.size() && "light curve array inconsistent length!");
    LightCurveData data;

    data.nu = nu * con::Hz;
    data.t = xt::eval(xt::adapt(t));
    data.Fv_obs = xt::eval(xt::adapt(Fv_obs));
    data.Fv_err = xt::eval(xt::adapt(Fv_err));
    data.Fv_model = xt::zeros_like(data.Fv_obs);

    sort_synchronized(data.t, data.Fv_obs, data.Fv_err);

    for (auto& t : data.t) {
        t *= con::sec;
    }
    for (auto& Fv_obs : data.Fv_obs) {
        Fv_obs *= con::erg / con::sec / con::cm2 / con::Hz;
    }
    for (auto& Fv_err : data.Fv_err) {
        Fv_err *= con::erg / con::sec / con::cm2 / con::Hz;
    }

    light_curve.push_back(std::move(data));
}

void MultiBandData::addObsSpectrum(double t, List const& nu, List const& Fv_obs, List const& Fv_err) {
    assert(nu.size() == Fv_obs.size() && nu.size() == Fv_err.size() && "spectrum array inconsistent length!");
    SpectrumData data;

    data.t = t * con::sec;
    data.nu = xt::eval(xt::adapt(nu));
    data.Fv_obs = xt::eval(xt::adapt(Fv_obs));
    data.Fv_err = xt::eval(xt::adapt(Fv_err));

    sort_synchronized(data.nu, data.Fv_obs, data.Fv_err);

    for (auto& nu : data.nu) {
        nu *= con::Hz;
    }
    for (auto& Fv_obs : data.Fv_obs) {
        Fv_obs *= con::erg / con::sec / con::cm2 / con::Hz;
    }
    for (auto& Fv_err : data.Fv_err) {
        Fv_err *= con::erg / con::sec / con::cm2 / con::Hz;
    }
    spectrum.push_back(std::move(data));
}

MultiBandModel::MultiBandModel(MultiBandData const& obs_data) : obs_data(obs_data) {
    std::vector<Real> evaluate_t;
    evaluate_t.reserve(100);

    for (auto& data : obs_data.light_curve) {
        for (auto t : data.t) {
            evaluate_t.push_back(t);
        }
    }

    for (auto& data : obs_data.spectrum) {
        evaluate_t.push_back(data.t);
    }
    std::sort(evaluate_t.begin(), evaluate_t.end());
    this->t_eval = xt::eval(xt::adapt(evaluate_t));

    if (t_eval.size() == 0) {
        std::cerr << "Error: No observation time data provided!" << std::endl;
    }
}

void MultiBandModel::configure(ConfigParams const& param) { this->config = param; }

/*
std::vector<double> MultiBandModel::chiSquareBatch(std::vector<Params> const& param_batch) {
    const size_t N = param_batch.size();
    std::vector<double> results(N);

#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        results[i] = chiSquare(param_batch[i]);
    }

    return results;
}
*/

void MultiBandModel::buildSystem(Params const& param, Array const& t_eval, Observer& obs, SynElectronGrid& electrons,
                                 SynPhotonGrid& photons) {
    Real E_iso = param.E_iso * con::erg;
    Real Gamma0 = param.Gamma0;
    Real theta_c = param.theta_c;
    Real theta_v = param.theta_v;
    Real theta_w = param.theta_w;
    Real p = param.p;
    Real eps_e = param.eps_e;
    Real eps_B = param.eps_B;
    Real xi = param.xi;

    Real lumi_dist = config.lumi_dist * con::cm;
    Real z = config.z;

    // create model
    Medium medium;
    if (config.medium == "ism") {
        Real n_ism = param.n_ism / con::cm3;
        medium.rho = evn::ISM(n_ism);
    } else if (config.medium == "wind") {
        medium.rho = evn::wind(param.A_star);
    } else {
        std::cerr << "Error: Unknown medium type" << std::endl;
    }

    Ejecta jet;
    if (config.jet == "tophat") {
        jet.dE0dOmega = math::tophat(theta_c, E_iso / (4 * con::pi));
        jet.Gamma0 = math::tophat(theta_c, Gamma0);
    } else if (config.jet == "gaussian") {
        jet.dE0dOmega = math::gaussian(theta_c, E_iso / (4 * con::pi));
        jet.Gamma0 = math::gaussian(theta_c, Gamma0);
    } else if (config.jet == "powerlaw") {
        jet.dE0dOmega = math::powerLaw(theta_c, E_iso / (4 * con::pi), param.k_jet);
        jet.Gamma0 = math::powerLaw(theta_c, Gamma0, param.k_jet);
    } else {
        std::cerr << "Error: Unknown jet type" << std::endl;
    }

    size_t t_num = config.t_grid;
    size_t theta_num = config.theta_grid;
    size_t phi_num = config.phi_grid;

    auto coord = autoGrid(jet, t_eval, theta_w, theta_v, phi_num, theta_num, t_num);

    auto shock = genForwardShock(coord, medium, jet, eps_e, eps_B, config.rtol);

    obs.observeAt(t_eval, coord, shock, theta_v, lumi_dist, z);

    electrons = genSynElectrons(shock, p, xi);

    photons = genSynPhotons(shock, electrons);
}

double MultiBandModel::chiSquare(Params const& param) {
    Observer obs;
    SynElectronGrid electrons;
    SynPhotonGrid photons;

    buildSystem(param, t_eval, obs, electrons, photons);

    for (auto& data : obs_data.light_curve) {
        data.Fv_model = obs.specificFlux(data.t, data.nu, photons);
    }

    for (auto& data : obs_data.spectrum) {
        data.Fv_model = obs.spectrum(data.nu, data.t, photons);
    }

    return obs_data.calcChiSquare();
}

std::vector<std::vector<double>> F_nu2Grid(MeshGrid const& grid) {
    std::vector<std::vector<double>> out;
    out.reserve(grid.shape()[0]);
    double unit = con::erg / con::sec / con::cm2 / con::Hz;

    for (size_t i = 0; i < grid.shape()[0]; ++i) {
        std::vector<double> row;
        row.reserve(grid.shape()[1]);
        for (size_t j = 0; j < grid.shape()[1]; ++j) {
            row.push_back(grid(i, j) / unit);
        }
        out.push_back(std::move(row));
    }
    return out;
}

auto MultiBandModel::lightCurve(Params const& param, List const& t, List const& nu) -> Grid {
    Observer obs;
    SynElectronGrid electrons;
    SynPhotonGrid photons;

    Array t_bins = xt::adapt(t) * con::sec;
    buildSystem(param, t_bins, obs, electrons, photons);
    auto nu_bins = xt::adapt(nu) * con::Hz;
    auto F_nu = obs.specificFlux(t_bins, nu_bins, photons);

    return F_nu2Grid(F_nu);
}

auto MultiBandModel::spectrum(Params const& param, List const& nu, List const& t) -> Grid {
    Observer obs;
    SynElectronGrid electrons;
    SynPhotonGrid photons;

    Array t_bins = xt::adapt(t) * con::sec;
    buildSystem(param, t_bins, obs, electrons, photons);
    auto nu_bins = xt::adapt(nu) * con::Hz;
    auto F_nu = obs.spectrum(nu_bins, t_bins, photons);

    return F_nu2Grid(F_nu);
}