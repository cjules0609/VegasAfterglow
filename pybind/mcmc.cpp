//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "mcmc.h"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

template <typename T>
void sort_synchronized(T& a, T& b, T& c) {
    std::vector<size_t> idx(a.size());
    for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;

    std::sort(idx.begin(), idx.end(), [&](size_t i, size_t j) { return a[i] < a[j]; });

    std::vector<double> a_sorted, b_sorted, c_sorted;
    a_sorted.reserve(a.size());
    b_sorted.reserve(b.size());
    c_sorted.reserve(c.size());

    for (size_t i : idx) {
        a_sorted.push_back(a[i]);
        b_sorted.push_back(b[i]);
        c_sorted.push_back(c[i]);
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
void LightCurveData::resize(size_t size) {
    t.resize(boost::extents[size]);
    Fv_obs.resize(boost::extents[size]);
    Fv_err.resize(boost::extents[size]);
    Fv_model.resize(boost::extents[size]);
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

void SpectrumData::resize(size_t size) {
    nu.resize(boost::extents[size]);
    Fv_obs.resize(boost::extents[size]);
    Fv_err.resize(boost::extents[size]);
    Fv_model.resize(boost::extents[size]);
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

void MultiBandData::addObsLightCurve(double nu, List& t, List& Fv_obs, List& Fv_err) {
    assert(t.size() == Fv_obs.size() && t.size() == Fv_err.size() && "light curve array inconsistent length!");
    LightCurveData data;
    sort_synchronized(t, Fv_obs, Fv_err);
    size_t size = t.size();
    data.resize(size);
    data.nu = nu * con::Hz;

    for (size_t i = 0; i < size; ++i) {
        data.t[i] = t[i] * con::sec;
    }
    for (size_t i = 0; i < size; ++i) {
        data.Fv_obs[i] = Fv_obs[i] * con::erg / con::sec / con::cm2 / con::Hz;
    }
    for (size_t i = 0; i < size; ++i) {
        data.Fv_err[i] = Fv_err[i] * con::erg / con::sec / con::cm2 / con::Hz;
    }
    light_curve.push_back(std::move(data));
}

void MultiBandData::addObsSpectrum(double t, List& nu, List& Fv_obs, List& Fv_err) {
    assert(nu.size() == Fv_obs.size() && nu.size() == Fv_err.size() && "spectrum array inconsistent length!");
    SpectrumData data;
    sort_synchronized(nu, Fv_obs, Fv_err);
    size_t size = nu.size();
    data.resize(size);
    data.t = t * con::sec;

    for (size_t i = 0; i < size; ++i) {
        data.nu[i] = nu[i] * con::Hz;
    }
    for (size_t i = 0; i < size; ++i) {
        data.Fv_obs[i] = Fv_obs[i] * con::erg / con::sec / con::cm2 / con::Hz;
    }
    for (size_t i = 0; i < size; ++i) {
        data.Fv_err[i] = Fv_err[i] * con::erg / con::sec / con::cm2 / con::Hz;
    }
    spectrum.push_back(std::move(data));
}

MultiBandModel::MultiBandModel(MultiBandData const& obs_data)
    : obs_data(obs_data),
      coord(zeros(32), zeros(32), zeros(32)),
      shock(Shock(1, 32, 32, 0, 0)),
      electrons(createSynElectronGrid(1, 32, 32)),
      photons(createSynPhotonGrid(1, 32, 32)),
      obs(coord, shock, 0, 1, 0) {
    for (auto const& data : obs_data.light_curve) {
        for (auto t : data.t) {
            if (t_min == 0) {
                t_min = t / 2;
                t_max = t * 2;
            }
            if (t < t_min) t_min = t;
            if (t > t_max) t_max = t;
        }
    }

    if (t_min == 0 && t_max == 0) {
        std::cerr << "Error: No observation time data provided!" << std::endl;
    }
}

void MultiBandModel::configure(ConfigParams const& param) {
    this->config = param;
    this->coord.phi.resize(boost::extents[config.phi_grid]);
    this->coord.theta.resize(boost::extents[config.theta_grid]);
    this->coord.t.resize(boost::extents[config.t_grid]);
    this->shock.resize(1, config.theta_grid, config.t_grid);
    this->electrons.resize(boost::extents[1][config.theta_grid][config.t_grid]);
    this->photons.resize(boost::extents[1][config.theta_grid][config.t_grid]);
    this->obs.resize(config.phi_grid, config.theta_grid, config.t_grid);
    this->obs.lumi_dist = config.lumi_dist * con::cm;
    this->obs.z = config.z;
}
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
std::vector<double> MultiBandModel::chiSquareBatch(std::vector<Params> const& param_batch) {
    const size_t N = param_batch.size();
    std::vector<double> results(N);

#pragma omp parallel
    {
        MultiBandModel model_local = *this;

#pragma omp for schedule(dynamic)
        for (int i = 0; i < N; ++i) {
            results[i] = model_local.chiSquare(param_batch[i]);
        }
    }

    return results;
}

double MultiBandModel::chiSquare(Params const& param) {
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
        return -1;
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
        return -1;
    }

    size_t t_num = config.t_grid;
    size_t theta_num = config.theta_grid;
    size_t phi_num = config.phi_grid;

    Array t_bins = logspace(t_min, t_max, t_num);

    autoGrid(this->coord, jet, t_bins, theta_w, theta_v, phi_num, theta_num, t_num);

    genForwardShock(this->shock, this->coord, medium, jet, eps_e, eps_B, config.rtol);

    genSynElectrons(this->electrons, this->shock, p, xi);

    genSynPhotons(this->photons, this->shock, this->electrons);

    this->obs.changeViewingAngle(theta_v);

    for (auto& data : obs_data.light_curve) {
        data.Fv_model = obs.specificFlux(data.t, data.nu, this->photons);
    }

    for (auto& data : obs_data.spectrum) {
        data.Fv_model = obs.spectrum(data.nu, data.t, this->photons);
    }
    /*auto coord = autoGrid(jet, t_bins, theta_w, theta_v, phi_num, theta_num, t_num);

    auto shock = genForwardShock(coord, medium, jet, eps_e, eps_B, config.rtol);

    auto electrons = genSynElectrons(shock, p, xi);

    auto photons = genSynPhotons(shock, electrons);

    Observer obs(coord, shock, theta_v, lumi_dist, z);

    for (auto& data : obs_data.light_curve) {
        data.Fv_model = obs.specificFlux(data.t, data.nu, photons);
    }

    for (auto& data : obs_data.spectrum) {
        data.Fv_model = obs.spectrum(data.nu, data.t, photons);
    }*/

    return obs_data.calcChiSquare();
}