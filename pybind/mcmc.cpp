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

#include "afterglow.h"
void sort_synchronized(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c) {
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

double SpectrumData::calcChiSquare() const {
    double chi_square = 0;
    for (size_t i = 0; i < nu.size(); ++i) {
        if (Fv_err[i] == 0) continue;
        double diff = Fv_obs[i] - Fv_model[i];
        chi_square += (diff * diff) / (Fv_err[i] * Fv_err[i]);
    }
    return chi_square;
}

void MultiBandData::genEstimatePoints() {
    t_grid.clear();
    lc_band.clear();
    t_grid.reserve(light_curve.size() * 10);
    lc_band.reserve(light_curve.size());
    for (auto& data : light_curve) {
        t_grid.insert(t_grid.end(), data.t.begin(), data.t.end());
        lc_band.push_back(data.nu);
        data.Fv_model.resize(data.t.size());
    }
    std::sort(t_grid.begin(), t_grid.end());
    t_grid.erase(std::unique(t_grid.begin(), t_grid.end()), t_grid.end());

    nu_grid.clear();
    spectrum_t.clear();
    nu_grid.reserve(spectrum.size() * 10);
    spectrum_t.reserve(spectrum.size());
    for (auto& data : spectrum) {
        nu_grid.insert(nu_grid.end(), data.nu.begin(), data.nu.end());
        spectrum_t.push_back(data.t);
        data.Fv_model.resize(data.nu.size());
    }
    std::sort(nu_grid.begin(), nu_grid.end());
    nu_grid.erase(std::unique(nu_grid.begin(), nu_grid.end()), nu_grid.end());

    if (t_grid.empty()) {
        std::cerr << "Error: No light curve data loaded!\n" << std::endl;
        std::exit(EXIT_FAILURE);
    }
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

void MultiBandData::addModelLightCurve(MeshGrid const& specific_flux) {
    for (size_t i = 0; i < light_curve.size(); ++i) {
        auto& data = light_curve[i];
        for (size_t j = 0, k = 0; j < t_grid.size() && k < data.t.size(); ++j) {
            if (t_grid[j] == data.t[k]) {
                data.Fv_model[k] = specific_flux[i][j];
                ++k;
            }
        }
    }
}

void MultiBandData::addModelSpectrum(MeshGrid const& specific_flux) {
    for (size_t i = 0; i < spectrum.size(); ++i) {
        auto& data = spectrum[i];
        for (size_t j = 0, k = 0; j < nu_grid.size() && k < data.nu.size(); ++j) {
            if (nu_grid[j] == data.nu[k]) {
                data.Fv_model[k] = specific_flux[i][j];
                ++k;
            }
        }
    }
}

void MultiBandData::addObsLightCurve(double nu, List const& t, List const& Fv_obs, List const& Fv_err) {
    assert(t.size() == Fv_obs.size() && t.size() == Fv_err.size() && "light curve array inconsistent length!");
    LightCurveData data;
    data.nu = nu * con::Hz;
    data.t = t;
    data.Fv_obs = Fv_obs;
    data.Fv_err = Fv_err;
    sort_synchronized(data.t, data.Fv_obs, data.Fv_err);

    for (auto& tt : data.t) {
        tt *= con::sec;
    }
    for (auto& f : data.Fv_obs) {
        f *= con::erg / con::sec / con::cm2 / con::Hz;
    }
    for (auto& f : data.Fv_err) {
        f *= con::erg / con::sec / con::cm2 / con::Hz;
    }
    data.Fv_model.resize(t.size());
    light_curve.push_back(std::move(data));
}
void MultiBandData::addObsSpectrum(double t, List const& nu, List const& Fv_obs, List const& Fv_err) {
    assert(nu.size() == Fv_obs.size() && nu.size() == Fv_err.size() && "spectrum array inconsistent length!");
    SpectrumData data;
    data.t = t * con::sec;
    data.nu = nu;
    data.Fv_obs = Fv_obs;
    data.Fv_err = Fv_err;
    sort_synchronized(data.nu, data.Fv_obs, data.Fv_err);

    for (auto& v : data.nu) {
        v *= con::Hz;
    }
    for (auto& f : data.Fv_obs) {
        f *= con::erg / con::sec / con::cm2 / con::Hz;
    }
    for (auto& f : data.Fv_err) {
        f *= con::erg / con::sec / con::cm2 / con::Hz;
    }
    data.Fv_model.resize(nu.size());
    spectrum.push_back(std::move(data));
}

void MultiBandModel::configure(ConfigParams const& param) { this->config = param; }

std::vector<double> MultiBandModel::chiSquareBatch(std::vector<Params> const& param_batch) {
    const size_t N = param_batch.size();
    std::vector<double> results(N);

#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        results[i] = chiSquare(param_batch[i]);
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

    Coord coord = autoGrid(jet, obs_data.t_grid, theta_w, theta_v, phi_num, theta_num, t_num);

    Shock f_shock = genForwardShock(coord, medium, jet, eps_e, eps_B, config.rtol);

    auto syn_e = genSynElectrons(f_shock, p, xi);

    auto syn_ph = genSynPhotons(f_shock, syn_e);

    Observer obs(coord, f_shock, theta_v, lumi_dist, z);

    auto F_syn = obs.specificFlux(obs_data.t_grid, obs_data.lc_band, syn_ph);

    auto F_spec = obs.spectrum(obs_data.nu_grid, obs_data.spectrum_t, syn_ph);

    obs_data.addModelLightCurve(F_syn);

    obs_data.addModelSpectrum(F_spec);

    return obs_data.calcChiSquare();
}