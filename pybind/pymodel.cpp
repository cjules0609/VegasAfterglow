//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "pymodel.h"

#include <algorithm>
#include <numeric>

#include "afterglow.h"
#include "error_handling.h"
#include "xtensor/misc/xsort.hpp"

Ejecta PyTophatJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading, Real duration,
                   std::optional<PyMagnetar> magnetar) {
    Ejecta jet;
    jet.eps_k = math::tophat(theta_c, E_iso);
    jet.Gamma0 = math::tophat_plus_one(theta_c, Gamma0 - 1);
    jet.spreading = spreading;
    jet.T0 = duration;

    if (magnetar) {
        jet.deps_dt = math::magnetar_injection(magnetar->t0, magnetar->q, magnetar->L0, theta_c);
    }

    return jet;
}

Ejecta PyGaussianJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading, Real duration,
                     std::optional<PyMagnetar> magnetar) {
    Ejecta jet;
    jet.eps_k = math::gaussian(theta_c, E_iso);
    jet.Gamma0 = math::gaussian_plus_one(theta_c, Gamma0 - 1);
    jet.spreading = spreading;
    jet.T0 = duration;

    if (magnetar) {
        jet.deps_dt = math::magnetar_injection(magnetar->t0, magnetar->q, magnetar->L0, theta_c);
    }

    return jet;
}

Ejecta PyPowerLawJet(Real theta_c, Real E_iso, Real Gamma0, Real k_e, Real k_g, bool spreading, Real duration,
                     std::optional<PyMagnetar> magnetar) {
    Ejecta jet;
    jet.eps_k = math::powerlaw(theta_c, E_iso, k_e);
    jet.Gamma0 = math::powerlaw_plus_one(theta_c, Gamma0 - 1, k_g);
    jet.spreading = spreading;
    jet.T0 = duration;

    if (magnetar) {
        jet.deps_dt = math::magnetar_injection(magnetar->t0, magnetar->q, magnetar->L0, theta_c);
    }

    return jet;
}

Ejecta PyPowerLawWing(Real theta_c, Real E_w, Real Gamma_w, Real k_e, Real k_g, bool spreading, Real duration) {
    Ejecta jet;
    jet.eps_k = math::powerlaw_wing(theta_c, E_w, k_e);
    jet.Gamma0 = math::powerlaw_wing_plus_one(theta_c, Gamma_w - 1, k_g);
    jet.spreading = spreading;
    jet.T0 = duration;

    return jet;
}

Ejecta PyStepPowerLawJet(Real theta_c, Real E_c, Real Gamma_c, Real E_w, Real Gamma_w, Real k_e, Real k_g,
                         bool spreading, Real duration, std::optional<PyMagnetar> magnetar) {
    Ejecta jet;
    jet.eps_k = math::step_powerlaw(theta_c, E_c, E_w, k_e);
    jet.Gamma0 = math::step_powerlaw_plus_one(theta_c, Gamma_c - 1, Gamma_w - 1, k_g);

    jet.spreading = spreading;
    jet.T0 = duration;

    if (magnetar) {
        jet.deps_dt = math::magnetar_injection(magnetar->t0, magnetar->q, magnetar->L0, theta_c);
    }

    return jet;
}

Ejecta PyTwoComponentJet(Real theta_c, Real E_iso_c, Real Gamma0_c, Real theta_w, Real E_iso_w, Real Gamma0_w,
                         bool spreading, Real duration, std::optional<PyMagnetar> magnetar) {
    Ejecta jet;
    jet.eps_k = math::two_component(theta_c, theta_w, E_iso_c, E_iso_w);

    jet.Gamma0 = math::two_component_plus_one(theta_c, theta_w, Gamma0_c - 1, Gamma0_w - 1);

    jet.spreading = spreading;
    jet.T0 = duration;

    if (magnetar) {
        jet.deps_dt = math::magnetar_injection(magnetar->t0, magnetar->q, magnetar->L0, theta_c);
    }

    return jet;
}

Medium PyISM(Real n_ism) {
    Medium medium;
    medium.rho = [=](Real phi, Real theta, Real r) { return n_ism * 1.67e-24; };

    return medium;
}

Medium PyWind(Real A_star, Real n_ism, Real n_0) {
    Medium medium;

    Real rho_ism = n_ism * 1.67e-24;
    Real r02 = A_star * 5e11 / (n_0 * 1.67e-24);

    medium.rho = [=](Real phi, Real theta, Real r) { return A_star * 5e11 / (r02 + r * r) + rho_ism; };
    return medium;
}

void convert_unit(Ejecta& jet, Medium& medium) {
    auto eps_k_cgs = jet.eps_k;
    jet.eps_k = [=](Real phi, Real theta) { return eps_k_cgs(phi, theta) * (unit::erg / (4 * con::pi)); };

    auto deps_dt_cgs = jet.deps_dt;
    jet.deps_dt = [=](Real phi, Real theta, Real t) {
        return deps_dt_cgs(phi, theta, t / unit::sec) * (unit::erg / (4 * con::pi * unit::sec));
    };

    auto dm_dt_cgs = jet.dm_dt;
    jet.dm_dt = [=](Real phi, Real theta, Real t) {
        return dm_dt_cgs(phi, theta, t / unit::sec) * (unit::g / (4 * con::pi * unit::sec));
    };

    jet.T0 *= unit::sec;

    auto rho_cgs = medium.rho; // number density from python side
    medium.rho = [=](Real phi, Real theta, Real r) {
        return rho_cgs(phi, theta, r / unit::cm) * (unit::g / unit::cm3); // convert to density
    };
}

void PyModel::single_shock_emission(Shock const& shock, Coord const& coord, Array const& t_obs, Array const& nu_obs,
                                    Observer& obs, PyRadiation rad, ArrayDict& flux_dict, std::string suffix,
                                    bool serilized) {
    obs.observe(coord, shock, obs_setup.lumi_dist, obs_setup.z);

    auto syn_e = generate_syn_electrons(shock);

    auto syn_ph = generate_syn_photons(shock, syn_e);

    if (rad.IC_cooling) {
        if (rad.KN) {
            KN_cooling(syn_e, syn_ph, shock);
        } else {
            Thomson_cooling(syn_e, syn_ph, shock);
        }
    }

    if (rad.SSC) {
        auto IC_ph = generate_IC_photons(syn_e, syn_ph, rad.KN);

        if (serilized) {
            flux_dict["IC" + suffix] = obs.specific_flux_series(t_obs, nu_obs, IC_ph) / unit::flux_den_cgs;
        } else {
            flux_dict["IC" + suffix] = obs.specific_flux(t_obs, nu_obs, IC_ph) / unit::flux_den_cgs;
        }
    }

    if (serilized) {
        flux_dict["syn" + suffix] = obs.specific_flux_series(t_obs, nu_obs, syn_ph) / unit::flux_den_cgs;
    } else {
        flux_dict["syn" + suffix] = obs.specific_flux(t_obs, nu_obs, syn_ph) / unit::flux_den_cgs;
    }
}

auto PyModel::compute_specific_flux(Array const& t_obs, Array const& nu_obs, bool serilized) -> ArrayDict {
    Coord coord = auto_grid(jet, t_obs, this->theta_w, obs_setup.theta_obs, obs_setup.z, phi_resol, theta_resol,
                            t_resol, axisymmetric);

    ArrayDict flux_dict;

    Observer observer;

    if (!rvs_rad_opt) {
        auto fwd_shock = generate_fwd_shock(coord, medium, jet, fwd_rad.rad, rtol);

        single_shock_emission(fwd_shock, coord, t_obs, nu_obs, observer, fwd_rad, flux_dict, "", serilized);

        return flux_dict;
    } else {
        auto rvs_rad = *rvs_rad_opt;
        auto [fwd_shock, rvs_shock] = generate_shock_pair(coord, medium, jet, fwd_rad.rad, rvs_rad.rad, rtol);

        single_shock_emission(fwd_shock, coord, t_obs, nu_obs, observer, fwd_rad, flux_dict, "", serilized);

        single_shock_emission(rvs_shock, coord, t_obs, nu_obs, observer, rvs_rad, flux_dict, "_rvs", serilized);

        return flux_dict;
    }
}

void save_shock_details(Shock const& shock, ArrayDict& detail_dict, std::string suffix) {
    detail_dict["Gamma" + suffix] = shock.Gamma;
    detail_dict["Gamma_th" + suffix] = shock.Gamma_th;
    detail_dict["r" + suffix] = shock.r / unit::cm;
    detail_dict["t_comv" + suffix] = shock.t_comv / unit::sec;
    detail_dict["B" + suffix] = shock.B / unit::Gauss;
    detail_dict["N_p" + suffix] = shock.N_p;
    detail_dict["theta" + suffix] = shock.theta;
}

template <typename ElectronGrid>
void save_electron_details(ElectronGrid const& electrons, ArrayDict& detail_dict, std::string suffix) {
    auto shape = electrons.shape();
    detail_dict["gamma_m" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    detail_dict["gamma_c" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    detail_dict["gamma_a" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    detail_dict["gamma_M" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    detail_dict["N_e" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    for (size_t i = 0; i < shape[0]; ++i) {
        for (size_t j = 0; j < shape[1]; ++j) {
            for (size_t k = 0; k < shape[2]; ++k) {
                detail_dict["gamma_a" + suffix](i, j, k) = electrons(i, j, k).gamma_a;
                detail_dict["gamma_m" + suffix](i, j, k) = electrons(i, j, k).gamma_m;
                detail_dict["gamma_c" + suffix](i, j, k) = electrons(i, j, k).gamma_c;
                detail_dict["gamma_M" + suffix](i, j, k) = electrons(i, j, k).gamma_M;
                detail_dict["N_e" + suffix](i, j, k) = electrons(i, j, k).N_e;
            }
        }
    }
}
template <typename PhotonGrid>
void save_photon_details(PhotonGrid const& photons, ArrayDict& detail_dict, std::string suffix) {
    auto shape = photons.shape();
    detail_dict["nu_m" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    detail_dict["nu_c" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    detail_dict["nu_a" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    detail_dict["nu_M" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    detail_dict["I_nu_max" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    for (size_t i = 0; i < shape[0]; ++i) {
        for (size_t j = 0; j < shape[1]; ++j) {
            for (size_t k = 0; k < shape[2]; ++k) {
                detail_dict["nu_a" + suffix](i, j, k) = photons(i, j, k).nu_a;
                detail_dict["nu_m" + suffix](i, j, k) = photons(i, j, k).nu_m;
                detail_dict["nu_c" + suffix](i, j, k) = photons(i, j, k).nu_c;
                detail_dict["nu_M" + suffix](i, j, k) = photons(i, j, k).nu_M;
                detail_dict["I_nu_max" + suffix](i, j, k) = photons(i, j, k).I_nu_max;
            }
        }
    }
    detail_dict["nu_m" + suffix] /= unit::Hz;
    detail_dict["nu_c" + suffix] /= unit::Hz;
    detail_dict["nu_a" + suffix] /= unit::Hz;
    detail_dict["nu_M" + suffix] /= unit::Hz;
    detail_dict["I_nu_max" + suffix] /= (unit::erg / (unit::Hz * unit::sec * unit::cm2));
}

void PyModel::single_evo_details(Shock const& shock, Coord const& coord, Array const& t_obs, Observer& obs,
                                 PyRadiation rad, ArrayDict& detail_dict, std::string suffix) {
    obs.observe(coord, shock, obs_setup.lumi_dist, obs_setup.z);

    detail_dict["t_obs" + suffix] = obs.time / unit::sec;
    detail_dict["Doppler" + suffix] = xt::exp2(obs.lg2_doppler);
    // detail_dict["Omega" + suffix] = xt::exp2(obs.lg2_emission_area);

    auto syn_e = generate_syn_electrons(shock);

    auto syn_ph = generate_syn_photons(shock, syn_e);

    if (rad.IC_cooling) {
        if (rad.KN) {
            KN_cooling(syn_e, syn_ph, shock);
        } else {
            Thomson_cooling(syn_e, syn_ph, shock);
        }
    }
    save_electron_details(syn_e, detail_dict, suffix);
    save_photon_details(syn_ph, detail_dict, suffix);
}

auto PyModel::details(Real t_min, Real t_max) -> ArrayDict {
    Array t_obs = xt::logspace(std::log10(t_min * unit::sec), std::log10(t_max * unit::sec), 10);
    Coord coord = auto_grid(jet, t_obs, this->theta_w, obs_setup.theta_obs, obs_setup.z, phi_resol, theta_resol,
                            t_resol, axisymmetric);

    ArrayDict details_dict;

    details_dict["phi"] = coord.phi;
    details_dict["theta"] = coord.theta;
    details_dict["t_src"] = coord.t / unit::sec;

    Observer observer;

    if (!rvs_rad_opt) {
        auto fwd_shock = generate_fwd_shock(coord, medium, jet, fwd_rad.rad, rtol);

        save_shock_details(fwd_shock, details_dict, "_fwd");

        single_evo_details(fwd_shock, coord, t_obs, observer, fwd_rad, details_dict, "_fwd");

        return details_dict;
    } else {
        auto rvs_rad = *rvs_rad_opt;
        auto [fwd_shock, rvs_shock] = generate_shock_pair(coord, medium, jet, fwd_rad.rad, rvs_rad.rad, rtol);

        save_shock_details(fwd_shock, details_dict, "_fwd");

        save_shock_details(rvs_shock, details_dict, "_rvs");

        single_evo_details(fwd_shock, coord, t_obs, observer, fwd_rad, details_dict, "_fwd");

        single_evo_details(rvs_shock, coord, t_obs, observer, rvs_rad, details_dict, "_rvs");

        return details_dict;
    }
}

auto PyModel::specific_flux_series(PyArray const& t, PyArray const& nu) -> ArrayDict {
    AFTERGLOW_REQUIRE(
        t.size() == nu.size(),
        "time and frequency arrays must have the same size\nIf you intend to get matrix-like output, use the "
        "generic `specific_flux` instead");
    AFTERGLOW_REQUIRE(is_ascending(t), "time array must be in ascending order");

    Array t_obs = t * unit::sec;
    Array nu_obs = nu * unit::Hz;
    bool serialized = true;

    return compute_specific_flux(t_obs, nu_obs, serialized);
}

auto PyModel::specific_flux_series_with_expo(PyArray const& t, PyArray const& nu, PyArray const& expo_time,
                                             size_t num_points) -> ArrayDict {
    AFTERGLOW_REQUIRE(t.size() == nu.size() && t.size() == expo_time.size(),
                      "time, frequency, and exposure time arrays must have the same size");
    AFTERGLOW_REQUIRE(num_points >= 2, "num_points must be at least 2 to sample within each exposure time");

    size_t total_points = t.size() * num_points;
    Array t_obs = Array::from_shape({total_points});
    Array nu_obs = Array::from_shape({total_points});
    std::vector<size_t> idx(total_points);

    for (size_t i = 0, j = 0; i < t.size() && j < total_points; ++i) {
        Real t_start = t(i);
        Real dt = expo_time(i) / (num_points - 1);

        for (size_t k = 0; k < num_points && j < total_points; ++k, ++j) {
            t_obs(j) = t_start + k * dt;
            nu_obs(j) = nu(i);
            idx[j] = i;
        }
    }

    std::vector<size_t> sort_indices(total_points);
    std::iota(sort_indices.begin(), sort_indices.end(), 0);
    std::sort(sort_indices.begin(), sort_indices.end(), [&t_obs](size_t i, size_t j) { return t_obs(i) < t_obs(j); });

    Array t_obs_sorted = Array::from_shape({total_points});
    Array nu_obs_sorted = Array::from_shape({total_points});
    std::vector<size_t> idx_sorted(idx.size());

    for (size_t i = 0; i < sort_indices.size(); ++i) {
        size_t orig_idx = sort_indices[i];
        t_obs_sorted(i) = t_obs(orig_idx);
        nu_obs_sorted(i) = nu_obs(orig_idx);
        idx_sorted[i] = idx[orig_idx];
    }

    t_obs_sorted *= unit::sec;
    nu_obs_sorted *= unit::Hz;

    bool serialized = true;

    auto result = compute_specific_flux(t_obs_sorted, nu_obs_sorted, serialized);

    for (auto& [key, val] : result) {
        Array summed = xt::zeros<Real>({t.size()});

        for (size_t j = 0; j < val.size(); j++) {
            size_t orig_time_idx = idx_sorted[j];
            summed(orig_time_idx) += val(j);
        }

        summed /= static_cast<Real>(num_points);
        result[key] = summed;
    }
    return result;
}

auto PyModel::specific_flux(PyArray const& t, PyArray const& nu) -> ArrayDict {
    AFTERGLOW_REQUIRE(is_ascending(t), "time array must be in ascending order");

    Array t_obs = t * unit::sec;
    Array nu_obs = nu * unit::Hz;
    bool return_trace = false;

    return compute_specific_flux(t_obs, nu_obs, return_trace);
}
