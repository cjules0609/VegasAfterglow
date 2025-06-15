//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "pymodel.h"

#include "afterglow.h"

Ejecta PyTophatJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading, Real T0, std::optional<PyMagnetar> magnetar) {
    Ejecta jet;
    jet.eps_k = math::tophat(theta_c, E_iso * unit::erg / (4 * con::pi));
    jet.Gamma0 = math::tophat(theta_c, Gamma0);
    jet.spreading = spreading;
    jet.T0 = T0 * unit::sec;

    if (magnetar) {
        jet.deps_dt = [=](Real phi, Real theta, Real t) {
            if (theta <= theta_c) {
                Real tt = 1 + t / magnetar->t_0;
                return magnetar->L_0 / (tt * tt);
            } else {
                return 0.;
            }
        };
    }

    return jet;
}

Ejecta PyGaussianJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading, Real T0,
                     std::optional<PyMagnetar> magnetar) {
    Ejecta jet;
    jet.eps_k = math::gaussian(theta_c, E_iso * unit::erg / (4 * con::pi));
    jet.Gamma0 = math::gaussian(theta_c, Gamma0);
    jet.spreading = spreading;
    jet.T0 = T0 * unit::sec;

    if (magnetar) {
        jet.deps_dt = [=](Real phi, Real theta, Real t) {
            if (theta <= theta_c) {
                Real tt = 1 + t / magnetar->t_0;
                return magnetar->L_0 / (tt * tt);
            } else {
                return 0.;
            }
        };
    }

    return jet;
}

Ejecta PyPowerLawJet(Real theta_c, Real E_iso, Real Gamma0, Real k, bool spreading, Real T0,
                     std::optional<PyMagnetar> magnetar) {
    Ejecta jet;
    jet.eps_k = math::powerlaw(theta_c, E_iso * unit::erg / (4 * con::pi), k);
    jet.Gamma0 = math::powerlaw(theta_c, Gamma0, k);
    jet.spreading = spreading;
    jet.T0 = T0 * unit::sec;

    if (magnetar) {
        jet.deps_dt = [=](Real phi, Real theta, Real t) {
            if (theta <= theta_c) {
                Real tt = 1 + t / magnetar->t_0;
                return magnetar->L_0 / (tt * tt);
            } else {
                return 0.;
            }
        };
    }

    return jet;
}

Medium PyISM(Real n_ism) {
    Medium medium;
    std::tie(medium.rho, medium.mass) = evn::ISM(n_ism / unit::cm3);
    return medium;
}

Medium PyWind(Real A_star) {
    Medium medium;
    std::tie(medium.rho, medium.mass) = evn::wind(A_star);
    return medium;
}

void PyModel::single_shock_emission(Shock const& shock, Coord const& coord, Array const& t_obs, Array const& nu_obs,
                                    Observer& obs, PyRadiation rad, FluxDict& flux_dict, std::string suffix,
                                    bool return_trace) {
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

        if (return_trace) {
            flux_dict["IC" + suffix] = obs.point_specific_flux(t_obs, nu_obs, IC_ph) / unit::flux_den_cgs;
        } else {
            flux_dict["IC" + suffix] = obs.specific_flux(t_obs, nu_obs, IC_ph) / unit::flux_den_cgs;
        }
    }
    if (return_trace) {
        flux_dict["syn" + suffix] = obs.point_specific_flux(t_obs, nu_obs, syn_ph) / unit::flux_den_cgs;
    } else {
        flux_dict["syn" + suffix] = obs.specific_flux(t_obs, nu_obs, syn_ph) / unit::flux_den_cgs;
    }
}

auto PyModel::compute_specific_flux(Array const& t_obs, Array const& nu_obs, bool return_trace) -> FluxDict {
    Coord coord = auto_grid(jet, t_obs, this->theta_w, obs_setup.theta_obs, obs_setup.z, phi_resol, theta_resol,
                            t_resol, axisymmetric);

    FluxDict flux_dict;

    Observer observer;

    if (!rvs_rad_opt) {
        auto fwd_shock = generate_fwd_shock(coord, medium, jet, fwd_rad.rad, rtol);

        single_shock_emission(fwd_shock, coord, t_obs, nu_obs, observer, fwd_rad, flux_dict, "", return_trace);

        return flux_dict;
    } else {
        auto rvs_rad = *rvs_rad_opt;
        auto [fwd_shock, rvs_shock] = generate_shock_pair(coord, medium, jet, fwd_rad.rad, rvs_rad.rad, rtol);

        single_shock_emission(fwd_shock, coord, t_obs, nu_obs, observer, fwd_rad, flux_dict, "", return_trace);

        single_shock_emission(rvs_shock, coord, t_obs, nu_obs, observer, rvs_rad, flux_dict, "_rvs", return_trace);

        return flux_dict;
    }
}

auto PyModel::specific_flux(PyArray const& t, PyArray const& nu) -> FluxDict {
    Array t_obs = t * unit::sec;
    Array nu_obs = nu * unit::Hz;
    bool return_trace = true;

    if (t_obs.size() != nu_obs.size()) {
        throw std::invalid_argument(
            "time and frequency arrays must have the same size\n"
            "If you intend to get matrix output, use `specific_flux_matrix` instead");
    }

    return compute_specific_flux(t_obs, nu_obs, return_trace);
}

auto PyModel::specific_flux_matrix(PyArray const& t, PyArray const& nu) -> FluxDict {
    Array t_obs = t * unit::sec;
    Array nu_obs = nu * unit::Hz;
    bool return_trace = false;

    return compute_specific_flux(t_obs, nu_obs, return_trace);
}
