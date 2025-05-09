//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "pymodel.h"

#include "afterglow.h"

Ejecta PyTophatJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading, Real T0) {
    Ejecta jet;
    jet.eps_k = math::tophat(theta_c, E_iso * unit::erg / (4 * con::pi));
    jet.Gamma0 = math::tophat(theta_c, Gamma0);
    jet.spreading = spreading;
    jet.T0 = T0 * unit::sec;
    return jet;
}

Ejecta PyGaussianJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading, Real T0) {
    Ejecta jet;
    jet.eps_k = math::gaussian(theta_c, E_iso * unit::erg / (4 * con::pi));
    jet.Gamma0 = math::gaussian(theta_c, Gamma0);
    jet.spreading = spreading;
    jet.T0 = T0 * unit::sec;
    return jet;
}

Ejecta PyPowerLawJet(Real theta_c, Real E_iso, Real Gamma0, Real k, bool spreading, Real T0) {
    Ejecta jet;
    jet.eps_k = math::powerlaw(theta_c, E_iso * unit::erg / (4 * con::pi), k);
    jet.Gamma0 = math::powerlaw(theta_c, Gamma0, k);
    jet.spreading = spreading;
    jet.T0 = T0 * unit::sec;
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

void PyModel::specific_flux_for(Shock const& shock, Coord const& coord, Array const& t_obs, Array const& nu_obs,
                                Observer& obs, PyRadiation rad, FluxDict& flux_dict, std::string suffix) {
    obs.observe(coord, shock, obs_setup.lumi_dist, obs_setup.z);

    auto syn_e = generate_syn_electrons(shock, rad.p, rad.xi_e);

    auto syn_ph = generate_syn_photons(shock, syn_e);

    if (rad.SSC) {
        if (rad.Klein_Nishina) {
            KN_cooling(syn_e, syn_ph, shock);
        } else {
            Thomson_cooling(syn_e, syn_ph, shock);
        }

        auto IC_ph = gen_IC_photons(syn_e, syn_ph);

        flux_dict["IC" + suffix] = obs.specific_flux(t_obs, nu_obs, IC_ph) / unit::flux_den_cgs;
    }

    flux_dict["syn" + suffix] = obs.specific_flux(t_obs, nu_obs, syn_ph) / unit::flux_den_cgs;
}

auto PyModel::specific_flux_(Array const& t_obs, Array const& nu_obs) -> FluxDict {
    Coord coord =
        auto_grid(jet, t_obs, this->theta_w, obs_setup.theta_obs, obs_setup.z, phi_num, theta_num, t_num, axisymmetric);

    FluxDict flux_dict;

    Observer obs;

    if (!rvs_rad_opt) {
        auto f_shock = generate_fwd_shock(coord, medium, jet, fwd_rad.eps_e, fwd_rad.eps_B, rtol);

        specific_flux_for(f_shock, coord, t_obs, nu_obs, obs, fwd_rad, flux_dict, "");

        return flux_dict;
    } else {
        auto rvs_rad = *rvs_rad_opt;
        auto [f_shock, r_shock] =
            generate_shock_pair(coord, medium, jet, fwd_rad.eps_e, fwd_rad.eps_B, rvs_rad.eps_e, rvs_rad.eps_B, rtol);

        specific_flux_for(f_shock, coord, t_obs, nu_obs, obs, fwd_rad, flux_dict, "");

        specific_flux_for(r_shock, coord, t_obs, nu_obs, obs, rvs_rad, flux_dict, "_rvs");

        return flux_dict;
    }
}

auto PyModel::specific_flux(PyArray const& t, PyArray const& nu) -> FluxDict {
    Array t_obs = t * unit::sec;
    Array nu_obs = nu * unit::Hz;

    return specific_flux_(t_obs, nu_obs);
}

auto PyModel::spectra(PyArray const& nu, PyArray const& t) -> FluxDict {
    Array nu_obs = nu * unit::Hz;
    Array t_obs = t * unit::sec;

    FluxDict flux_dict = specific_flux_(t_obs, nu_obs);

    for (auto const& [key, value] : flux_dict) {
        flux_dict[key] = xt::transpose(value);
    }

    return flux_dict;
}
