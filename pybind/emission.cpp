#include "emission.h"

Emission::Emission(const Params& param, const ConfigParams& config) {
    this->param = param;
    this->config = config;

    if (config.medium == "ism") {
        std::tie(medium.rho, medium.mass) = evn::ISM(param.n_ism / unit::cm3);
    } else if (config.medium == "wind") {
        std::tie(medium.rho, medium.mass) = evn::wind(param.A_star);
    } else {
        std::cerr << "Unknown medium: " << config.medium << std::endl;
        return;
    }

    Real E_iso = param.E_iso * unit::erg;
    Real Gamma0 = param.Gamma0;

    if (config.jet == "tophat") {
        jet.eps_k = math::tophat(param.theta_c, E_iso / (4 * con::pi));
        jet.Gamma0 = math::tophat(param.theta_c, Gamma0);
    } else if (config.jet == "gaussian") {
        jet.eps_k = math::gaussian(param.theta_c, E_iso / (4 * con::pi));
        jet.Gamma0 = math::gaussian(param.theta_c, Gamma0);
    } else if (config.jet == "powerlaw") {
        jet.eps_k = math::powerlaw(param.theta_c, E_iso / (4 * con::pi), param.k_jet);
        jet.Gamma0 = math::powerlaw(param.theta_c, Gamma0, param.k_jet);
    } else {
        std::cerr << "Unknown jet structure: " << config.jet << std::endl;
        return;
    }
    jet.spreading = config.spreading;
}

void Emission::observe(const Params& param, const ConfigParams& config, const Array& t) {
    auto coord = auto_grid(jet, t, param.theta_w, param.theta_v, config.z, config.phi_grid, config.theta_grid, config.t_grid);
    // auto shock = generate_fwd_shock(coord, medium, jet, param.eps_e, param.eps_B, config.rtol);
    auto [f_shock, r_shock] = generate_shock_pair(coord, medium, jet, param.eps_e, param.eps_B, param.eps_e_r, param.eps_B_r);

    obs.observe(coord, f_shock, config.lumi_dist * unit::cm, config.z);
    obs.observe(coord, r_shock, config.lumi_dist * unit::cm, config.z);

    electrons = generate_syn_electrons(f_shock, param.p, param.xi);
    photons = generate_syn_photons(f_shock, electrons);

    electrons_r = generate_syn_electrons(r_shock, param.p_r, param.xi_r);
    photons_r = generate_syn_photons(r_shock, electrons);
}

std::vector<Real> Emission::lc(Real nu, const std::vector<Real>& t) {
    observe(param, config, xt::adapt(t) * unit::sec);  // Perform observation step
    auto Fv = obs.specific_flux(xt::adapt(t) * unit::sec, xt::adapt(std::vector<Real>{nu}) * unit::Hz, photons) / unit::Jy;
    return std::vector<Real>(Fv.begin(), Fv.end());
}

std::vector<Real> Emission::lc_r(Real nu, const std::vector<Real>& t) {
    observe(param, config, xt::adapt(t) * unit::sec);  // Perform observation step
    auto Fv_r = obs.specific_flux(xt::adapt(t) * unit::sec, xt::adapt(std::vector<Real>{nu}) * unit::Hz, photons_r) / unit::Jy;
    return std::vector<Real>(Fv_r.begin(), Fv_r.end());
}

std::vector<Real> Emission::spec(const std::vector<Real>& nu, Real t) {
    observe(param, config, xt::adapt(std::vector<Real>{t}) * unit::sec);  // Perform observation step
    auto Fv = obs.spectra(xt::adapt(nu) * unit::Hz, xt::adapt(std::vector<Real>{t}) * unit::sec, photons) / unit::Jy;
    
    return std::vector<Real>(Fv.begin(), Fv.end());
}