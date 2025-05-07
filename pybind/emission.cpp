#include "emission.h"

// Light curve
Emission::Emission(Real nu, const std::vector<Real>& t)
    : nu(xt::adapt(std::vector<Real>{nu}) * unit::Hz),
      t(xt::adapt(t) * unit::sec),
      F_nu(xt::adapt(std::vector<Real>(t.size(), 0.0))) {}

// Spectrum
Emission::Emission(const std::vector<Real>& nu, Real t)
    : t(xt::adapt(std::vector<Real>{t}) * unit::sec),
      nu(xt::adapt(nu) * unit::Hz),
      F_nu(xt::adapt(std::vector<Real>(nu.size(), 0.0))) {}
      

void Emission::generate(const Params& param, const ConfigParams& config) {
    Observer obs;
    SynElectronGrid electrons;
    SynPhotonGrid photons;

    Medium medium;
    if (config.medium == "ism") {
        std::tie(medium.rho, medium.mass) = evn::ISM(param.n_ism / unit::cm3);
    } else if (config.medium == "wind") {
        std::tie(medium.rho, medium.mass) = evn::wind(param.A_star);
    } else {
        std::cerr << "Unknown medium: " << config.medium << std::endl;
        return;
    }

    Ejecta jet;
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

    auto coord = auto_grid(jet, t, param.theta_w, param.theta_v, config.z, config.phi_grid, config.theta_grid, config.t_grid);

    auto shock = generate_fwd_shock(coord, medium, jet, param.eps_e, param.eps_B, config.rtol);

    obs.observe(coord, shock, config.lumi_dist * unit::cm, config.z);

    electrons = generate_syn_electrons(shock, param.p, param.xi);
    photons = generate_syn_photons(shock, electrons);

    if (t.size() > 1 && nu.size() == 1) {
        F_nu = xt::view(obs.specific_flux(t, nu, photons), 0, xt::all()) / unit::Jy;
    } else {
        F_nu = xt::view(obs.spectra(nu, t, photons), 0, xt::all()) / unit::Jy;
    }
}