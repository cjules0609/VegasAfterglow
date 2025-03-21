
#include <boost/numeric/odeint.hpp>
#include <filesystem>
#include <fstream>

#include "afterglow.h"

void GCN36236(Real theta_c) {
    Real E_iso = 1e53 * con::erg * (0.088 / theta_c) * (0.088 / theta_c);
    Real lumi_dist = 1.2e27 * con::cm;
    Real z = 0.01;
    Real theta_w = 0.6;
    Real Gamma0 = 300;
    Real n_ism = 0.0199526231496888 / con::cm3;
    Real eps_e = 0.01;
    Real eps_B = 0.00019952623149688788;
    Real p = 2.139;

    // create model
    Medium medium;

    medium.rho = evn::ISM(n_ism);

    Ejecta jet;
    jet.dE0dOmega = math::gaussian(theta_c, E_iso / (4 * con::pi));
    jet.Gamma0 = math::gaussian(theta_c, Gamma0);

    size_t r_num = 64;
    size_t theta_num = 64;
    size_t phi_num = 64;

    Array t_bins = logspace(1e-1 * con::day, 1e3 * con::day, 100);

    Coord coord = autoGrid(jet, t_bins, theta_w, phi_num, theta_num, r_num);

    Shock f_shock = genForwardShock(coord, medium, jet, eps_e, eps_B);

    auto syn_e = genSynElectrons(f_shock, p);

    auto syn_ph = genSynPhotons(f_shock, syn_e);

    Array theta_obs = linspace(5 * con::deg, 60 * con::deg, 56);

    Observer obs(coord, f_shock, 0, lumi_dist, z);

    std::vector<std::string> prefix = {"Optical"};

    Array band_pass = logspace(eVtoHz(1.5 * con::eV), eVtoHz(3.5 * con::eV), 5);
    for (auto theta_v : theta_obs) {
        obs.changeViewingAngle(theta_v);

        auto to_suffix = [](Real theta) { return std::to_string(int((theta / con::deg))); };

        char buff[100] = {0};

        sprintf(buff, "data/F_nu_syn_%.0lf_%.0lf_", (theta_v / con::deg), (theta_c / con::deg));

        std::string fname = buff + prefix[0];

        std::cout << fname << std::endl;

        // Array F_syn = obs.flux(t_bins, band_pass, syn_ph);
        Array F_syn = obs.specificFlux(t_bins, eVtoHz(2.25 * con::eV), syn_ph);

        output(F_syn, fname, con::erg / con::sec / con::cm / con::cm / con::Hz);
    }

    output(t_bins, "data/t_obs", con::sec);
}

void std_afterglow(Real theta_c) {
    Real E_iso = 1e53 * con::erg * (0.088 / theta_c) * (0.088 / theta_c);
    Real lumi_dist = 3e26 * con::cm;
    Real z = 0.01;
    Real theta_w = 0.6;
    std::vector<Real> theta_obs{0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};
    // 0.6;
    Real Gamma0 = 300;
    Real n_ism = 0.0199526231496888 / con::cm3;
    Real eps_e = 0.01;
    Real eps_B = 0.00019952623149688788;
    Real p = 2.139;

    // create model
    Medium medium;
    medium.rho = evn::ISM(n_ism);

    Ejecta jet;
    jet.dE0dOmega = math::gaussian(theta_c, E_iso / (4 * con::pi));
    jet.Gamma0 = math::gaussian(theta_c, Gamma0);

    size_t r_num = 256;
    size_t theta_num = 64;
    size_t phi_num = 64;

    Array t_bins = logspace(1e-1 * con::day, 3e3 * con::day, 100);

    Real theta_obs_max = (*std::max_element(theta_obs.begin(), theta_obs.end())) * con::deg;

    Coord coord = autoGrid(jet, t_bins, theta_w, theta_obs_max, phi_num, theta_num, r_num);

    Shock f_shock = genForwardShock(coord, medium, jet, eps_e, eps_B);

    // output(f_shock, "shock");

    auto syn_e = genSynElectrons(f_shock, p);

    // output(syn_e, "syn_e");

    auto syn_ph = genSynPhotons(f_shock, syn_e);
    // output(syn_ph, "syn_ph");

    Observer obs(coord, f_shock, 0, lumi_dist, z);
    // output(obs.t_obs_grid, "t_obs", con::sec);
    // output(obs.doppler, "doppler");

    std::array<double, 5> lo = {eVtoHz(0.3 * con::keV)};

    std::array<double, 5> hi = {eVtoHz(10 * con::keV)};

    std::vector<std::string> prefix = {"EP"};

    Array band_pass = logspace(lo[0], hi[0], 10);

    for (auto theta_v : theta_obs) {
        obs.changeViewingAngle(theta_v * con::deg);

        char buff[100] = {0};

        sprintf(buff, "data/F_std_%.1lf_%.1lf", std::floor(theta_c / con::deg), (theta_v));

        std::string fname = buff + prefix[0];

        std::cout << fname << std::endl;

        Array F_syn = obs.flux(t_bins, band_pass, syn_ph);

        output(F_syn, fname, con::erg / con::sec / con::cm / con::cm);
    }
    output(t_bins, "data/t_std", con::sec);
}

int main() {
    Real theta_c[] = {5 * con::deg, 7 * con::deg, 10 * con::deg};

    for (auto theta : theta_c) {
        GCN36236(theta);
    }

    // std_afterglow(5 * con::deg);
    // std_afterglow(10 * con::deg);
    // std_afterglow(15.01 * con::deg);

    return 0;
}