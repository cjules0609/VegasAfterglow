#include <boost/numeric/odeint.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>

#include "afterglow.h"

void tests(size_t r_num, size_t theta_num, size_t phi_num, Real n_ism, Real eps_e, Real eps_B, Real p, Real E_iso,
           Real Gamma0, Real theta_c, Real theta_v, bool verbose = false) {
    Real z = 0.009;

    Real lumi_dist = 1.23e26 * con::cm;

    Array t_obs = xt::logspace(std::log10(1e3 * con::sec), std::log10(1e7 * con::sec), 50);
    // create model
    Medium medium;

    medium.rho = evn::ISM(n_ism);

    Ejecta jet;

    jet.dE0dOmega = math::gaussian(theta_c, E_iso / (4 * con::pi));

    jet.Gamma0 = math::gaussian(theta_c, Gamma0);

    Coord coord = autoGrid(jet, t_obs, 0.6, theta_v, phi_num, theta_num, r_num);

    // output(coord, "coord");

    Shock f_shock = genForwardShock(coord, medium, jet, eps_e, eps_B);

    // output(f_shock, "shock");

    auto syn_e = genSynElectrons(f_shock, p);

    auto syn_ph = genSynPhotons(f_shock, syn_e);

    Observer obs;

    obs.observe(coord, f_shock, theta_v, lumi_dist, z);

    // output(obs.t_obs_grid, "t_obs", con::sec);

    Real nu_obs = eVtoHz(1 * con::keV);

    Array F_nu = obs.specificFlux(t_obs, nu_obs, syn_ph);
    // Array F_nu = obs.flux(t_bins, linspace(eVtoHz(0.1 * con::keV), eVtoHz(10 * con::keV), 5), syn_ph);

    if (verbose) {
        output(t_obs, "t_c", con::sec);
        output(F_nu, "F_nu" + std::to_string(phi_num) + "-" + std::to_string(theta_num) + "-" + std::to_string(r_num),
               (con::erg / con::cm / con::cm / con::sec));
    }

    return;
}

int main() {
    size_t r_num = 100;
    size_t theta_num = 100;
    size_t phi_num = 100;

    Real n_ism = 2 / con::cm3;
    Real eps_e = 1e-2;
    Real eps_B = 1e-4;
    Real p = 2.1;
    Real Gamma0 = 300;

    Array E_iso = xt::logspace(std::log10(1e48 * con::erg), std::log10(1e52 * con::erg), 100);
    Array theta_c = xt::linspace(0.01, 0.1, 100);
    Array theta_v = xt::linspace(0.01, 0.5, 5);

    tests(256, 256, 256, n_ism, eps_e, eps_B, p, 1e52 * con::erg, Gamma0, 0.1, 0.3, true);
    tests(128, 128, 128, n_ism, eps_e, eps_B, p, 1e52 * con::erg, Gamma0, 0.1, 0.3, true);
    tests(64, 64, 64, n_ism, eps_e, eps_B, p, 1e52 * con::erg, Gamma0, 0.1, 0.3, true);
    tests(32, 32, 32, n_ism, eps_e, eps_B, p, 1e52 * con::erg, Gamma0, 0.1, 0.3, true);
    tests(30, 30, 30, n_ism, eps_e, eps_B, p, 1e52 * con::erg, Gamma0, 0.1, 0.3, true);
    tests(28, 28, 28, n_ism, eps_e, eps_B, p, 1e52 * con::erg, Gamma0, 0.1, 0.3, true);
    tests(16, 16, 16, n_ism, eps_e, eps_B, p, 1e52 * con::erg, Gamma0, 0.1, 0.3, true);
    tests(8, 8, 8, n_ism, eps_e, eps_B, p, 1e52 * con::erg, Gamma0, 0.1, 0.3, true);
    // tests(32, 32, 32, n_ism, eps_e, eps_B, p, 1e52 * con::erg, Gamma0, 0.1, 0.3, true);

    // return 0;
    size_t resolu[] = {32};

    for (auto r : resolu) {
        r_num = theta_num = phi_num = r;
        std::ofstream file("benchmark" + std::to_string(phi_num) + "-" + std::to_string(theta_num) + "-" +
                           std::to_string(r_num) + ".txt");

        for (size_t i = 0; i < E_iso.size(); ++i) {
            for (size_t j = 0; j < theta_c.size(); ++j) {
                auto start = std::chrono::high_resolution_clock::now();
                for (size_t k = 0; k < theta_v.size(); ++k) {
                    tests(r_num, theta_num, phi_num, n_ism, eps_e, eps_B, p, E_iso[i], Gamma0, theta_c[j], theta_v[k]);
                    // tests(r_num, theta_num, phi_num, n_ism, eps_e, eps_B, p, E_iso[i], Gamma0, theta_c[j], 0);
                }
                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                file << duration.count() / 1000000. / theta_v.size() << std::endl;
            }
        }
    }
    return 0;
}
