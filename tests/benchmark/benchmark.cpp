#include <boost/numeric/odeint.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>

#include "afterglow.h"

void tests(size_t r_num, size_t theta_num, size_t phi_num, Real n_ism, Real eps_e, Real eps_B, Real p, Real E_iso,
           Real Gamma0, Real theta_c, Real theta_v, bool verbose = false) {
    Real z = 0.009;

    Real lumi_dist = 1.23e26 * unit::cm;

    Array t_obs = xt::logspace(std::log10(1e3 * unit::sec), std::log10(1e7 * unit::sec), 30);

    ISM medium(n_ism);

    GaussianJet jet(theta_c, E_iso, Gamma0);

    Coord coord = auto_grid(jet, t_obs, 0.6, theta_v, z, phi_num, theta_num, r_num);

    Shock f_shock = generate_fwd_shock(coord, medium, jet, eps_e, eps_B);

    Observer obs;

    // obs.observe_at(t_obs, coord, f_shock, lumi_dist, z);

    obs.observe(coord, f_shock, lumi_dist, z);

    auto syn_e = generate_syn_electrons(f_shock, p);

    auto syn_ph = generate_syn_photons(f_shock, syn_e);

    Real nu_obs = eVtoHz(1 * unit::keV);

    Array F_nu = obs.specific_flux(t_obs, nu_obs, syn_ph);

    if (verbose) {
        write_npz("F_nu" + std::to_string(phi_num) + "-" + std::to_string(theta_num) + "-" + std::to_string(r_num),
                  "F_nu", xt::eval(F_nu / unit::Jy), "t_obs", xt::eval(t_obs / unit::sec));
    }

    return;
}

int main() {
    Real n_ism = 2 / unit::cm3;
    Real eps_e = 1e-2;
    Real eps_B = 1e-4;
    Real p = 2.1;
    Real Gamma0 = 300;

    Array E_iso = xt::logspace(std::log10(1e48 * unit::erg), std::log10(1e52 * unit::erg), 10);
    Array theta_c = xt::linspace(0.01, 0.1, 100);
    Array theta_v = xt::linspace(0.01, 0.5, 10);

    size_t resolu[] = {256, 128, 64, 32, 30, 28, 25, 24, 16, 8};

    for (auto r : resolu) {
        tests(r, r, r, n_ism, eps_e, eps_B, p, 1e52 * unit::erg, Gamma0, 0.1, 0.3, true);
    }
    return 0;

    size_t benchmark_resolu[] = {16, 24, 32, 64};

    for (auto r : benchmark_resolu) {
        std::ofstream file("benchmark" + std::to_string(r) + "-" + std::to_string(r) + "-" + std::to_string(r) +
                           ".txt");

        for (size_t i = 0; i < E_iso.size(); ++i) {
            for (size_t j = 0; j < theta_c.size(); ++j) {
                auto start = std::chrono::high_resolution_clock::now();
                for (size_t k = 0; k < theta_v.size(); ++k) {
                    tests(r, r, r, n_ism, eps_e, eps_B, p, E_iso[i], Gamma0, theta_c[j], theta_v[k]);
                    // tests(r, r, r, n_ism, eps_e, eps_B, p, E_iso[i], Gamma0, theta_c[j], 0);
                }
                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                file << duration.count() / 1000000. / theta_v.size() << std::endl;
            }
        }
    }
    return 0;
}
