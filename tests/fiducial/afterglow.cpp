#include "afterglow.h"

#include <boost/numeric/odeint.hpp>
#include <filesystem>
#include <fstream>

#include "json.hpp"

void tests() {
    std::string data_folder = "data";
    if (std::filesystem::create_directory(data_folder)) {
    } else {
        std::cout << "Failed to create directory: " << data_folder << std::endl;
    }

    Real n_ism = 1e-3 / con::cm3;
    Real eps_e = 0.15;
    Real eps_B = 0.02;
    Real p = 2.1;
    Real E_iso = 1e53 * con::erg;
    Real Gamma0 = 300;
    Real theta_c = 0.1;

    // create model
    Medium medium;

    medium.rho = evn::ISM(n_ism);

    Ejecta jet;
    jet.dE0dOmega = math::tophat(theta_c, E_iso / (4 * con::pi));
    jet.Gamma0 = math::tophat(theta_c, Gamma0);

    size_t r_num = 500;
    size_t theta_num = 250;
    size_t phi_num = 37;

    Array t_obs = xt::logspace(std::log10(1e2 * con::sec), std::log10(1e8 * con::sec), 100);
    Coord coord = autoGrid(jet, t_obs, 0.5, 0, phi_num, theta_num, r_num);

    output(coord, data_folder + "/coord");

    // solve dynamics
    // auto [r_shock, f_shock] = gen_shocks(coord, jet, medium);
    Shock f_shock = genForwardShock(coord, medium, jet, eps_e, eps_B);
    // Shock r_shock(coord, eps_e_r, eps_B_r, p);

    output(f_shock, data_folder + "/f_shock");

    auto syn_e = genSynElectrons(f_shock, p);

    auto syn_ph = genSynPhotons(f_shock, syn_e);

    output(syn_ph, data_folder + "/syn_ph");

    size_t j = 0;
    size_t k = 200;
    size_t resol = 60;
    Real gamma_min = 6e3;
    Real gamma_max = 4e6;
    Real nu_min = 1e12 * con::Hz;
    Real nu_max = 1e20 * con::Hz;

    /*auto syn_e_sptr = co_moving_e_spectrum(resol, gamma_min, gamma_max, syn_e[j][k]);

    output(syn_e_sptr, data_folder + "/syn_e_sptr");

    auto syn_ph_sptr = co_moving_spectrum(resol, nu_min, nu_max, syn_ph[j][k]);

    output(syn_ph_sptr, data_folder + "/syn_ph_sptr");

    eCoolingThomson(syn_e, syn_ph, f_shock);

    auto ICT_ph = genSynPhotons(syn_e, coord, f_shock);

    output(syn_e, data_folder + "/ICT_e");

    output(ICT_ph, data_folder + "/ICT_ph");

    auto ict_e_sptr = co_moving_e_spectrum(resol, gamma_min, gamma_max, syn_e[j][k]);

    output(ict_e_sptr, data_folder + "/ict_e_sptr");

    auto ict_ph_sptr = co_moving_spectrum(resol, nu_min, nu_max, ICT_ph[j][k]);

    output(ict_ph_sptr, data_folder + "/ict_ph_sptr");

    eCoolingKleinNishina(syn_e, syn_ph, f_shock);

    auto ICKN_ph = genSynPhotons(syn_e, coord, f_shock);

    auto& ee = syn_e[j][k];

    std::cout << ee.gamma_c << " " << ee.gamma_m << " " << ee.Ys[0].gamma_hat_m << " " << ee.Ys[0].gamma_hat_c << " "
              << ee.Ys[0].Y_T << '\n';

    std::cout << ee.Ys[0].nu_hat_m << " " << ee.Ys[0].nu_hat_c << '\n';

    output(syn_e, data_folder + "/ICKN_e");

    output(ICKN_ph, data_folder + "/ICKN_ph");

    auto ickn_e_sptr = co_moving_e_spectrum(resol, gamma_min, gamma_max, syn_e[j][k]);

    output(ickn_e_sptr, data_folder + "/ickn_e_sptr");

    auto ickn_ph_sptr = co_moving_spectrum(resol, nu_min, nu_max, ICKN_ph[j][k]);

    output(ickn_ph_sptr, data_folder + "/ickn_ph_sptr");*/
}

int main() {
    tests();
    return 0;
}