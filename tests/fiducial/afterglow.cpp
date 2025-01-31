#include "../../afterglow.h"

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

    double n_ism = 1e-3 / con::cm3;
    double eps_e = 0.15;
    double eps_B = 0.02;
    double p = 2.1;
    double E_iso = 1e53 * con::erg;
    double Gamma0 = 300;
    double theta_c = 0.1;

    // create model
    auto medium = createISM(n_ism);
    auto jet = TophatJet(theta_c, E_iso, Gamma0);
    // auto jet = create_gaussian_jet(E_iso, Gamma0, theta_c, 1 * con::sec);

    size_t r_num = 500;
    size_t theta_num = 250;
    size_t phi_num = 37;

    double R_thin = thinShellDecRadius(E_iso, n_ism, Gamma0);
    double R_dec = decRadius(E_iso, n_ism, Gamma0, jet.duration);
    double R_cross = RShockCrossingRadius(E_iso, n_ism, Gamma0, jet.duration);
    double R_spread = shellSpreadingRadius(Gamma0, jet.duration);
    double R_N = RSTransitionRadius(E_iso, n_ism, Gamma0, jet.duration);

    std::cout << R_thin / con::cm << ' ' << R_cross / con::cm << ' ' << R_spread / con::cm << ' ' << R_N / con::cm
              << '\n';

    auto r = logspace(R_dec / 100, R_dec * 100, r_num);
    auto theta = adaptiveThetaSpace(theta_num, jet.Gamma0_profile);
    auto phi = linspace(0, 2 * con::pi, phi_num);

    // Coord coord{r_min, r_max, con::pi / 2, r_num, theta_num, phi_num};
    Coord coord{r, theta, phi};

    output(coord, data_folder + "/coord");

    // solve dynamics
    // auto [r_shock, f_shock] = gen_shocks(coord, jet, medium);
    Shock f_shock = genForwardShock(coord, jet, medium, eps_e, eps_B, p);
    // Shock r_shock(coord, eps_e_r, eps_B_r, p);

    output(f_shock, data_folder + "/f_shock");

    auto syn_e = genSynElectrons(coord, f_shock);

    auto syn_ph = genSynPhotons(syn_e, coord, f_shock);

    output(syn_e, data_folder + "/syn_e");

    output(syn_ph, data_folder + "/syn_ph");

    size_t j = 0;
    size_t k = 200;
    size_t resol = 60;
    double gamma_min = 6e3;
    double gamma_max = 4e6;
    double nu_min = 1e12 * con::Hz;
    double nu_max = 1e20 * con::Hz;

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