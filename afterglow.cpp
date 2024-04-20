#include "afterglow.h"

#include <boost/numeric/odeint.hpp>

void afterglow_gen() {
    std::string prefix = "tests/";
    double E_iso = 1e53 * con::erg;
    double Gamma0 = 300;
    double n_ism = 0.1 / con::cm / con::cm / con::cm;
    double eps_e = 0.3;
    double eps_B = 0.001;
    double p = 2.3;
    double theta_j = 1.0 * con::deg;

    // create model
    auto medium = create_ISM(n_ism, eps_e, eps_B);
    auto inj = create_iso_power_law_injection(0 * con::erg / con::sec, 1000 * con::sec, 1 * con::sec, 2);
    auto jet = create_tophat_jet(E_iso, Gamma0, theta_j, inj);
    // auto blast = create_power_law_jet(E_iso, Gamma0, theta_j, 4, inj);
    // auto jet = create_gaussian_jet(E_iso, Gamma0, theta_j / 6, inj);

    // generate grid
    double M0 = E_iso / (Gamma0 * con::c * con::c);
    double R_ES = pow(3 * M0 / (4 * con::pi * n_ism * con::mp * Gamma0), 1.0 / 3);
    size_t r_num = 500;
    size_t theta_num = 30;
    size_t phi_num = 37;

    double r_min = R_ES / 100;
    double r_max = R_ES * 100;

    Coord coord{r_min, r_max, theta_j, r_num, theta_num, phi_num};

    write2file(coord, prefix + "coord");
    // solve dynamics
    auto [shock_r, shock_f] = gen_shocks(coord, jet, medium);
    write2file(shock_f, prefix + "shock");

    auto syn_e = gen_syn_electrons(p, coord, shock_f);
    auto syn_ph = gen_syn_photons(syn_e, coord, shock_f);
    write2file(syn_ph, prefix + "syn");

    auto Y_eff = create_grid_like(shock_f.Gamma, 0);
    write2file(Y_eff, prefix + "Y");

    Y_eff = solve_IC_Y_Thomson(syn_e, shock_f, medium);
    auto syn_e_IC = gen_syn_electrons(p, coord, shock_f, Y_eff);
    auto syn_ph_IC = gen_syn_photons(syn_e_IC, coord, shock_f);
    write2file(syn_ph_IC, prefix + "syn_IC");
    write2file(Y_eff, prefix + "Y_IC");

    Y_eff = solve_IC_Y_KN(syn_e, shock_f, medium);
    auto syn_e_IC_KN = gen_syn_electrons(p, coord, shock_f, Y_eff);
    auto syn_ph_IC_KN = gen_syn_photons(syn_e_IC_KN, coord, shock_f);
    write2file(syn_ph_IC_KN, prefix + "syn_ICKN");
    write2file(syn_e_IC_KN, prefix + "electron_ICKN");
    write2file(Y_eff, prefix + "Y_ICKN");

    auto IC_ph = gen_IC_photons(syn_e_IC_KN, syn_ph_IC_KN, shock_f);

    auto syn_spectrum = co_moving_spectrums(50, 1e4 * con::Hz, 1e25 * con::Hz, syn_ph_IC_KN[0]);
    write2file(syn_spectrum, prefix + "syn_spectrum");

    auto IC_spectrum = co_moving_spectrums(50, 1e4 * con::Hz, 1e43 * con::Hz, IC_ph[0]);
    write2file(IC_spectrum, prefix + "IC_spectrum");

    auto n_specturm = co_moving_n_spectrums(50, 1, 1e9, syn_e_IC_KN[0]);
    write2file(n_specturm, prefix + "n_spectrum");

    Observer obs;

    double theta_obs = 0 * con::deg;
    double z = 0.003;
    obs.observe(coord, shock_f, theta_obs, z);
    write2file(obs.t_obs, prefix + "t_obs");
    write2file(obs.doppler, prefix + "doppler");

    // specify observables
    Array nu_obs{1e14 * con::Hz, 1e15 * con::Hz, 1e16 * con::Hz, 1e17 * con::Hz, 1e18 * con::Hz,
                 1e19 * con::Hz, 1e20 * con::Hz, 1e21 * con::Hz, 1e22 * con::Hz, 1e23 * con::Hz,
                 1e24 * con::Hz, 1e25 * con::Hz, 1e26 * con::Hz, 1e27 * con::Hz, 1e28 * con::Hz};

    size_t time_resol = 50;

    MeshGrid F_nu_syn = obs.gen_light_curve(time_resol, nu_obs, syn_ph_IC_KN);
    write2file(F_nu_syn, prefix + "F_nu_syn");

    MeshGrid F_nu_IC = obs.gen_light_curve(time_resol, nu_obs, IC_ph);
    write2file(F_nu_IC, prefix + "F_nu_IC");

    MeshGrid F_nu_tot = obs.gen_light_curve(time_resol, nu_obs, IC_ph, syn_ph_IC_KN);
    write2file(F_nu_tot, prefix + "F_nu_tot");

    for (size_t i = 0; i < nu_obs.size(); ++i) {
        MeshGrid3d j_syn_obs = obs.gen_j_nu_grid(nu_obs[i], syn_ph_IC_KN);
        write2file(j_syn_obs, prefix + "j_nu_syn" + std::to_string(int(log10(nu_obs[i] / 500))));
    }

    for (size_t i = 0; i < nu_obs.size(); ++i) {
        MeshGrid3d j_IC_obs = obs.gen_j_nu_grid(nu_obs[i], IC_ph);
        write2file(j_IC_obs, prefix + "j_nu_IC" + std::to_string(int(log10(nu_obs[i] / 500))));
    }
}

int main() {
    afterglow_gen();
    return 0;
}