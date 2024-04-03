#include "afterglow.h"

#include <boost/numeric/odeint.hpp>

void afterglow_gen() {
    std::string prefix = "tests/";
    double E_iso = 1e53 * con::erg;
    double Gamma0 = 300;
    double n_ism = 1 / con::cm / con::cm / con::cm;
    double eps_rad = 1;
    double eps_e = 0.1;
    double eps_B = 0.01;
    double pel = 2.3;
    double theta_j = 35.0 * con::deg;

    // create model
    auto medium = create_ISM(n_ism, eps_e, eps_B, pel);
    auto inj = create_iso_power_law_injection(0 / 9e20 / 2e33, 1000, 1, 2);
    auto jet = create_tophat_jet(E_iso, Gamma0, theta_j, inj);
    //  auto blast = create_power_law_jet(E_iso, Gamma0, theta_j, 4, inj);
    // auto jet = create_gaussian_jet(E_iso, Gamma0, theta_j / 6, inj);

    // generate grid
    double M0 = E_iso / (Gamma0 * con::c * con::c);
    double R_ES = pow(3 * M0 / (4 * con::pi * n_ism * con::mp * Gamma0), 1.0 / 3);
    size_t grid_num = 300;
    Array r = logspace(R_ES / 100, R_ES * 100, grid_num);
    Array theta = linspace(0, theta_j, 30);
    Array phi = linspace(0, 2 * con::pi, 2);
    Coord coord{r, theta, phi};

    write2file(coord.r, prefix + "r");
    write2file(coord.theta, prefix + "theta");
    write2file(coord.phi, prefix + "phi");

    // solve dynamics
    auto shock_f = gen_forward_shock(coord, jet, medium);
    write2file(shock_f.Gamma, prefix + "Gamma");

    auto syn_e = gen_syn_electrons(coord, shock_f, medium);
    auto IC_ph = gen_IC_photons(coord, shock_f, syn_e, syn_e, medium);

    auto Y_eff = create_grid_like(shock_f.Gamma);

    write2file(syn_e, prefix + "syn");
    write2file(Y_eff, prefix + "Y");

    auto syn0 = syn_e;
    Y_eff = IC_cooling_KN(shock_f, syn_e, IC_ph, medium);
    write2file(syn_e, prefix + "IC_syn");
    write2file(Y_eff, prefix + "IC_Y");

    syn_e = syn0;
    Y_eff = IC_cooling_Thomson(shock_f, syn_e, IC_ph, medium);
    write2file(syn_e, prefix + "ICnoKN_syn");
    write2file(Y_eff, prefix + "ICnoKN_Y");

    Observer obs;

    double theta_obs = 0.0 * con::deg;
    obs.observe(coord, shock_f, theta_obs);
    write2file(obs.t_obs, prefix + "t_obs");
    write2file(obs.doppler, prefix + "doppler");

    // specify observables
    double nu_obs_R = 1e9 * 500;
    double nu_obs_O = 1e14 * 500;
    double nu_obs_X = 1e17 * 500;
    Array nu_obs{nu_obs_R, nu_obs_O, nu_obs_X};

    size_t data_points = 100;

    for (size_t i = 0; i < nu_obs.size(); ++i) {
        MeshGrid3d I_syn_obs = obs.gen_I_nu_grid(nu_obs[i], syn_e);
        MeshGrid L_nu = obs.gen_light_curve(data_points, nu_obs[i], syn_e);

        write2file(I_syn_obs, prefix + "I_nu_" + std::to_string(int(log10(nu_obs[i] / 500))));
        write2file(L_nu, prefix + "L_nu_" + std::to_string(int(log10(nu_obs[i] / 500))));
    }

    // write to file
}

int main() {
    afterglow_gen();
    return 0;
}