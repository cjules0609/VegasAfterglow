#include "afterglow.h"

#include <boost/numeric/odeint.hpp>

void afterglow_gen() {
    std::string prefix = "tests/";
    double E_iso = 1e51 * con::erg;
    double Gamma0 = 100;
    double n_ism = 1 / con::cm / con::cm / con::cm;
    double eps_e = 0.1;
    double eps_B = 0.001;
    double p = 2.3;
    double theta_j = 10.0 * con::deg;

    // create model
    auto medium = create_ISM(n_ism, eps_e, eps_B);
    auto inj = create_iso_power_law_injection(0 / 9e20 / 2e33, 1000, 1, 2);
    auto jet = create_tophat_jet(E_iso, Gamma0, theta_j, inj);
    //  auto blast = create_power_law_jet(E_iso, Gamma0, theta_j, 4, inj);
    // auto jet = create_gaussian_jet(E_iso, Gamma0, theta_j / 6, inj);

    // generate grid
    double M0 = E_iso / (Gamma0 * con::c * con::c);
    double R_ES = pow(3 * M0 / (4 * con::pi * n_ism * con::mp * Gamma0), 1.0 / 3);
    size_t grid_num = 600;
    Array r = logspace(R_ES / 100, R_ES * 50, grid_num);
    Array theta = linspace(0, theta_j, 50);
    Array phi = linspace(0, 2 * con::pi, 37);
    Coord coord{r, theta, phi};

    write2file(coord.r, prefix + "r");
    write2file(coord.theta, prefix + "theta");
    write2file(coord.phi, prefix + "phi");

    // solve dynamics
    auto shock_f = gen_forward_shock(coord, jet, medium);
    write2file(shock_f.Gamma, prefix + "Gamma");
    write2file(shock_f.B, prefix + "B");
    write2file(shock_f.D_com, prefix + "D_com");
    write2file(shock_f.t_com, prefix + "t_com");

    auto syn_e = gen_syn_electrons(p, coord, shock_f, medium);
    auto syn_ph = gen_syn_photons(coord, syn_e, shock_f, medium);
    auto Y_eff = create_grid_like(shock_f.Gamma, 0);

    write2file(syn_ph, prefix + "syn");
    write2file(Y_eff, prefix + "Y");

    Y_eff = solve_IC_Y_Thomson(shock_f, syn_e, medium);
    auto syn_e_IC = gen_syn_electrons(p, coord, shock_f, medium, Y_eff);
    auto syn_ph_IC = gen_syn_photons(coord, syn_e_IC, shock_f, medium);
    write2file(syn_ph_IC, prefix + "syn_IC");
    write2file(Y_eff, prefix + "Y_IC");

    Y_eff = solve_IC_Y_KN(shock_f, syn_e, medium);
    auto syn_e_IC_KN = gen_syn_electrons(p, coord, shock_f, medium, Y_eff);
    auto syn_ph_IC_KN = gen_syn_photons(coord, syn_e_IC_KN, shock_f, medium);
    write2file(syn_ph_IC_KN, prefix + "syn_ICKN");
    write2file(Y_eff, prefix + "Y_ICKN");

    /*auto IC_ph = gen_IC_photons(coord, shock_f, syn_e, syn_ph, medium);

    auto IC_ph_IC = gen_IC_photons(coord, shock_f, syn_e_IC, syn_ph_IC, medium);

    auto IC_ph_ICKN = gen_IC_photons(coord, shock_f, syn_e_IC_KN, syn_ph_IC_KN, medium);*/

    Observer obs;

    double theta_obs = 0 * con::deg;
    obs.observe(coord, shock_f, theta_obs);
    write2file(obs.t_obs, prefix + "t_obs");
    write2file(obs.doppler, prefix + "doppler");

    // specify observables

    Array nu_obs{1e9 * con::Hz,  1e10 * con::Hz, 1e11 * con::Hz, 1e12 * con::Hz, 1e13 * con::Hz,
                 1e14 * con::Hz, 1e15 * con::Hz, 1e16 * con::Hz, 1e17 * con::Hz};

    size_t time_resol = 100;

    MeshGrid L_nu = obs.gen_light_curve(time_resol, nu_obs, syn_ph);
    write2file(L_nu, prefix + "L_nu");

    MeshGrid L_nu_IC = obs.gen_light_curve(time_resol, nu_obs, syn_ph_IC);
    write2file(L_nu_IC, prefix + "L_nu_IC");

    /*for (size_t i = 0; i < nu_obs.size(); ++i) {
        MeshGrid3d I_syn_obs = obs.gen_I_nu_grid(nu_obs[i], syn_ph);
        write2file(I_syn_obs, prefix + "I_nu_" + std::to_string(int(log10(nu_obs[i] / 500))));
    }

    for (size_t i = 0; i < nu_obs.size(); ++i) {
        MeshGrid3d I_syn_obs = obs.gen_I_nu_grid(nu_obs[i], syn_ph_IC);
        write2file(I_syn_obs, prefix + "I_nu_IC_" + std::to_string(int(log10(nu_obs[i] / 500))));
    }*/

    /*for (size_t i = 0; i < nu_obs.size(); ++i) {
        MeshGrid3d I_syn_obs = obs.gen_I_nu_grid(nu_obs[i], IC_ph);
        MeshGrid L_nu = obs.gen_light_curve(time_resol, nu_obs[i], IC_ph);

        write2file(I_syn_obs, prefix + "I_nu_SSC_" + std::to_string(int(log10(nu_obs[i] / 500))));
        write2file(L_nu, prefix + "L_nu_SSC_" + std::to_string(int(log10(nu_obs[i] / 500))));
    }*/

    /*auto syn_spectrum = full_spectrum(syn_ph[0], 1e4 * con::Hz, 1e25 * con::Hz);
    write2file(syn_spectrum, prefix + "syn_spectrum");

    auto syn_IC_spectrum = full_spectrum(syn_ph_IC[0], 1e4 * con::Hz, 1e25 * con::Hz);
    write2file(syn_IC_spectrum, prefix + "syn_IC_spectrum");

    auto syn_ICKN_spectrum = full_spectrum(syn_ph_IC_KN[0], 1e4 * con::Hz, 1e25 * con::Hz);
    write2file(syn_ICKN_spectrum, prefix + "syn_ICKN_spectrum");*/

    /*auto IC_spectrum = full_spectrum(IC_ph[0], 1e4 * con::Hz, 1e35 * con::Hz);
    write2file(IC_spectrum, prefix + "IC_spectrum");

    auto IC_IC_spectrum = full_spectrum(IC_ph_IC[0], 1e4 * con::Hz, 1e35 * con::Hz);
    write2file(IC_IC_spectrum, prefix + "IC_IC_spectrum");

    auto IC_ICKN_spectrum = full_spectrum(IC_ph_ICKN[0], 1e4 * con::Hz, 1e35 * con::Hz);
    write2file(IC_ICKN_spectrum, prefix + "IC_ICKN_spectrum");*/

    // write to file

    // write2file(IC_ph.nu_IC_, prefix + "IC_spectrum_test_t");
    // write2file(IC_ph.I_nu_, prefix + "IC_spectrum_test_I");
    /*size_t is[4] = {0, 10, 100, 200};

     for (auto i : is) {
         auto IC_ph = ICPhoton(syn_e_IC[0][i], syn_ph_IC[0][i]);
         auto IC_specturm = spectrum(IC_ph, 1e4 * con::Hz, 1e35 * con::Hz, 100);

         write2file(IC_specturm, prefix + "IC_spectrum_test" + std::to_string(i));
     }*/
}

int main() {
    afterglow_gen();
    return 0;
}