#include "afterglow.h"

#include <boost/numeric/odeint.hpp>

void GCN36236(std::string prefix, double Gamma0, double theta_j, Array const& theta_obs, double eps_e) {
    double E_iso = 1e53 * con::erg * (5 * con::deg * 5 * con::deg / theta_j / theta_j);
    // std::string prefix = "GCN36236/";
    // double E_iso = 1e53 * con::erg;
    // double Gamma0 = 300;
    // std::cout << eVtoHz(0.3 * con::keV) / con::Hz << eVtoHz(10 * con::keV) / con::Hz << '\n';
    double n_ism = pow(10, -1.7) / con::cm / con::cm / con::cm;
    // double n_ism = 1 / con::cm / con::cm / con::cm;
    double eps_B = pow(10, -3.7);
    double p = 2.139;
    //  double theta_j = 1.0 * con::deg;

    // create model
    auto medium = create_ISM(n_ism, eps_e, eps_B);
    // auto jet = create_tophat_jet(E_iso, Gamma0, theta_j, 2 * con::sec);

    auto jet = create_gaussian_jet(E_iso, Gamma0, theta_j, 2 * con::sec);

    // generate grid
    double M0 = E_iso / (Gamma0 * con::c * con::c);
    double R_ES = pow(3 * M0 / (4 * con::pi * n_ism * con::mp * Gamma0), 1.0 / 3);
    size_t r_num = 1500;
    size_t theta_num = 150;
    size_t phi_num = 37;

    double r_min = R_ES / 100;
    double r_max = R_ES * 100;

    auto r = logspace(r_min, r_max, r_num);
    // auto theta = adaptive_theta_space(theta_num, jet.Gamma0);
    auto theta = linspace(0, 0.6, theta_num);
    auto phi = linspace(0, 2 * con::pi, phi_num);

    // Coord coord{r_min, r_max, con::pi / 2, r_num, theta_num, phi_num};
    Coord coord{r, theta, phi};
    write2file(coord, prefix + "coord");

    // solve dynamics
    auto [shock_r, shock_f] = gen_shocks(coord, jet, medium);
    write2file(shock_f, prefix + "shock");

    auto syn_e = gen_syn_electrons_w_IC_cooling(p, coord, shock_f, medium);
    auto syn_ph = gen_syn_photons(syn_e, coord, shock_f);

    write2file(syn_ph, prefix + "syn_ph");
    write2file(syn_e, prefix + "syn_e");

    // Array theta_obs = {15 * con::deg, 20 * con::deg, 22 * con::deg, 24 * con::deg,
    //                   26 * con::deg, 28 * con::deg, 30 * con::deg, 35 * con::deg};

    Array t_bins = logspace(5e-2 * con::day, 5e2 * con::day, 100);

    for (auto theta_v : theta_obs) {
        std::cout << theta_j / con::deg << ' ' << theta_v / con::deg << '\n';
        Observer obs;
        double z = 0.046;
        obs.observe(coord, shock_f, theta_v, z);
        // std::cout << obs.lumi_dist / con::cm << '\n';
        obs.lumi_dist = 6.6e26 * con::cm;
        Array band_pass = logspace(eVtoHz(0.3 * con::keV), eVtoHz(10 * con::keV), 5);

        Array F_nu_syn = obs.gen_flux(t_bins, band_pass, syn_ph);
        write2file(F_nu_syn,
                   prefix + "Vegas_" + std::to_string(int(ceil(theta_v / con::deg))) + "_" +
                       std::to_string(int(ceil(theta_j / con::deg))),
                   con::erg / con::sec / con::cm / con::cm);
    }

    /* Observer obs;
     double z = 0.01;
     double theta_v = 0.54;
     obs.observe(coord, shock_f, theta_v, z);
     std::cout << obs.lumi_dist / con::cm << '\n';
     obs.lumi_dist = 1.23e26 * con::cm;
     Array band_pass = logspace(eVtoHz(0.3 * con::keV), eVtoHz(10 * con::keV), 5);

     Array F_nu_syn = obs.gen_flux(t_bins, band_pass, syn_ph);
     write2file(F_nu_syn, prefix + "F_nu_syn", con::erg / con::sec / con::cm / con::cm);*/

    // specify observables

    write2file(boundary2centerlog(t_bins), prefix + "t_obs", con::sec);
}
int main() {
    Array theta_obs;
    for (size_t i = 5; i < 30; i++) {
        theta_obs.push_back(i * con::deg);
    }
    double eps_e = 0.1;
    GCN36236("prediction-n1/", 300, 8 * con::deg, theta_obs, eps_e);

    GCN36236("prediction-n1/", 300, 7 * con::deg, theta_obs, eps_e);

    GCN36236("prediction-n1/", 300, 6 * con::deg, theta_obs, eps_e);

    GCN36236("prediction-n1/", 300, 5 * con::deg, theta_obs, eps_e);

    GCN36236("prediction-n1/", 300, 4 * con::deg, theta_obs, eps_e);

    GCN36236("prediction-n1/", 300, 3 * con::deg, theta_obs, eps_e);

    // GCN36236("GCN/", 1e53 * con::erg, 300, 0.088);
    //  GCN36236("code-comp/", 1e53 * con::erg, 300, 0.088);

    //  GCN36236("GCN36236-2/", pow(10, 52.4) * con::erg, 300, 8.0 * con::deg);
    // afterglow_gen();
    return 0;
}