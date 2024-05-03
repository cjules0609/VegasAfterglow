#include "afterglow.h"

#include <boost/numeric/odeint.hpp>

void GCN36236(std::string prefix, double Gamma0, double theta_j, Array const& theta_obs, double eps_e) {
    double E_iso = 1e53 * con::erg * (5 * con::deg * 5 * con::deg / theta_j / theta_j);

    double n_ism = pow(10, -1.7) / con::cm / con::cm / con::cm;
    // double n_ism = 1 / con::cm / con::cm / con::cm;
    double eps_B = pow(10, -3.7);
    double p = 2.139;
    //  double theta_j = 1.0 * con::deg;

    // create model
    auto medium = create_ISM(n_ism, eps_e, eps_B);
    // auto jet = create_tophat_jet(E_iso, Gamma0, theta_j, 2 * con::sec);
    auto jet = create_gaussian_jet(E_iso, Gamma0, theta_j, 1 * con::sec);

    size_t r_num = 800;
    size_t theta_num = 50;
    size_t phi_num = 37;
    double R_dec = dec_radius(E_iso, n_ism, Gamma0, jet.duration);

    auto r = logspace(R_dec / 1000, R_dec * 1000, r_num);
    auto theta = adaptive_theta_space(theta_num, jet.Gamma0_profile, 0.6);
    auto phi = linspace(0, 2 * con::pi, phi_num);

    // Coord coord{r_min, r_max, con::pi / 2, r_num, theta_num, phi_num};
    Coord coord{r, theta, phi};
    write2file(coord, prefix + "coord");

    // solve dynamics
    auto [shock_r, shock_f] = gen_shocks(coord, jet, medium);
    write2file(shock_f, prefix + "shock");

    auto syn_e = gen_syn_electrons_w_IC_cooling(p, coord, shock_f, medium);
    // auto syn_e = gen_syn_electrons(p, coord, shock_f);
    auto syn_ph = gen_syn_photons(syn_e, coord, shock_f);

    write2file(syn_ph, prefix + "syn_ph");
    write2file(syn_e, prefix + "syn_e");

    //

    Array t_bins = logspace(5e-2 * con::day, 5e3 * con::day, 100);
    // double lumi_dist = 6.6e26 * con::cm;
    double lumi_dist = 1.23e26 * con::cm;
    double z = luminosity_distance_to_z(lumi_dist);

    for (auto theta_v : theta_obs) {
        std::cout << theta_j / con::deg << ' ' << theta_v / con::deg << '\n';
        Observer obs;

        obs.observe(coord, shock_f, theta_v, lumi_dist, z);

        Array band_pass = logspace(eVtoHz(0.3 * con::keV), eVtoHz(10 * con::keV), 5);

        auto to_suffix = [](double theta) { return std::to_string(int(ceil(theta / con::deg))); };

        Array F_syn = obs.flux(t_bins, band_pass, syn_ph);
        std::string fname = prefix + "F_nu_syn_" + to_suffix(theta_v) + "_" + to_suffix(theta_j);
        write2file(F_syn, fname, con::erg / con::sec / con::cm / con::cm);
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
    double theta = 5 * con::deg, eps_e = 0.01;
    GCN36236("tests/", 300, theta, {0.54}, eps_e);
    /* Array theta_obs = linspace(5 * con::deg, 30 * con::deg, 26);

     double eps_e = 0.1;
     Array theta_c = {3 * con::deg, 4 * con::deg, 5 * con::deg, 6 * con::deg, 7 * con::deg, 8 * con::deg};
     std::vector<std::thread> threads;

     for (auto theta : theta_c) {
         threads.emplace_back(GCN36236, "prediction-n1/", 300, theta, theta_obs, eps_e);
     }

     for (auto& thread : threads) {
         thread.join();
     }*/

    // afterglow_gen();
    return 0;
}