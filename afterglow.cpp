#include "afterglow.h"

#include <boost/numeric/odeint.hpp>
#include <filesystem>
#include <fstream>

#include "json.hpp"

void tests() {
    double n_ism = 1e-3 / con::cm3;
    double eps_e = 0.15;
    double eps_B = 0.02;
    double p = 2.1;
    double E_iso = 1e54 * con::erg;
    double Gamma0 = 300;
    double theta_c = 0.1;

    // create model
    auto medium = create_ISM(n_ism);
    auto jet = create_tophat_jet(E_iso, Gamma0, 0, theta_c, 1 * con::sec);
    // auto jet = create_gaussian_jet(E_iso, Gamma0, theta_c, 1 * con::sec);

    size_t r_num = 500;
    size_t theta_num = 250;
    size_t phi_num = 37;

    double R_thin = thin_shell_dec_radius(E_iso, n_ism, Gamma0);
    double R_dec = dec_radius(E_iso, n_ism, Gamma0, jet.duration);
    double R_cross = RS_crossing_radius(E_iso, n_ism, Gamma0, jet.duration);
    double R_spread = shell_spreading_radius(Gamma0, jet.duration);
    double R_N = RS_transition_radius(E_iso, n_ism, Gamma0, jet.duration);

    std::cout << R_thin / con::cm << ' ' << R_cross / con::cm << ' ' << R_spread / con::cm << ' ' << R_N / con::cm
              << '\n';

    auto r = logspace(R_dec / 100, R_dec * 100, r_num);
    auto theta = adaptive_theta_space(theta_num, jet.Gamma0_profile);
    auto phi = linspace(0, 2 * con::pi, phi_num);

    // Coord coord{r_min, r_max, con::pi / 2, r_num, theta_num, phi_num};
    Coord coord{r, theta, phi};

    write2file(coord, "tests/coord");

    // solve dynamics
    // auto [r_shock, f_shock] = gen_shocks(coord, jet, medium);
    Shock f_shock(coord, eps_e, eps_B, p);
    // Shock r_shock(coord, eps_e_r, eps_B_r, p);

    solve_shocks(coord, jet, medium, f_shock);

    write2file(f_shock, "tests/f_shock");

    auto syn_e = gen_syn_electrons(coord, f_shock);

    auto syn_ph = gen_syn_photons(syn_e, coord, f_shock);

    write2file(syn_e, "tests/syn_e");

    write2file(syn_ph, "tests/syn_ph");

    size_t j = 0;
    size_t k = 200;
    size_t resol = 60;
    double gamma_min = 6e3;
    double gamma_max = 4e6;
    double nu_min = 1e12 * con::Hz;
    double nu_max = 1e20 * con::Hz;

    auto syn_e_sptr = co_moving_e_spectrum(resol, gamma_min, gamma_max, syn_e[j][k]);

    write2file(syn_e_sptr, "tests/syn_e_sptr");

    auto syn_ph_sptr = co_moving_spectrum(resol, nu_min, nu_max, syn_ph[j][k]);

    write2file(syn_ph_sptr, "tests/syn_ph_sptr");

    IC_cooling_Thomson(syn_e, syn_ph, f_shock);

    auto ICT_ph = gen_syn_photons(syn_e, coord, f_shock);

    write2file(syn_e, "tests/ICT_e");

    write2file(ICT_ph, "tests/ICT_ph");

    auto ict_e_sptr = co_moving_e_spectrum(resol, gamma_min, gamma_max, syn_e[j][k]);

    write2file(ict_e_sptr, "tests/ict_e_sptr");

    auto ict_ph_sptr = co_moving_spectrum(resol, nu_min, nu_max, ICT_ph[j][k]);

    write2file(ict_ph_sptr, "tests/ict_ph_sptr");

    IC_cooling_KN(syn_e, syn_ph, f_shock);

    auto ICKN_ph = gen_syn_photons(syn_e, coord, f_shock);

    auto& ee = syn_e[j][k];

    std::cout << ee.gamma_c << " " << ee.gamma_m << " " << ee.Ys[0].gamma_hat_m << " " << ee.Ys[0].gamma_hat_c << " "
              << ee.Ys[0].Y_T << '\n';

    std::cout << ee.Ys[0].nu_hat_m << " " << ee.Ys[0].nu_hat_c << '\n';

    write2file(syn_e, "tests/ICKN_e");

    write2file(ICKN_ph, "tests/ICKN_ph");

    auto ickn_e_sptr = co_moving_e_spectrum(resol, gamma_min, gamma_max, syn_e[j][k]);

    write2file(ickn_e_sptr, "tests/ickn_e_sptr");

    auto ickn_ph_sptr = co_moving_spectrum(resol, nu_min, nu_max, ICKN_ph[j][k]);

    write2file(ickn_ph_sptr, "tests/ickn_ph_sptr");
}

/*
void lc_gen(std::string folder_name) {
    using json = nlohmann::json;

    std::ifstream f(folder_name + "/problem-setups.json");
    json data = json::parse(f);

    double E_iso = data["E_iso"];
    E_iso *= con::erg;

    double lumi_dist = data["luminosity distance"];
    lumi_dist *= con::cm;
    double z = data["z"];
    std::string jet_type = data["jet type"];
    double theta_c = data["theta_core"];
    double theta_w = data["theta_wing"];
    double Gamma0 = data["Gamma0"];

    bool ic_cool = data["inverse compton cooling"];

    double n_ism = data["n_ism"];
    n_ism /= (con::cm * con::cm * con::cm);

    double eps_e = data["epsilon_e"];
    double eps_B = data["epsilon_B"];
    double p = data["p"];

    double theta_view = data["theta_view"];

    Array t_obs = data["t_obs"];

    Array band_pass_ = data["band pass (kev)"];

    std::cout << E_iso / con::erg << ' ' << lumi_dist / con::cm << ' ' << z << ' ' << jet_type << ' ' << theta_c << ' '
              << theta_w << ' ' << Gamma0 << ' ' << n_ism * con::cm * con::cm * con::cm << ' ' << eps_e << ' ' << eps_B
              << ' ' << p << ' ' << theta_view << ' ' << t_obs[0] << ' ' << t_obs[1] << ' ' << band_pass_[0] << ' '
              << band_pass_[1] << '\n';

    // create model
    auto medium = create_ISM(n_ism);

    auto jet = create_tophat_jet(E_iso, Gamma0, theta_c, 1 * con::sec);
    if (jet_type == "Gaussian") {
        jet = create_gaussian_jet(E_iso, Gamma0, theta_c, 1 * con::sec);
    } else if (jet_type == "tophat") {
    } else {
        throw std::runtime_error("Jet type not recognized");
    }

    size_t r_num = 1000;
    size_t theta_num = 250;
    size_t phi_num = 100;

    double R_dec = dec_radius(E_iso, n_ism, Gamma0, jet.duration);

    auto r = logspace(R_dec / 1000, R_dec * 500, r_num);
    auto theta = adaptive_theta_space(theta_num, jet.Gamma0_profile, theta_w);
    // auto theta = adaptive_theta_space(theta_num, jet.Gamma0_profile);
    auto phi = linspace(0, 2 * con::pi, phi_num);

    Coord coord{r, theta, phi};

    // solve dynamics
    Shock f_shock(coord, eps_e, eps_B, p);
    Shock r_shock(coord, eps_e, eps_B, p);

    solve_shocks(coord, jet, medium, f_shock, r_shock);

    auto syn_e = gen_syn_electrons_w_IC_cooling(coord, f_shock);
    auto syn_ph = gen_syn_photons(syn_e, coord, f_shock);

    auto syn_e_no_cooling = gen_syn_electrons(coord, f_shock);
    auto syn_ph_no_cooling = gen_syn_photons(syn_e_no_cooling, coord, f_shock);

    Array t_bins = logspace(t_obs[0] * con::sec / 10, t_obs[1] * con::sec, 200);

    Observer obs;

    obs.observe(coord, f_shock, theta_view, lumi_dist, z);

    Array band_pass = logspace(eVtoHz(band_pass_[0] * con::keV), eVtoHz(band_pass_[1] * con::keV), 10);

    Array F_nu_syn = obs.flux(t_bins, band_pass, syn_ph);
    Array F_nu_syn_no_cool = obs.flux(t_bins, band_pass, syn_ph_no_cooling);

    auto t_c = boundary2centerlog(t_bins);

    namespace fs = std::filesystem;

    std::string working_dir = folder_name + "/" + "VegasAfterglow";

    fs::create_directory(working_dir);

    std::ofstream file(working_dir + "/flux.csv");

    if (ic_cool) {
        for (size_t i = 0; i < t_c.size(); ++i) {
            file << t_c[i] / con::sec << ',' << F_nu_syn[i] / (con::erg / con::cm / con::cm / con::sec) << '\n';
        }
    } else {
        for (size_t i = 0; i < t_c.size(); ++i) {
            file << t_c[i] / con::sec << ',' << F_nu_syn_no_cool[i] / (con::erg / con::cm / con::cm / con::sec) << '\n';
        }
    }
    std::cout << "finish" + working_dir << '\n';
    // specify observables
}*/
int main() {
    tests();
    return 0;
}