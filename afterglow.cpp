#include "afterglow.h"

#include <boost/numeric/odeint.hpp>
#include <filesystem>
#include <fstream>

#include "json.hpp"

void test() {
    double n_ism = 1 / con::cm3;
    double eps_e = 0.1;
    double eps_B = 0.001;
    double eps_e_r = 0.3;
    double eps_B_r = 0.1;
    double p = 2.2;
    double E_iso = 1e53 * con::erg;
    double Gamma0 = 300;
    double theta_c = 0.1;

    // create model
    auto medium = create_ISM(n_ism);
    auto jet = create_tophat_jet(E_iso, Gamma0, theta_c, 1 * con::sec);
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
    Shock r_shock(coord, eps_e_r, eps_B_r, p);

    solve_shocks(coord, jet, medium, f_shock, r_shock);

    write2file(r_shock, "tests/r_shock");
    write2file(f_shock, "tests/f_shock");

    auto syn_ph = gen_syn_photons_w_IC_cooling(coord, f_shock);

    auto syn_ph_r = gen_syn_photons_w_IC_cooling(coord, r_shock);

    Array t_bins = logspace(1e-5 * con::day, 1e2 * con::day, 100);

    double theta_v = 0.;
    double lumi_dist = 1.23e26 * con::cm;
    double z = 0.01;

    Observer obs;

    obs.observe(coord, f_shock, theta_v, lumi_dist, z);

    // Array band_pass = logspace(eVtoHz(0.3 * con::keV), eVtoHz(10 * con::keV), 5);

    // Array band_pass = logspace(eVtoHz(1e2 * con::keV), eVtoHz(1e3 * con::keV), 10);

    Array band_pass = {1e11 * con::Hz, 1e12 * con::Hz, 1e13 * con::Hz, 1e14 * con::Hz, 1e15 * con::Hz, 1e16 * con::Hz};

    auto to_suffix = [](double theta) { return std::to_string(int(ceil(theta / con::deg))); };

    auto F_syn = obs.specific_flux(t_bins, band_pass, syn_ph);

    write2file(syn_ph, "tests/ph_f");

    write2file(F_syn, "tests/Flux_f", con::erg / con::sec / con::cm / con::cm);

    auto F_syn_r = obs.specific_flux(t_bins, band_pass, syn_ph_r);

    write2file(syn_ph_r, "tests/ph_r");

    write2file(F_syn_r, "tests/Flux_r", con::erg / con::sec / con::cm / con::cm);

    write2file(boundary2centerlog(t_bins), "tests/t_obs", con::sec);
}

void GCN36236(std::string folder_name) {
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

    double n_ism = data["n_ism"];
    n_ism /= (con::cm * con::cm * con::cm);

    double eps_e = data["epsilon_e"];
    double eps_B = data["epsilon_B"];
    double p = data["p"];

    double theta_view = data["theta_view"];

    Array t_obs = data["t_obs"];

    Array band_pass_ = data["band pass (kev)"];

    // create model
    auto medium = create_ISM(n_ism);
    // auto jet = create_tophat_jet(E_iso, Gamma0, theta_j, 2 * con::sec);
    auto jet = create_gaussian_jet(E_iso, Gamma0, theta_c, 100 * con::sec);

    size_t r_num = 800;
    size_t theta_num = 150;
    size_t phi_num = 37;
    double R_thin = thin_shell_dec_radius(E_iso, n_ism, Gamma0);
    double R_dec = dec_radius(E_iso, n_ism, Gamma0, jet.duration);
    double R_cross = RS_crossing_radius(E_iso, n_ism, Gamma0, jet.duration);
    double R_spread = shell_spreading_radius(Gamma0, jet.duration);
    double R_N = RS_transition_radius(E_iso, n_ism, Gamma0, jet.duration);

    std::cout << R_thin / con::cm << ' ' << R_cross / con::cm << ' ' << R_spread / con::cm << ' ' << R_N / con::cm
              << '\n';

    auto r = logspace(R_dec / 1000, R_dec * 1000, r_num);
    auto theta = adaptive_theta_space(theta_num, jet.Gamma0_profile, 0.6);
    auto phi = linspace(0, 2 * con::pi, phi_num);

    // Coord coord{r_min, r_max, con::pi / 2, r_num, theta_num, phi_num};
    Coord coord{r, theta, phi};

    // solve dynamics
    Shock f_shock(coord, eps_e, eps_B, p);
    Shock r_shock(coord, eps_e, eps_B, p);

    solve_shocks(coord, jet, medium, f_shock, r_shock);

    auto syn_e = gen_syn_electrons_w_IC_cooling(coord, f_shock);
    // auto syn_e = gen_syn_electrons(coord, shock_f);
    auto syn_ph = gen_syn_photons(syn_e, coord, f_shock);

    Array t_bins = logspace(5e-4 * con::day, 5e3 * con::day, 100);

    Array theta_obs = {0, 5 * con::deg, 10 * con::deg, 15 * con::deg, 20 * con::deg, 25 * con::deg, 30 * con::deg};

    for (auto theta_v : theta_obs) {
        Observer obs;

        obs.observe(coord, f_shock, theta_v, lumi_dist, z);

        Array band_pass = logspace(eVtoHz(0.3 * con::keV), eVtoHz(10 * con::keV), 5);

        auto to_suffix = [](double theta) { return std::to_string(int(ceil(theta / con::deg))); };

        Array F_syn = obs.flux(t_bins, band_pass, syn_ph);

        std::string fname = "F_nu_syn_" + to_suffix(theta_v) + "_" + to_suffix(theta_v);

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

    write2file(boundary2centerlog(t_bins), "t_obs", con::sec);
}

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
}
int main() {
    // double theta = 5 * con::deg, eps_e = 0.01;
    // GCN36236("tests/", 300, theta, {0.54}, eps_e);
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
    // GCN36236("/Users/yihanwang/Projects/afterglow-code-comparison/tests/case6");
    //  afterglow_gen();
    /*lc_gen("/Users/yihanwang/Projects/afterglow-code-comparison/tests/case1");
    lc_gen("/Users/yihanwang/Projects/afterglow-code-comparison/tests/case2");
    lc_gen("/Users/yihanwang/Projects/afterglow-code-comparison/tests/case3");
    lc_gen("/Users/yihanwang/Projects/afterglow-code-comparison/tests/case4");
    lc_gen("/Users/yihanwang/Projects/afterglow-code-comparison/tests/case5");*/
    test();
    return 0;
}