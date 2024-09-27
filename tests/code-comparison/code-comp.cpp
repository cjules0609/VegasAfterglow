#include <boost/numeric/odeint.hpp>
#include <filesystem>
#include <fstream>

#include "../../afterglow.h"
#include "json.hpp"

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

    double sigma = 0;

    auto jet_top = TophatJet(theta_c, E_iso, Gamma0);

    auto jet_gauss = GaussianJet(theta_c, E_iso, Gamma0, 1);

    Jet* jet = &jet_top;

    if (jet_type == "Gaussian") {
        jet = &jet_gauss;
    } else if (jet_type == "tophat") {
    } else {
        throw std::runtime_error("Jet type not recognized");
    }

    size_t r_num = 800;
    size_t theta_num = 150;
    size_t phi_num = 100;

    double R_dec = dec_radius(E_iso, n_ism, Gamma0, jet->duration);

    auto r = logspace(R_dec / 1000, R_dec * 500, r_num);
    auto theta = adaptive_theta_space(theta_num, jet->Gamma0_profile, theta_w);
    // auto theta = adaptive_theta_space(theta_num, jet.Gamma0_profile);
    auto phi = linspace(0, 2 * con::pi, phi_num);

    Coord coord{r, theta, phi};

    // solve dynamics
    Shock f_shock(coord, eps_e, eps_B, p);

    solve_shocks(coord, *jet, medium, f_shock);

    auto syn_e = gen_syn_electrons(coord, f_shock);

    auto syn_ph = gen_syn_photons(syn_e, coord, f_shock);

    IC_cooling_KN(syn_e, syn_ph, f_shock);

    Array t_bins = logspace(t_obs[0] * con::sec / 10, t_obs[1] * con::sec, 200);

    Observer obs;

    obs.observe(coord, f_shock, theta_view, lumi_dist, z);

    Array band_pass = logspace(eVtoHz(band_pass_[0] * con::keV), eVtoHz(band_pass_[1] * con::keV), 10);

    Array F_nu_syn = obs.flux(t_bins, band_pass, syn_ph);
    Array F_nu_syn_no_cool = obs.flux(t_bins, band_pass, syn_ph);

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
    lc_gen("/Users/yihanwang/Projects/afterglow-code-comparison/tests/case1");
    lc_gen("/Users/yihanwang/Projects/afterglow-code-comparison/tests/case2");
    lc_gen("/Users/yihanwang/Projects/afterglow-code-comparison/tests/case3");
    lc_gen("/Users/yihanwang/Projects/afterglow-code-comparison/tests/case4");
    lc_gen("/Users/yihanwang/Projects/afterglow-code-comparison/tests/case5");
    return 0;
}