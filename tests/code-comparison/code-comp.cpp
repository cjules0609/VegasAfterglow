#include <boost/numeric/odeint.hpp>
#include <filesystem>
#include <fstream>

#include "afterglow.h"
#include "json.hpp"

void lc_gen(std::string folder_name) {
    using json = nlohmann::json;

    std::ifstream f(folder_name + "/problem-setups.json");
    json data = json::parse(f);

    Real E_iso = data["E_iso"];
    E_iso *= con::erg;

    Real lumi_dist = data["luminosity distance"];
    lumi_dist *= con::cm;
    Real z = data["z"];
    std::string jet_type = data["jet type"];
    Real theta_c = data["theta_core"];
    Real theta_w = data["theta_wing"];
    Real Gamma0 = data["Gamma0"];

    bool ic_cool = data["inverse compton cooling"];

    Real n_ism = data["n_ism"];
    n_ism /= (con::cm * con::cm * con::cm);

    Real eps_e = data["epsilon_e"];
    Real eps_B = data["epsilon_B"];
    Real p = data["p"];

    Real theta_view = data["theta_view"];

    std::vector<Real> t_obs = data["t_obs"];

    std::vector<Real> band_pass_ = data["band pass (kev)"];

    Array t_bins = logspace(t_obs[0] * con::sec / 10, t_obs[1] * con::sec, 100);
    // create model
    auto medium = createISM(n_ism);

    Real sigma = 0;

    Ejecta jet;

    if (jet_type == "Gaussian") {
        jet.dEdOmega = math::gaussian(theta_c, E_iso / (4 * con::pi));
        jet.Gamma0 = math::gaussian(theta_c, Gamma0);
    } else if (jet_type == "tophat") {
        jet.dEdOmega = math::tophat(theta_c, E_iso / (4 * con::pi));
        jet.Gamma0 = math::tophat(theta_c, Gamma0);
    } else {
        throw std::runtime_error("Jet type not recognized");
    }

    size_t r_num = 128;
    size_t theta_num = 128;
    size_t phi_num = 128;

    // Coord coord = adaptiveGrid(medium, jet, inject::none, t_bins, theta_w, phi_num, theta_num, r_num);
    Coord coord = adaptiveGrid(medium, jet, inject::none, t_bins, theta_w);

    // solve dynamics
    Shock f_shock = genForwardShock(coord, medium, jet, inject::none, eps_e, eps_B);

    auto syn_e = genSynElectrons(f_shock, p);

    auto syn_ph = genSynPhotons(f_shock, syn_e);

    Observer obs(coord, f_shock, theta_view, lumi_dist, z);

    Array band_pass = logspace(eVtoHz(band_pass_[0] * con::keV), eVtoHz(band_pass_[1] * con::keV), 50);

    namespace fs = std::filesystem;

    std::string working_dir = folder_name + "/" + "VegasAfterglow";

    fs::create_directory(working_dir);

    std::ofstream file(working_dir + "/flux.csv");

    if (ic_cool) {
        Array F_nu_syn = obs.flux(t_bins, band_pass, syn_ph);
        for (size_t i = 0; i < t_bins.size(); ++i) {
            file << t_bins[i] / con::sec << ',' << F_nu_syn[i] / (con::erg / con::cm / con::cm / con::sec) << '\n';
        }
    } else {
        Array F_nu_syn_no_cool = obs.flux(t_bins, band_pass, syn_ph);
        for (size_t i = 0; i < t_bins.size(); ++i) {
            file << t_bins[i] / con::sec << ',' << F_nu_syn_no_cool[i] / (con::erg / con::cm / con::cm / con::sec)
                 << '\n';
        }
    }
    std::cout << "finish" + working_dir << '\n';
    // specify observables
}

int main() {
    std::vector<std::thread> threads;

    threads.emplace_back(std::thread(lc_gen, "/Users/yihanwang/Projects/afterglow-code-comparison/tests/case1"));
    threads.emplace_back(std::thread(lc_gen, "/Users/yihanwang/Projects/afterglow-code-comparison/tests/case2"));
    threads.emplace_back(std::thread(lc_gen, "/Users/yihanwang/Projects/afterglow-code-comparison/tests/case3"));
    threads.emplace_back(std::thread(lc_gen, "/Users/yihanwang/Projects/afterglow-code-comparison/tests/case4"));
    threads.emplace_back(std::thread(lc_gen, "/Users/yihanwang/Projects/afterglow-code-comparison/tests/case5"));

    for (auto& t : threads) {
        t.join();
    }
    /*lc_gen("/Users/yihanwang/Projects/afterglow-code-comparison/tests/case1");
    lc_gen("/Users/yihanwang/Projects/afterglow-code-comparison/tests/case2");
    lc_gen("/Users/yihanwang/Projects/afterglow-code-comparison/tests/case3");
    lc_gen("/Users/yihanwang/Projects/afterglow-code-comparison/tests/case4");
    lc_gen("/Users/yihanwang/Projects/afterglow-code-comparison/tests/case5");*/
    return 0;
}