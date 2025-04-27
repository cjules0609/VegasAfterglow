#include <boost/numeric/odeint.hpp>
#include <filesystem>
#include <fstream>

#include "afterglow.h"
#include "json.hpp"

void lc_gen(std::string folder_name, bool out = false) {
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

    Array t_bins = xt::logspace(std::log10(t_obs[0] * con::sec / 10), std::log10(t_obs[1] * con::sec), 100);
    // create model
    ISM medium(n_ism);

    Real sigma = 0;

    Ejecta jet;

    if (jet_type == "Gaussian") {
        jet.dE0dOmega = math::gaussian(theta_c, E_iso / (4 * con::pi));
        jet.Gamma0 = math::gaussian(theta_c, Gamma0);
    } else if (jet_type == "tophat") {
        jet.dE0dOmega = math::tophat(theta_c, E_iso / (4 * con::pi));
        jet.Gamma0 = math::tophat(theta_c, Gamma0);
    } else {
        throw std::runtime_error("Jet type not recognized");
    }
    jet.spreading = true;

    size_t t_num = 32;
    size_t theta_num = 32;
    size_t phi_num = 32;

    Coord coord = autoGrid(jet, t_bins, con::pi / 2, theta_view, z, phi_num, theta_num, t_num);

    // solve dynamics
    Shock f_shock = genForwardShock(coord, medium, jet, eps_e, eps_B);

    Observer obs;

    obs.observeAt(t_bins, coord, f_shock, lumi_dist, z);

    auto syn_e = genSynElectrons(f_shock, p);

    auto syn_ph = genSynPhotons(f_shock, syn_e);

    if (out) {
        write_npz("coord", coord);
        write_npz("shock", f_shock);
        write_npz("syn_e", syn_e);
        write_npz("syn_ph", syn_ph);
        write_npy("t_grid", obs.t_obs_grid, con::sec);
    }

    Array band_pass =
        xt::logspace(std::log10(eVtoHz(band_pass_[0] * con::keV)), std::log10(eVtoHz(band_pass_[1] * con::keV)), 10);

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

    threads.emplace_back(std::thread(lc_gen, "/Users/yihanwang/Projects/afterglow-code-comparison/tests/case1", true));
    threads.emplace_back(std::thread(lc_gen, "/Users/yihanwang/Projects/afterglow-code-comparison/tests/case2", false));
    threads.emplace_back(std::thread(lc_gen, "/Users/yihanwang/Projects/afterglow-code-comparison/tests/case3", false));
    threads.emplace_back(std::thread(lc_gen, "/Users/yihanwang/Projects/afterglow-code-comparison/tests/case4", false));
    threads.emplace_back(std::thread(lc_gen, "/Users/yihanwang/Projects/afterglow-code-comparison/tests/case5", false));

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