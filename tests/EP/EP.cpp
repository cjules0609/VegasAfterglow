
#include <boost/numeric/odeint.hpp>
#include <filesystem>
#include <fstream>

#include "afterglow.h"

void GCN36236(Real theta_c) {
    Real E_iso = 1e53 * con::erg * (0.088 / theta_c) * (0.088 / theta_c);
    Real lumi_dist = 1.2e27 * con::cm;
    Real z = 0.01;
    Real theta_w = 0.6;
    Real Gamma0 = 300;
    Real n_ism = 0.0199526231496888 / con::cm3;
    Real eps_e = 0.1;
    Real eps_B = 0.00019952623149688788;
    Real p = 2.139;

    // create model
    auto medium = createISM(n_ism);

    auto jet = GaussianJet(theta_c, E_iso, Gamma0, 1);

    size_t r_num = 64;
    size_t theta_num = 64;
    size_t phi_num = 64;

    Array t_bins = logspace(1e-1 * con::day, 1e3 * con::day, 100);

    Coord coord = adaptiveGrid(medium, jet, inject::none, t_bins, theta_w, phi_num, theta_num, r_num);

    Shock f_shock = genForwardShock(coord, medium, jet, inject::none, eps_e, eps_B);

    auto syn_e = genSynElectrons(f_shock, p);

    auto syn_ph = genSynPhotons(f_shock, syn_e);

    Array theta_obs = linspace(5 * con::deg, 60 * con::deg, 56);

    Observer obs(coord, f_shock, 0, lumi_dist, z);

    std::array<double, 5> lo = {eVtoHz(0.2 * con::keV)};

    std::array<double, 5> hi = {eVtoHz(10 * con::keV)};

    std::vector<std::string> prefix = {"EP"};

    for (size_t i = 0; i < 1; i++) {
        Array band_pass = logspace(lo[i], hi[i], 5);
        for (auto theta_v : theta_obs) {
            obs.changeViewingAngle(theta_v);

            auto to_suffix = [](Real theta) { return std::to_string(int((theta / con::deg))); };

            char buff[100] = {0};

            sprintf(buff, "data/F_nu_syn_%.0lf_%.0lf_", (theta_v / con::deg), (theta_c / con::deg));

            std::string fname = buff + prefix[i];

            std::cout << fname << std::endl;

            Array F_syn = obs.flux(t_bins, band_pass, syn_ph);

            output(F_syn, fname, con::erg / con::sec / con::cm / con::cm);
        }
    }
    output(t_bins, "data/t_obs", con::sec);
}

int main() {
    Real theta_c[] = {5 * con::deg, 7 * con::deg, 10 * con::deg};

    for (auto theta : theta_c) {
        GCN36236(theta);
    }

    return 0;
}