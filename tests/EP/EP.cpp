
#include <boost/numeric/odeint.hpp>
#include <filesystem>
#include <fstream>

#include "afterglow.h"

void GCN36236(Real theta_c) {
    Real E_iso = 1e53 * con::erg * (0.088 / theta_c) * (0.088 / theta_c);
    Real lumi_dist = 6.6e26 * con::cm;
    Real z = 0.009;
    Real theta_w = 0.6;
    Real Gamma0 = 300;
    Real n_ism = 0.0199526231496888 / con::cm3;
    Real eps_e = 0.1;
    Real eps_B = 0.00019952623149688788;
    Real p = 2.139;

    // create model
    auto medium = createISM(n_ism);

    auto jet = GaussianJet(theta_c, E_iso, Gamma0);

    size_t r_num = 500;
    size_t theta_num = 150;
    size_t phi_num = 50;

    Array t_bins = logspace(1e-1 * con::day, 1e3 * con::day, 100);

    Coord coord = adaptiveGrid(medium, jet, inject::none, t_bins, theta_w, phi_num, theta_num, r_num);

    // solve dynamics
    // Shock f_shock(coord, eps_e, eps_B, p);

    // solve_shocks(coord, jet, medium, f_shock);
    Shock f_shock = genForwardShock(coord, medium, jet, inject::none, eps_e, eps_B);

    auto syn_e = genSynElectrons(f_shock, p);

    auto syn_ph = genSynPhotons(f_shock, syn_e);

    Array theta_obs = linspace(5 * con::deg, 60 * con::deg, 56);

    Observer obs(coord, f_shock, 0, lumi_dist, z);

    for (auto theta_v : theta_obs) {
        obs.changeViewingAngle(theta_v);

        Array band_pass = logspace(eVtoHz(0.3 * con::keV), eVtoHz(10 * con::keV), 5);

        auto to_suffix = [](Real theta) { return std::to_string(int((theta / con::deg))); };

        Array F_syn = obs.flux(t_bins, band_pass, syn_ph);

        char buff[100] = {0};

        sprintf(buff, "ep/F_nu_syn_%.0lf_%.0lf", (theta_v / con::deg), (theta_c / con::deg));

        std::string fname = buff;

        output(F_syn, fname, con::erg / con::sec / con::cm / con::cm);
        std::cout << theta_v / con::deg << std::endl;
    }

    output(boundaryToCenterLog(t_bins), "ep/t_obs", con::sec);
}

auto solve_u2s1(Real sigma, Real gamma_max, size_t size) {
    Array gamma_rel = logspace(1e-5, gamma_max, size);
    Array u2s = zeros(size);

    for (size_t i = 0; i < size; ++i) {
        Real gamma = gamma_rel[i] + 1;
        Real ad_idx = adiabaticIndex(gamma);
        Real A = ad_idx * (2 - ad_idx) * (gamma - 1) + 2;
        Real B = -(gamma + 1) * ((2 - ad_idx) * (ad_idx * gamma * gamma + 1) + ad_idx * (ad_idx - 1) * gamma) * sigma -
                 (gamma - 1) * (ad_idx * (2 - ad_idx) * (gamma * gamma - 2) + 2 * gamma + 3);
        Real C = (gamma + 1) * (ad_idx * (1 - ad_idx / 4) * (gamma * gamma - 1) + 1) * sigma * sigma +
                 (gamma * gamma - 1) * (2 * gamma - (2 - ad_idx) * (ad_idx * gamma - 1)) * sigma +
                 (gamma + 1) * (gamma - 1) * (gamma - 1) * (ad_idx - 1) * (ad_idx - 1);
        Real D = -(gamma - 1) * (gamma + 1) * (gamma + 1) * (2 - ad_idx) * (2 - ad_idx) * sigma * sigma / 4;

        Real x0 = (-B - sqrt(B * B - 3 * A * C)) / 3 / A;
        Real x1 = (-B + sqrt(B * B - 3 * A * C)) / 3 / A;
        u2s[i] =
            sqrt(rootBisection([=](Real x) -> Real { return A * x * x * x + B * x * x + C * x + D; }, x0, x1, 1e-13));

        Real exp = sqrt((gamma - 1) * (ad_idx - 1) * (ad_idx - 1) / (ad_idx * (2 - ad_idx) * (gamma - 1) + 2));
        Real y = exp * exp;
        Real z = u2s[i] * u2s[i];
        std::cout << gamma << ", " << ad_idx << ", " << A << ", " << B << ", " << C << ", " << D << ", " << u2s[i]
                  << ", " << sqrt(gamma * gamma - 1) << ", " << exp << "," << A * y * y * y + B * y * y + C * y + D
                  << "," << A * z * z * z + B * z * z + C * z + D << std::endl;
    }

    return std::make_pair(gamma_rel, u2s);
}

int main() {
    Array theta_c = linspace(1 * con::deg, 30 * con::deg, 30);

    for (auto theta : theta_c) {
        GCN36236(theta);
    }

    return 0;
}