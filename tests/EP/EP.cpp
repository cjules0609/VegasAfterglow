
#include <boost/numeric/odeint.hpp>
#include <filesystem>
#include <fstream>

#include "afterglow.h"

void GCN36236(double theta_c) {
    double E_iso = 1e53 * con::erg * (0.088 / theta_c) * (0.088 / theta_c);
    double lumi_dist = 6.6e26 * con::cm;
    double z = 0.009;
    double theta_w = 0.6;
    double Gamma0 = 300;
    double n_ism = 0.0199526231496888 / con::cm3;
    double eps_e = 0.1;
    double eps_B = 0.00019952623149688788;
    double p = 2.139;

    // create model
    auto medium = create_ISM(n_ism);
    auto jet = GaussianJet(theta_c, E_iso, Gamma0, 1);

    size_t r_num = 500;
    size_t theta_num = 150;
    size_t phi_num = 50;

    double R_dec = dec_radius(E_iso, n_ism, Gamma0, jet.duration);

    auto r = logspace(R_dec / 100, R_dec * 100, r_num);
    auto theta = adaptive_theta_space(theta_num, jet.Gamma0_profile, 0.6);
    auto phi = linspace(0, 2 * con::pi, phi_num);

    ;
    Coord coord{r, theta, phi};

    // solve dynamics
    // Shock f_shock(coord, eps_e, eps_B, p);

    // solve_shocks(coord, jet, medium, f_shock);
    Shock f_shock = gen_forward_shock(coord, jet, medium, eps_e, eps_B, p);
    auto syn_e = gen_syn_electrons(coord, f_shock);
    auto syn_ph = gen_syn_photons(syn_e, coord, f_shock);

    Array t_bins = logspace(1e-1 * con::day, 1e3 * con::day, 100);

    Array theta_obs = linspace(5 * con::deg, 60 * con::deg, 56);

    Observer obs(coord);

    for (auto theta_v : theta_obs) {
        obs.observe(f_shock.Gamma, theta_v, lumi_dist, z);

        Array band_pass = logspace(eVtoHz(0.3 * con::keV), eVtoHz(10 * con::keV), 5);

        auto to_suffix = [](double theta) { return std::to_string(int((theta / con::deg))); };

        Array F_syn = obs.flux(t_bins, band_pass, syn_ph);

        char buff[100] = {0};

        sprintf(buff, "ep/F_nu_syn_%.0lf_%.0lf", (theta_v / con::deg), (theta_c / con::deg));

        std::string fname = buff;

        output(F_syn, fname, con::erg / con::sec / con::cm / con::cm);
        std::cout << theta_v / con::deg << std::endl;
    }

    output(boundary2centerlog(t_bins), "ep/t_obs", con::sec);
}

auto solve_u2s1(double sigma, double gamma_max, size_t size) {
    Array gamma_rel = logspace(1e-5, gamma_max, size);
    Array u2s = zeros(size);

    for (size_t i = 0; i < size; ++i) {
        double gamma = gamma_rel[i] + 1;
        double ad_idx = adiabatic_index(gamma);
        double A = ad_idx * (2 - ad_idx) * (gamma - 1) + 2;
        double B =
            -(gamma + 1) * ((2 - ad_idx) * (ad_idx * gamma * gamma + 1) + ad_idx * (ad_idx - 1) * gamma) * sigma -
            (gamma - 1) * (ad_idx * (2 - ad_idx) * (gamma * gamma - 2) + 2 * gamma + 3);
        double C = (gamma + 1) * (ad_idx * (1 - ad_idx / 4) * (gamma * gamma - 1) + 1) * sigma * sigma +
                   (gamma * gamma - 1) * (2 * gamma - (2 - ad_idx) * (ad_idx * gamma - 1)) * sigma +
                   (gamma + 1) * (gamma - 1) * (gamma - 1) * (ad_idx - 1) * (ad_idx - 1);
        double D = -(gamma - 1) * (gamma + 1) * (gamma + 1) * (2 - ad_idx) * (2 - ad_idx) * sigma * sigma / 4;

        double x0 = (-B - sqrt(B * B - 3 * A * C)) / 3 / A;
        double x1 = (-B + sqrt(B * B - 3 * A * C)) / 3 / A;
        u2s[i] = sqrt(
            root_bisection([=](double x) -> double { return A * x * x * x + B * x * x + C * x + D; }, x0, x1, 1e-13));

        double exp = sqrt((gamma - 1) * (ad_idx - 1) * (ad_idx - 1) / (ad_idx * (2 - ad_idx) * (gamma - 1) + 2));
        double y = exp * exp;
        double z = u2s[i] * u2s[i];
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