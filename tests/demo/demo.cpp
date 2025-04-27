#include "afterglow.h"

void test_reverse_shock(double xi, double sigma) {
    Real E_iso = 1e51 * con::erg;
    Real theta_c = 0.1;
    Real theta_v = 0;

    Real n_ism = 100 / con::cm3;
    Real eps_e = 1e-2;
    Real eps_B = 1e-4;
    Real p = 2.1;
    Real Gamma0 = 100;
    Real z = 0;

    Array t_obs = xt::logspace(std::log10(0.1 * con::sec), std::log10(1e10 * con::sec), 130);

    ISM medium(n_ism);

    Ejecta jet;

    jet.eps_k = math::tophat(theta_c, E_iso);
    jet.Gamma0 = math::tophat(theta_c, Gamma0);
    jet.sigma0 = math::tophat(theta_c, sigma);

    jet.T0 = calc_engine_duration(E_iso, n_ism, Gamma0, xi);

    size_t t_num = 512;
    size_t theta_num = 2;
    size_t phi_num = 1;

    Coord coord = auto_grid(jet, t_obs, 0.6, theta_v, z, phi_num, theta_num, t_num);

    auto [f_shock, r_shock] = generate_shock_pair(coord, medium, jet, eps_e, eps_B, eps_e, eps_B);
    // auto f_shock = generate_fwd_shock(coord, medium, jet, eps_e, eps_B);
    // auto r_shock = generate_fwd_shock(coord, medium, jet, eps_e, eps_B);

    write_npz("rshock-data/coord" + std::to_string(xi) + "-" + std::to_string(sigma), coord);
    write_npz("rshock-data/f_shock" + std::to_string(xi) + "-" + std::to_string(sigma), f_shock);
    write_npz("rshock-data/r_shock" + std::to_string(xi) + "-" + std::to_string(sigma), r_shock);

    return;
}

int main() {
    double xi[] = {0.001, 0.01, 0.1, 1, 10, 100};
    double sigma[] = {0, 0.0001, 0.01, 1, 100};

    for (auto x : xi) {
        for (auto s : sigma) {
            test_reverse_shock(x, s);
        }
    }
    return 0;
}
