
#include "afterglow.h"
void test_reverse_shock(double xi, double sigma) {
    Real E_iso = 1e49 * unit::erg;
    Real theta_c = 0.1;
    Real theta_v = 0;

    Real n_ism = 100 / unit::cm3;
    Real eps_e = 1e-2;
    Real eps_B = 1e-4;
    Real Gamma0 = 200;
    Real z = 0;

    Array t_obs = xt::logspace(std::log10(0.03 * unit::sec), std::log10(1e9 * unit::sec), 130);

    ISM medium(n_ism);

    Ejecta jet;

    jet.eps_k = math::tophat(theta_c, E_iso);
    jet.Gamma0 = math::tophat(theta_c, Gamma0);
    jet.sigma0 = math::tophat(theta_c, sigma);

    jet.T0 = calc_engine_duration(E_iso, n_ism, Gamma0, xi);
    std::cout << "T0: " << jet.T0 / unit::sec << ' ' << xi << ' ' << sigma << std::endl;

    Coord coord = auto_grid(jet, t_obs, 0.6, theta_v, z, 0.001, 1, 50);

    auto [f_shock, r_shock] = generate_shock_pair(coord, medium, jet, eps_e, eps_B, eps_e, eps_B);
    // auto f_shock = generate_fwd_shock(coord, medium, jet, eps_e, eps_B);
    // auto r_shock = generate_fwd_shock(coord, medium, jet, eps_e, eps_B);

    write_npz("rshock-data/coord" + std::to_string(xi) + "-" + std::to_string(sigma), coord);
    write_npz("rshock-data/f_shock" + std::to_string(xi) + "-" + std::to_string(sigma), f_shock);
    write_npz("rshock-data/r_shock" + std::to_string(xi) + "-" + std::to_string(sigma), r_shock);

    return;
}

void test_spreading() {
    Real E_iso = 1e51 * unit::erg;
    Real theta_c = 10 * unit::deg;
    Real theta_v = 0;

    Real n_ism = 1 / unit::cm3;
    Real eps_e = 1e-2;
    Real eps_B = 1e-3;
    Real Gamma0 = 300;
    Real z = 0;

    Array t_obs = xt::logspace(std::log10(0.1 * unit::sec), std::log10(1e8 * unit::sec), 130);

    ISM medium(n_ism);

    Ejecta jet;

    jet.eps_k = math::tophat(theta_c, E_iso);
    jet.Gamma0 = math::tophat(theta_c, Gamma0);

    jet.spreading = true;

    Coord coord = auto_grid(jet, t_obs, con::pi / 2, theta_v, z);

    auto shock = generate_fwd_shock(coord, medium, jet, eps_e, eps_B);

    write_npz("spreading-data/shock", shock);
}

void test_grid() {
    Real E_iso = 1e51 * unit::erg;
    Real theta_c = 10 * unit::deg;
    Real theta_v = 20 * unit::deg;

    Real n_ism = 1 / unit::cm3;
    Real eps_e = 1e-2;
    Real eps_B = 1e-3;
    Real Gamma0 = 300;
    Real z = 0;

    Array t_obs = xt::logspace(std::log10(0.1 * unit::sec), std::log10(1e8 * unit::sec), 130);

    ISM medium(n_ism);

    Ejecta jet;

    jet.eps_k = math::tophat(theta_c, E_iso);
    jet.Gamma0 = math::tophat(theta_c, Gamma0);

    // jet.eps_k = math::gaussian(theta_c, E_iso);
    // jet.Gamma0 = math::gaussian(theta_c, Gamma0);

    jet.spreading = true;

    Coord coord = auto_grid(jet, t_obs, con::pi / 2, theta_v, z);

    write_npz("coord", coord);
}

int main() {
    double xi[] = {0.001, 0.01, 0.1, 1, 2, 3, 5, 10, 100};
    double sigma[] = {0, 0.01, 0.05, 0.1, 1, 100};

    for (auto x : xi) {
        for (auto s : sigma) {
            test_reverse_shock(x, s);
        }
    }

    test_spreading();

    return 0;
}
