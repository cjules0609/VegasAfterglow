#include "afterglow.h"

#include <boost/numeric/odeint.hpp>

void afterglow_gen() {
    std::string prefix = "off-axis-tophat/";
    double E_iso = 1e53 * con::erg;
    double Gamma0 = 300;
    double n_ism = 1 / con::cm / con::cm / con::cm;
    double eps_rad = 1;
    double eps_e = 0.1;
    double eps_B = 0.01;
    double pel = 2.3;
    double theta_j = 15.0 * con::deg;

    // create model
    auto medium = create_ISM(n_ism, eps_e, eps_B, pel);
    auto inj = createIsoPowerLawInjection(0 / 9e20 / 2e33, 1000, 1, 2);
    auto jet = createTopHatJet(E_iso, Gamma0, theta_j, inj);
    // auto blast = createPowerLawJet(E_iso, Gamma0, theta_j, 4, inj);
    // auto blast = createGaussianJet(E_iso, Gamma0, theta_j / 6, inj);

    // generate grid
    double M0 = E_iso / (Gamma0 * con::c * con::c);
    double R_ES = pow(3 * M0 / (4 * con::pi * n_ism * con::mp * Gamma0), 1.0 / 3);
    size_t grid_num = 300;
    Array r = logspace(R_ES / 100, R_ES * 100, grid_num);
    Array theta = linspace(0, theta_j, 30);
    Array phi = linspace(0, 2 * con::pi, 37);
    Coord coord{r, theta, phi};

    write2file(coord.r, prefix + "r");
    write2file(coord.theta, prefix + "theta");
    write2file(coord.phi, prefix + "phi");

    // solve dynamics
    auto [t_com_history, Gamma_history] = solve_dynamics(coord, jet, medium);
    write2file(Gamma_history, prefix + "Gamma");

    auto syn_rad_history = calc_syn_radiation(coord, Gamma_history, t_com_history, medium);
    write2file(syn_rad_history, prefix + "syn");

    Observer observer;

    double theta_obs = 10.0 * con::deg;
    observer.observe(coord, Gamma_history, theta_obs);
    write2file(observer.t_obs, prefix + "t_obs");
    write2file(observer.Doppler, prefix + "doppler");

    // specify observables
    double nu_obs_R = 1e9 * 500;
    double nu_obs_O = 1e14 * 500;
    double nu_obs_X = 1e17 * 500;
    Array nu_obs{nu_obs_R, nu_obs_O, nu_obs_X};

    size_t data_points = 100;

    for (size_t i = 0; i < nu_obs.size(); ++i) {
        MeshGrid3d I_syn_obs = observer.I_nu_history(nu_obs[i], syn_rad_history);
        MeshGrid L_nu = observer.light_curve(data_points, nu_obs[i], syn_rad_history);

        write2file(I_syn_obs, prefix + "I_nu_" + std::to_string(int(log10(nu_obs[i] / 500))));
        write2file(L_nu, prefix + "L_nu_" + std::to_string(int(log10(nu_obs[i] / 500))));
    }

    // write to file
}

int main() {
    afterglow_gen();
    return 0;
}