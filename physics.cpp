#include "physics.h"

#include <boost/numeric/odeint.hpp>

#include "shock.h"
#include "utilities.h"
double zToLuminosityDistance(double z) {
    using namespace boost::numeric::odeint;
    double atol = 0;
    double rtol = 1e-6;
    auto stepper = make_dense_output(atol, rtol, runge_kutta_dopri5<double>());

    auto eqn = [&](double const& y, double& dydz, double z0) {
        dydz = 1 / std::sqrt(con::Omega_m * (1 + z0) * (1 + z0) * (1 + z0) + con::Omega_L);
    };
    stepper.initialize(0, 0, 1e-6);

    for (; stepper.current_time() <= z;) {
        stepper.do_step(eqn);
    }
    double L = 0;
    stepper.calc_state(z, L);
    return (1 + z) * L * con::c / con::H0;
}

double luminosityDistanceToz(double L) {
    using namespace boost::numeric::odeint;
    double atol = 0;
    double rtol = 1e-6;
    auto stepper = make_dense_output(atol, rtol, runge_kutta_dopri5<double>());

    auto eqn = [&](double const& y, double& dydz, double z0) {
        dydz = 1 / std::sqrt(con::Omega_m * (1 + z0) * (1 + z0) * (1 + z0) + con::Omega_L);
    };
    stepper.initialize(0, 0, 1e-6);
    double L_current = 0;
    for (;;) {
        double z = stepper.current_time();
        double L_current = stepper.current_state();
        if ((1 + z) * L_current * con::c / con::H0 >= L) {
            return z;
        }
        stepper.do_step(eqn);
    }
}

double decRadius(double E_iso, double n_ism, double Gamma0, double engine_dura) {
    return std::max(thinShellDecRadius(E_iso, n_ism, Gamma0), thickShellDecRadius(E_iso, n_ism, Gamma0, engine_dura));
}
double thinShellDecRadius(double E_iso, double n_ism, double Gamma0) {
    return std::pow(3 * E_iso / (4 * con::pi * n_ism * con::mp * con::c2 * Gamma0 * Gamma0), 1.0 / 3);
}
double thickShellDecRadius(double E_iso, double n_ism, double Gamma0, double engine_dura) {
    return std::pow(3 * E_iso * engine_dura * con::c / (4 * con::pi * n_ism * con::mp * con::c2), 0.25);
}

double shellSpreadingRadius(double Gamma0, double engine_dura) { return Gamma0 * Gamma0 * con::c * engine_dura; }

double RSTransitionRadius(double E_iso, double n_ism, double Gamma0, double engine_dura) {
    return std::pow(SedovLength(E_iso, n_ism), 1.5) / std::sqrt(con::c * engine_dura) / Gamma0 / Gamma0;
}

std::tuple<double, double> findRadiusRange(Medium const& medium, Ejecta const& jet, Ejecta const& inj, double t_min,
                                           double t_max, double z) {
    double gamma0 = jet.Gamma0(0, 0, 0);

    double beta0 = std::sqrt(1 - 1 / (gamma0 * gamma0));

    double gamma_min = (gamma0 - 1) / 100 + 1;

    double beta_min = std::sqrt(1 - 1 / (gamma_min * gamma_min));

    double r_min = beta_min * con::c * t_min / (1 + beta_min);

    // find the on-axis r_max by solving theta =0, phi=0 case;
    auto eqn = ForwardShockEqn(medium, jet, inj, 0, 0, 0);
    double r_max = find_r_max(eqn, r_min, t_max);
    return {r_min, r_max};
}

double jetEdge(TernaryFunc const& gamma, double gamma_cut) {
    if (gamma(0, con::pi / 2, 0) > gamma_cut) {
        return con::pi / 2;
    }
    double low = 0;
    double hi = con::pi / 2;
    double eps = 1e-6;
    for (; hi - low > eps;) {
        double mid = 0.5 * (low + hi);
        if (gamma(0, mid, 0) > gamma_cut) {
            low = mid;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (low + hi);
}

Coord adaptiveGrid(Medium const& medium, Ejecta const& jet, Ejecta const& inj, Array const& t_obs, double theta_max,
                   size_t phi_num, size_t theta_num, size_t r_num) {
    double t_max = *std::max_element(t_obs.begin(), t_obs.end());
    double t_min = *std::min_element(t_obs.begin(), t_obs.end());
    auto [r_min, r_max] = findRadiusRange(medium, jet, inj, t_min, t_max);
    Array r = logspace(r_min, r_max, r_num);
    double jet_edge = jetEdge(jet.Gamma0, con::Gamma_cut);
    Array theta = uniform_cos(0, std::min(jet_edge, theta_max), theta_num);
    Array phi = linspace(0, 2 * con::pi, phi_num);
    Coord coord{phi, theta, r};
    coord.t_min = t_min;
    coord.t_max = t_max;
    return coord;
}