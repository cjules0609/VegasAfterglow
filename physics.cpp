#include "physics.h"

#include <boost/numeric/odeint.hpp>

#include "utilities.h"
double z_to_luminosity_distance(double z) {
    using namespace boost::numeric::odeint;
    double atol = 0;
    double rtol = 1e-6;
    auto stepper = make_dense_output(atol, rtol, runge_kutta_dopri5<double>());

    auto eqn = [&](double const& y, double& dydz, double z0) {
        dydz = 1 / sqrt(con::Omega_m * (1 + z0) * (1 + z0) * (1 + z0) + con::Omega_L);
    };
    stepper.initialize(0, 0, 1e-6);

    for (; stepper.current_time() <= z;) {
        stepper.do_step(eqn);
    }
    double L = 0;
    stepper.calc_state(z, L);
    return (1 + z) * L * con::c / con::H0;
}

double luminosity_distance_to_z(double L) {
    using namespace boost::numeric::odeint;
    double atol = 0;
    double rtol = 1e-6;
    auto stepper = make_dense_output(atol, rtol, runge_kutta_dopri5<double>());

    auto eqn = [&](double const& y, double& dydz, double z0) {
        dydz = 1 / sqrt(con::Omega_m * (1 + z0) * (1 + z0) * (1 + z0) + con::Omega_L);
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

double dec_radius(double E_iso, double n_ism, double Gamma0, double engine_dura) {
    return std::max(thin_shell_dec_radius(E_iso, n_ism, Gamma0),
                    thick_shell_dec_radius(E_iso, n_ism, Gamma0, engine_dura));
}
double thin_shell_dec_radius(double E_iso, double n_ism, double Gamma0) {
    return pow(3 * E_iso / (4 * con::pi * n_ism * con::mp * con::c2 * Gamma0 * Gamma0), 1.0 / 3);
}
double thick_shell_dec_radius(double E_iso, double n_ism, double Gamma0, double engine_dura) {
    return pow(3 * E_iso * engine_dura * con::c / (4 * con::pi * n_ism * con::mp * con::c2), 0.25);
}

double shell_spreading_radius(double Gamma0, double engine_dura) { return Gamma0 * Gamma0 * con::c * engine_dura; }

double RS_transition_radius(double E_iso, double n_ism, double Gamma0, double engine_dura) {
    return pow(Sedov_length(E_iso, n_ism), 1.5) / sqrt(con::c * engine_dura) / Gamma0 / Gamma0;
}

UDownStr::UDownStr(double sigma) : sigma(sigma), gamma_1(logspace(1e-6, 1e4, 100)), u2s(zeros(100)) {
    for (size_t i = 0; i < gamma_1.size(); ++i) {
        double gamma = gamma_1[i] + 1;
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
            root_bisection([=](double x) -> double { return A * x * x * x + B * x * x + C * x + D; }, x0, x1, 1e-9));
    }
}

double UDownStr::interp(double gamma) const { return loglog_interp_eq_spaced(gamma - 1, gamma_1, u2s, false, false); }