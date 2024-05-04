#include "physics.h"

#include <boost/numeric/odeint.hpp>

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