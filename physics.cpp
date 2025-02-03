//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

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
