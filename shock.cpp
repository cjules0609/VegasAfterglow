//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "shock.h"

#include <boost/numeric/odeint.hpp>

#include "afterglow.h"
#include "macros.h"
#include "mesh.h"
#include "physics.h"
#include "utilities.h"

double soundSpeed(double pressure, double ad_idx, double rho_rest) {
    return std::sqrt(ad_idx * pressure / (rho_rest * con::c2 + ad_idx / (ad_idx - 1) * pressure)) * con::c;
}

//----------------------------------------reverse shock---------------------------------------
double fa(double gamma34, double u3s_, double sigma) {
    return 1 - sigma * (gamma34 + 1) /
                   (u3s_ * u3s_ * gamma34 + u3s_ * std::sqrt((1 + u3s_ * u3s_) * (gamma34 * gamma34 - 1))) / 2;
}

double fc(double p2, double pB3) {
    double p3 = p2 - pB3;
    return p2 / p3;
}

//--------------------------------------------------------------------------------------------
Shock::Shock(size_t phi_size, size_t theta_size, size_t r_size, double eps_e, double eps_B)
    : t_com(create3DGrid(phi_size, theta_size, r_size, 0)),
      t_eng(create3DGrid(phi_size, theta_size, r_size, 0)),
      Gamma_rel(create3DGrid(phi_size, theta_size, r_size, 1)),
      B(create3DGrid(phi_size, theta_size, r_size, 0)),
      column_num_den(create3DGrid(phi_size, theta_size, r_size, 0)),
      eps_e(eps_e),
      eps_B(eps_B),
      phi_size(phi_size),
      theta_size(theta_size),
      r_size(r_size) {}

/*
void solve_single_shell(size_t j, Array const& r_b, Array const& r, Shock& f_shock, Shock& r_shock,
                        ForwardShockEqn& eqn) {
    // integrate the shell over r
    for (int k = 0, k1 = 1; stepper.current_time() <= r_b.back();) {
        stepper.do_step(eqn);

        if (eqn.jet.spreading) {  // check if this works for sigma!=0. should be sound speed in region 3;
            double dr = stepper.current_time_step();
            double r_last = stepper.current_time() - dr;
            double Gamma_axis = loglog_interp(r_last, r, f_shock.Gamma[0]);
            if (k == 0) {
                Gamma_axis = eqn.jet.Gamma0_profile(0);
            }

            double ad_idx = adiabatic_index(Gamma_axis);
            double n1 = eqn.medium.rho(r_last) / con::mp;
            double n2 = n_down_str(n1, Gamma_axis, 0);
            double p2 = (ad_idx - 1) * e_thermal_down_str(Gamma_axis, n2);
            double jet_cs = sound_speed(p2, ad_idx, n2 * con::mp);
            if (theta_c < con::pi / 2) {
                theta_c += jet_cs / con::c * dr / (r_last * Gamma_axis * Gamma_axis * theta_c);
                // theta_c += jet_cs / con::c * dr / (r_last * Gamma_axis);
            }
            double t_eng = state[2];
            double dEdOmega_init = eqn.jet.dEdOmega(eqn.theta, t_eng);
            double dEdOmega_spread = eqn.jet.dEdOmega_spread(eqn.theta, theta_c, t_eng);

            if (dEdOmega_init == 0) {
                eqn.spreading_factor = 1;
            } else {
                eqn.spreading_factor = dEdOmega_spread / dEdOmega_init;
            }
        }
    }
}*/