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

/********************************************************************************************************************
 * FUNCTION: soundSpeed
 * DESCRIPTION: Computes the sound speed given the pressure, adiabatic index (ad_idx), and rest mass density (rho_rest).
 *              The formula takes into account both the rest mass energy density (rho_rest * c^2) and the pressure,
 *              scaled by the adiabatic index.
 ********************************************************************************************************************/
double soundSpeed(double pressure, double ad_idx, double rho_rest) {
    // Compute sound speed and multiply by the speed of light (con::c)
    return std::sqrt(ad_idx * pressure / (rho_rest * con::c2 + ad_idx / (ad_idx - 1) * pressure)) * con::c;
}

/********************************************************************************************************************
 * FUNCTION: fa
 * DESCRIPTION: Computes the factor "fa" used in shock physics.
 *              This factor depends on gamma34 (a Lorentz factor ratio), u3s_ (a downstream fluid speed measure),
 *              and the magnetization sigma.
 ********************************************************************************************************************/
double fa(double gamma34, double u3s_, double sigma) {
    // Computes the factor with a correction term proportional to sigma and a denominator involving u3s_ and gamma34.
    return 1 - sigma * (gamma34 + 1) /
                   (u3s_ * u3s_ * gamma34 + u3s_ * std::sqrt((1 + u3s_ * u3s_) * (gamma34 * gamma34 - 1))) / 2;
}

/********************************************************************************************************************
 * FUNCTION: fc
 * DESCRIPTION: Computes the factor "fc" used in shock physics.
 *              It calculates the ratio between p2 (a pressure parameter) and p3 (p2 reduced by pB3, the magnetic
 *pressure).
 ********************************************************************************************************************/
double fc(double p2, double pB3) {
    double p3 = p2 - pB3;  // Effective pressure after subtracting magnetic pressure
    return p2 / p3;        // Ratio of original to effective pressure
}

/********************************************************************************************************************
 * CONSTRUCTOR: Shock::Shock
 * DESCRIPTION: Constructs a Shock object with the given grid dimensions (phi_size, theta_size, r_size)
 *              and energy fractions (eps_e and eps_B). The constructor initializes various 3D grids for storing
 *              comoving time (t_com), engine time (t_eng), relative Lorentz factor (Gamma_rel), magnetic field (B),
 *              and downstream proton column density (column_num_den). Default initial values are provided via
 *create3DGrid.
 ********************************************************************************************************************/
Shock::Shock(size_t phi_size, size_t theta_size, size_t r_size, double eps_e, double eps_B)
    : t_com(create3DGrid(phi_size, theta_size, r_size, 0)),           // Initialize comoving time grid with 0
      t_eng(create3DGrid(phi_size, theta_size, r_size, 0)),           // Initialize engine time grid with 0
      Gamma_rel(create3DGrid(phi_size, theta_size, r_size, 1)),       // Initialize Gamma_rel grid with 1
      B(create3DGrid(phi_size, theta_size, r_size, 0)),               // Initialize magnetic field grid with 0
      column_num_den(create3DGrid(phi_size, theta_size, r_size, 0)),  // Initialize column density grid with 0
      eps_e(eps_e),                                                   // Set electron energy fraction
      eps_B(eps_B),                                                   // Set magnetic energy fraction
      phi_size(phi_size),                                             // Store phi grid size
      theta_size(theta_size),                                         // Store theta grid size
      r_size(r_size) {}                                               // Store radial grid size

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