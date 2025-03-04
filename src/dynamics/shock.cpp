//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "shock.h"

#include "afterglow.h"
#include "macros.h"
#include "mesh.h"
#include "physics.h"
#include "utilities.h"

/********************************************************************************************************************
 * CONSTRUCTOR: Shock::Shock
 * DESCRIPTION: Constructs a Shock object with the given grid dimensions (phi_size, theta_size, t_size)
 *              and energy fractions (eps_e and eps_B). The constructor initializes various 3D grids for storing
 *              comoving time (t_com), engine time (t_eng), relative Lorentz factor (Gamma_rel), magnetic field (B),
 *              and downstream proton column density (column_num_den). Default initial values are provided via
 *create3DGrid.
 ********************************************************************************************************************/
Shock::Shock(size_t phi_size, size_t theta_size, size_t t_size, Real eps_e, Real eps_B)
    : t_com(boost::extents[phi_size][theta_size][t_size]),           // Initialize comoving time grid with 0
      r(boost::extents[phi_size][theta_size][t_size]),               // Initialize engine time grid with 0
      theta(boost::extents[phi_size][theta_size][t_size]),           // Initialize theta grid with 0
      Gamma_rel(boost::extents[phi_size][theta_size][t_size]),       // Initialize Gamma_rel grid with 1
      B(boost::extents[phi_size][theta_size][t_size]),               // Initialize magnetic field grid with 0
      column_num_den(boost::extents[phi_size][theta_size][t_size]),  // Initialize column density grid with 0
      eps_e(eps_e),                                                  // Set electron energy fraction
      eps_B(eps_B),                                                  // Set magnetic energy fraction
      phi_size(phi_size),                                            // Store phi grid size
      theta_size(theta_size),                                        // Store theta grid size
      t_size(t_size) {
    std::memset(t_com.data(), 0, t_com.num_elements() * sizeof(Real));
    std::memset(r.data(), 0, r.num_elements() * sizeof(Real));
    std::memset(theta.data(), 0, theta.num_elements() * sizeof(Real));
    std::fill(Gamma_rel.data(), Gamma_rel.data() + Gamma_rel.num_elements(), 1);  // Initialize Gamma_rel to 1
    std::memset(B.data(), 0, B.num_elements() * sizeof(Real));
    std::memset(column_num_den.data(), 0, column_num_den.num_elements() * sizeof(Real));
}

// Computes the downstream fluid velocity (u) for a given relative Lorentz factor (gamma_rel) and magnetization (sigma).
// This function uses the adiabatic index computed from gamma_rel and applies different formulas based on sigma.
Real u_DownStr(Real gamma_rel, Real sigma) {
    Real ad_idx = adiabaticIndex(gamma_rel);
    Real gamma_m_1 = gamma_rel - 1;  // (gamma_rel - 1)
    Real ad_idx_m_2 = ad_idx - 2;    // (ad_idx - 2)
    Real ad_idx_m_1 = ad_idx - 1;    // (ad_idx - 1)
    if (sigma == 0) {
        return std::sqrt(gamma_m_1 * ad_idx_m_1 * ad_idx_m_1 / (-ad_idx * ad_idx_m_2 * gamma_m_1 + 2));
    } else {
        Real gamma_sq = gamma_rel * gamma_rel;  // gamma_rel^2
        Real gamma_p_1 = gamma_rel + 1;         // (gamma_rel + 1)

        // Precompute common terms
        Real term1 = -ad_idx * ad_idx_m_2;
        Real term2 = gamma_sq - 1;
        Real term3 = gamma_sq - 2;
        Real term4 = gamma_p_1 * gamma_m_1;

        // Compute coefficients
        Real A = term1 * gamma_m_1 + 2;
        Real B = -gamma_p_1 * (-ad_idx_m_2 * (ad_idx * gamma_sq + 1) + ad_idx * ad_idx_m_1 * gamma_rel) * sigma -
                 gamma_m_1 * (term1 * term3 + 2 * gamma_rel + 3);
        Real C = gamma_p_1 * (ad_idx * (1 - ad_idx / 4) * term2 + 1) * sigma * sigma +
                 term2 * (2 * gamma_rel + ad_idx_m_2 * (ad_idx * gamma_rel - 1)) * sigma +
                 term4 * gamma_m_1 * ad_idx_m_1 * ad_idx_m_1;
        Real D = -gamma_m_1 * gamma_p_1 * gamma_p_1 * ad_idx_m_2 * ad_idx_m_2 * sigma * sigma / 4;

        Real b = B / A;
        Real c = C / A;
        Real d = D / A;
        Real P = c - b * b / 3;
        Real Q = 2 * b * b * b / 27 - b * c / 3 + d;
        Real u = std::sqrt(-P / 3);
        Real uds = 2 * u * std::cos(std::acos((3 * Q / (2 * P * u)) - 2 * con::pi) / 3) - b / 3;
        return std::sqrt(uds);
    }
}

Real u_UpStr2u_DownStr(Real gamma_rel, Real sigma) {
    Real u_down_s_ = u_DownStr(gamma_rel, sigma);
    Real u_up_s_ = u_UpStr(u_down_s_, gamma_rel);
    Real ratio_u = u_up_s_ / u_down_s_;
    if (u_down_s_ == 0) {
        ratio_u = (7 * gamma_rel + 1) / (gamma_rel + 1);  // (g_hat+1)/(g_hat-1)
    }
    return ratio_u;
}

Real soundSpeed(Real pressure, Real ad_idx, Real rho_rest) {
    // Compute sound speed and multiply by the speed of light (con::c)
    return std::sqrt(ad_idx * pressure / (rho_rest * con::c2 + ad_idx / (ad_idx - 1) * pressure)) * con::c;
}

void updateShockState(Shock& shock, size_t i, size_t j, size_t k, Real r, Real theta, Real Gamma_rel, Real t_com,
                      Real dNdOmega_up, Real n_up_str, Real sigma) {
    if (Gamma_rel > 1) {
        Real ratio_u = u_UpStr2u_DownStr(Gamma_rel, sigma);
        Real pB_up = calc_pB4(n_up_str, sigma);
        Real pB_down = pB_up * ratio_u * ratio_u;
        Real n_down_str = n_up_str * ratio_u;
        Real co_moving_width = dNdOmega_up / (r * r * n_down_str);
        Real e_th = e_ThermalDownStr(Gamma_rel, n_down_str);
        shock.Gamma_rel[i][j][k] = Gamma_rel;
        shock.t_com[i][j][k] = t_com;
        shock.r[i][j][k] = r;
        shock.theta[i][j][k] = theta;
        shock.column_num_den[i][j][k] = n_down_str * co_moving_width;
        shock.B[i][j][k] = coMovingWeibelB(shock.eps_B, e_th) + std::sqrt(pB_down * 8 * con::pi);
    } else {
        shock.Gamma_rel[i][j][k] = 1;
        shock.t_com[i][j][k] = 0;
        shock.r[i][j][k] = r;
        shock.theta[i][j][k] = theta;
        shock.column_num_den[i][j][k] = 0;
        shock.B[i][j][k] = 0;
    }
}
/*
void solve_single_shell(size_t j, Array const& r_b, Array const& r, Shock& f_shock, Shock& r_shock,
                        ForwardShockEqn& eqn) {
    // integrate the shell over r
    for (int k = 0, k1 = 1; stepper.current_time() <= r_b.back();) {
        stepper.do_step(eqn);

        if (eqn.jet.spreading) {  // check if this works for sigma!=0. should be sound speed in region 3;
            Real dr = stepper.current_time_step();
            Real r_last = stepper.current_time() - dr;
            Real Gamma_axis = loglog_interp(r_last, r, f_shock.Gamma[0]);
            if (k == 0) {
                Gamma_axis = eqn.jet.Gamma0_profile(0);
            }

            Real ad_idx = adiabatic_index(Gamma_axis);
            Real n1 = eqn.medium.rho(r_last) / con::mp;
            Real n2 = n_down_str(n1, Gamma_axis, 0);
            Real p2 = (ad_idx - 1) * e_thermal_down_str(Gamma_axis, n2);
            Real jet_cs = sound_speed(p2, ad_idx, n2 * con::mp);
            if (theta_c < con::pi / 2) {
                theta_c += jet_cs / con::c * dr / (r_last * Gamma_axis * Gamma_axis * theta_c);
                // theta_c += jet_cs / con::c * dr / (r_last * Gamma_axis);
            }
            Real t_eng = state[2];
            Real dEdOmega_init = eqn.jet.dEdOmega(eqn.theta, t_eng);
            Real dEdOmega_spread = eqn.jet.dEdOmega_spread(eqn.theta, theta_c, t_eng);

            if (dEdOmega_init == 0) {
                eqn.spreading_factor = 1;
            } else {
                eqn.spreading_factor = dEdOmega_spread / dEdOmega_init;
            }
        }
    }
}*/
