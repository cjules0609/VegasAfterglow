//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "shock.h"
/********************************************************************************************************************
 * @brief Constructs a Shock object with the given grid dimensions and energy fractions.
 * @details Initializes various 3D grids for storing physical properties of the shock, including comoving time,
 *          radius, Lorentz factors, magnetic fields, and downstream densities.
 * @param phi_size Number of grid points in phi direction
 * @param theta_size Number of grid points in theta direction
 * @param t_size Number of grid points in time direction
 * @param eps_e Electron energy fraction
 * @param eps_B Magnetic energy fraction
 ********************************************************************************************************************/
Shock::Shock(size_t phi_size, size_t theta_size, size_t t_size, Real eps_e, Real eps_B)
    : t_comv({phi_size, theta_size, t_size}, 0),          // Initialize comoving time grid with 0
      r({phi_size, theta_size, t_size}, 0),               // Initialize radius grid with 0
      theta({phi_size, theta_size, t_size}, 0),           // Initialize theta grid with 0
      Gamma({phi_size, theta_size, t_size}, 1),           // Initialize Gamma grid with 1
      Gamma_rel({phi_size, theta_size, t_size}, 1),       // Initialize Gamma_rel grid with 1
      B({phi_size, theta_size, t_size}, 0),               // Initialize magnetic field grid with 0
      column_num_den({phi_size, theta_size, t_size}, 0),  // Initialize column density grid with 0
      injection_idx({phi_size, theta_size}, t_size),      // Initialize cross index grid with t_size
      required({phi_size, theta_size, t_size}, 1),        // Initialize required grid with true
      eps_e(eps_e),                                       // Set electron energy fraction
      eps_B(eps_B),                                       // Set magnetic energy fraction
      phi_size(phi_size),                                 // Store phi grid size
      theta_size(theta_size),                             // Store theta grid size
      t_size(t_size) {}

/********************************************************************************************************************
 * @brief Resizes all grid components of the Shock object to new dimensions.
 * @param phi_size New number of grid points in phi direction
 * @param theta_size New number of grid points in theta direction
 * @param t_size New number of grid points in time direction
 ********************************************************************************************************************/
void Shock::resize(size_t phi_size, size_t theta_size, size_t t_size) {
    this->phi_size = phi_size;
    this->theta_size = theta_size;
    this->t_size = t_size;
    t_comv.resize({phi_size, theta_size, t_size});
    r.resize({phi_size, theta_size, t_size});
    theta.resize({phi_size, theta_size, t_size});
    Gamma.resize({phi_size, theta_size, t_size});
    Gamma_rel.resize({phi_size, theta_size, t_size});
    B.resize({phi_size, theta_size, t_size});
    column_num_den.resize({phi_size, theta_size, t_size});
    injection_idx.resize({phi_size, theta_size});
    injection_idx.fill(t_size);
    required.resize({phi_size, theta_size, t_size});
    required.fill(1);
}

/********************************************************************************************************************
 * @brief Computes the downstream fluid velocity for a given relative Lorentz factor and magnetization.
 * @param gamma_rel Relative Lorentz factor between upstream and downstream regions
 * @param sigma Magnetization parameter (magnetic to kinetic energy ratio)
 * @return The downstream four-velocity in the shock frame
 ********************************************************************************************************************/
Real compute_downstr_4vel(Real gamma_rel, Real sigma) {
    Real ad_idx = adiabatic_idx(gamma_rel);
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

        // Compute coefficients
        Real A = term1 * gamma_m_1 + 2;
        Real B = -gamma_p_1 * (-ad_idx_m_2 * (ad_idx * gamma_sq + 1) + ad_idx * ad_idx_m_1 * gamma_rel) * sigma -
                 gamma_m_1 * (term1 * (gamma_sq - 2) + 2 * gamma_rel + 3);
        Real C = gamma_p_1 * (ad_idx * (1 - ad_idx / 4) * term2 + 1) * sigma * sigma +
                 term2 * (2 * gamma_rel + ad_idx_m_2 * (ad_idx * gamma_rel - 1)) * sigma +
                 gamma_p_1 * gamma_m_1 * gamma_m_1 * ad_idx_m_1 * ad_idx_m_1;
        Real D = -gamma_m_1 * gamma_p_1 * gamma_p_1 * ad_idx_m_2 * ad_idx_m_2 * sigma * sigma / 4;

        Real b = B / A;
        Real c = C / A;
        Real d = D / A;
        Real P = c - b * b / 3;
        Real Q = 2 * b * b * b / 27 - b * c / 3 + d;
        Real u = std::sqrt(-P / 3);
        Real v = std::clamp(3 * Q / (2 * P * u), -1.0, 1.0);
        Real uds = 2 * u * std::cos((std::acos(v) - 2 * con::pi) / 3) - b / 3;
        return std::sqrt(uds);
    }
}

/********************************************************************************************************************
 * @brief Computes the ratio of upstream to downstream four-velocity across the shock front.
 * @param gamma_rel Relative Lorentz factor between upstream and downstream regions
 * @param sigma Magnetization parameter (magnetic to kinetic energy ratio)
 * @return The ratio of upstream to downstream four-velocity
 ********************************************************************************************************************/
Real compute_4vel_jump(Real gamma_rel, Real sigma) {
    Real u_down_s_ = compute_downstr_4vel(gamma_rel, sigma);
    Real u_up_s_ = compute_upstr_4vel(u_down_s_, gamma_rel);
    Real ratio_u = u_up_s_ / u_down_s_;
    if (u_down_s_ == 0) {
        ratio_u = 4 * gamma_rel;  // (g_hat*gamma_rel+1)/(g_hat-1)
    }
    return ratio_u;
}
