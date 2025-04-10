//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "shock.h"

#include "forward-shock.hpp"
#include "simple-shock.h"
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
      Gamma(boost::extents[phi_size][theta_size][t_size]),           // Initialize Gamma grid with 0
      Gamma_rel(boost::extents[phi_size][theta_size][t_size]),       // Initialize Gamma_rel grid with 1
      B(boost::extents[phi_size][theta_size][t_size]),               // Initialize magnetic field grid with 0
      column_num_den(boost::extents[phi_size][theta_size][t_size]),  // Initialize column density grid with 0
      injection_idx(boost::extents[phi_size][theta_size]),           // Initialize cross index grid with 0
      eps_e(eps_e),                                                  // Set electron energy fraction
      eps_B(eps_B),                                                  // Set magnetic energy fraction
      phi_size(phi_size),                                            // Store phi grid size
      theta_size(theta_size),                                        // Store theta grid size
      t_size(t_size) {
    std::memset(t_com.data(), 0, t_com.num_elements() * sizeof(Real));
    std::memset(r.data(), 0, r.num_elements() * sizeof(Real));
    std::memset(theta.data(), 0, theta.num_elements() * sizeof(Real));
    std::fill(Gamma.data(), Gamma.data() + Gamma.num_elements(), 1);              // Initialize Gamma to 0
    std::fill(Gamma_rel.data(), Gamma_rel.data() + Gamma_rel.num_elements(), 1);  // Initialize Gamma_rel to 1
    std::memset(B.data(), 0, B.num_elements() * sizeof(Real));
    std::memset(column_num_den.data(), 0, column_num_den.num_elements() * sizeof(Real));
    std::fill(injection_idx.data(), injection_idx.data() + injection_idx.num_elements(), t_size);  
}

// Computes the downstream fluid velocity (u) for a given relative Lorentz factor (gamma_rel) and magnetization (sigma).
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

Real u_UpStr2u_DownStr(Real gamma_rel, Real sigma) {
    Real u_down_s_ = u_DownStr(gamma_rel, sigma);
    Real u_up_s_ = u_UpStr(u_down_s_, gamma_rel);
    Real ratio_u = u_up_s_ / u_down_s_;
    if (u_down_s_ == 0) {
        ratio_u = 4 * gamma_rel;  // (g_hat*gamma_rel+1)/(g_hat-1)
    }
    return ratio_u;
}

/********************************************************************************************************************
 * FUNCTION: genForwardShock
 * DESCRIPTION: Generates a forward shock using the provided coordinates, medium, jet, and energy
 *              fractions. It creates a Shock object and iterates over phi, theta values, solving the shock
 *              evolution for each theta slice.
 ********************************************************************************************************************/
Shock genForwardShock(Coord const& coord, Medium const& medium, Ejecta const& jet, Real eps_e, Real eps_B, Real rtol,
                      bool is_axisymmetric) {
    auto [phi_size, theta_size, t_size] = coord.shape();  // Unpack coordinate dimensions
    size_t phi_size_needed = is_axisymmetric ? 1 : phi_size;
    Shock f_shock(phi_size_needed, theta_size, t_size, eps_e, eps_B);

    for (size_t i = 0; i < phi_size_needed; ++i) {
        Real theta_s =
            jetSpreadingEdge(jet, medium, coord.phi[i], coord.theta[0], coord.theta[theta_size - 1], coord.t[0]);
        for (size_t j = 0; j < theta_size; ++j) {
            // Create a ForwardShockEqn for each theta slice
            // auto eqn = ForwardShockEqn(medium, jet, coord.phi[i], coord.theta[j], eps_e, theta_s);
            auto eqn =
                SimpleShockEqn(medium, jet, coord.phi[i], j ? coord.theta[j - 1] : 0, coord.theta[j], eps_e, theta_s);
            //   Solve the shock shell for this theta slice
            solveForwardShell(i, j, coord.t, f_shock, eqn, rtol);
        }
    }

    return f_shock;
}

/********************************************************************************************************************
 * FUNCTION: genFRShocks
 * DESCRIPTION: Generates a pair of forward and reverse shocks using the provided coordinates, medium, jet,
 *              and energy fractions. It creates two Shock objects (one forward, one reverse) and solves
 *              the shock shells for each phi, theta slice.
 ********************************************************************************************************************/
ShockPair genFRShocks(Coord const& coord, Medium const& medium, Ejecta const& jet, Real eps_e_f, Real eps_B_f,
                      Real eps_e_r, Real eps_B_r, Real rtol, bool is_axisymmetric) {
    auto [phi_size, theta_size, t_size] = coord.shape();
    size_t phi_size_needed = is_axisymmetric ? 1 : phi_size;
    Shock f_shock(phi_size_needed, theta_size, t_size, eps_e_f, eps_B_f);
    Shock r_shock(phi_size_needed, theta_size, t_size, eps_e_r, eps_B_r);
    for (size_t i = 0; i < phi_size_needed; ++i) {
        Real theta_s =
            jetSpreadingEdge(jet, medium, coord.phi[i], coord.theta[0], coord.theta[theta_size - 1], coord.t[0]);
        for (size_t j = 0; j < theta_size; ++j) {
            // Create equations for forward and reverse shocks for each theta slice
            // auto eqn_f = ForwardShockEqn(medium, jet, inject, coord.phi[i], coord.theta[j], eps_e);
            auto eqn_f =
                SimpleShockEqn(medium, jet, coord.phi[i], j ? coord.theta[j - 1] : 0, coord.theta[j], eps_e_f, theta_s);
            auto eqn_r = FRShockEqn(medium, jet, coord.phi[i], coord.theta[j], eps_e_r);
            // Solve the forward-reverse shock shell
            solveFRShell(i, j, coord.t, f_shock, r_shock, eqn_f, eqn_r, rtol);
        }
    }

    return std::make_pair(std::move(f_shock), std::move(r_shock));
}