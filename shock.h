//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _FSDYNAMICS_
#define _FSDYNAMICS_

#include <boost/numeric/odeint.hpp>
#include <tuple>

#include "jet.h"
#include "medium.h"
#include "mesh.h"
#include "physics.h"

/********************************************************************************************************************
 * CLASS: Shock
 * DESCRIPTION: Represents a shock structure that stores grid-based data for the shock evolution, including
 *              comoving time, engine time, relative Lorentz factor, magnetic field, and proton column density.
 *              It also stores constant energy fractions (eps_e and eps_B) and provides a method to return
 *              the grid dimensions.
 ********************************************************************************************************************/
class Shock {
   public:
    Shock(size_t phi_size, size_t theta_size, size_t r_size, double eps_e, double eps_B);
    Shock();

    MeshGrid3d t_com;           // comoving time
    MeshGrid3d t_eng;           // engine time
    MeshGrid3d Gamma_rel;       // relative lorentz factor between down stream and up stream
    MeshGrid3d B;               // comoving magnetic field
    MeshGrid3d column_num_den;  // down stream proton column number density
    double eps_e{0};            // electron energy fraction
    double eps_B{0};            // magnetic energy fraction

    auto shape() const { return std::make_tuple(phi_size, theta_size, r_size); }  // Returns grid dimensions

   private:
    size_t phi_size{0};    // Number of grid points in phi direction
    size_t theta_size{0};  // Number of grid points in theta direction
    size_t r_size{0};      // Number of grid points in radial direction
};

/********************************************************************************************************************
 * CLASS: ForwardShockEqn
 * DESCRIPTION: Represents the forward shock equation for a given Jet and Injector. It defines a state vector
 *              (an array of 5 doubles) and overloads operator() to compute the derivatives of the state with
 *              respect to radius r. It also declares helper functions for the derivatives.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
class ForwardShockEqn {
   public:
    using State = std::array<double, 5>;  // State vector: typically [Gamma, u, t_eng, t_com, D_jet]

    ForwardShockEqn(Medium const& medium, Jet const& jet, Injector const& inject, double phi, double theta,
                    double eps_e);

    Medium const& medium;        // Reference to the medium properties
    Jet const& jet;              // Reference to the jet properties
    Injector const& inject;      // Reference to the injector properties
    double const phi{0};         // Angular coordinate phi
    double const theta{0};       // Angular coordinate theta
    double const eps_e{0};       // Electron energy fraction
    double const jet_sigma{0};   // Jet magnetization parameter
    double gamma4{1};            // Initial Lorentz factor (or a related parameter)
    double spreading_factor{1};  // Factor to account for jet spreading

    // Overloaded operator() to compute the derivatives of the state vector with respect to radius r.
    void operator()(State const& y, State& dydr, double r);

   private:
    // Helper function: computes the derivative of Gamma with respect to r.
    inline double dGammadr(double r, double Gamma, double u, double t_eng, double ad_idx, double rho, double dtdr);
    // Helper function: computes the derivative of u with respect to r.
    inline double dUdr(double r, double Gamma, double u, double t_eng, double ad_idx, double rho, double dGdr);

    double const jet_Gamma0{0};  // Initial Gamma from the jet model
    double const inj_Gamma0{0};  // Initial Gamma from the injector
    double const inj_sigma{0};   // Injector magnetization parameter
    double const dM0{0};         // Initial mass per unit solid angle
};

/********************************************************************************************************************
 * CLASS: FRShockEqn
 * DESCRIPTION: Represents the reverse shock (or forward-reverse shock) equation for a given Jet and Injector.
 *              It defines a state vector (an array of 5 doubles) and overloads operator() to compute the
 *              derivatives of the state with respect to radius r. It also declares a helper function to compute
 *              the derivative of N3 (number per solid angle) with respect to r.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
class FRShockEqn {
   public:
    using State = std::array<double, 5>;  // State vector for reverse shock variables

    FRShockEqn(Medium const& medium, Jet const& jet, Injector const& inject, double phi, double theta);

    Medium const& medium;       // Reference to the medium properties
    Jet const& jet;             // Reference to the jet properties
    Injector const& inject;     // Reference to the injector properties
    double const phi{0};        // Angular coordinate phi
    double const theta{0};      // Angular coordinate theta
    double const jet_sigma{0};  // Jet magnetization parameter
    double gamma4{1};           // Initial Gamma parameter from the jet

    // Overloaded operator() to compute the derivatives of the state vector with respect to radius r.
    void operator()(State const& y, State& dydr, double r);

   private:
    // Helper function: computes the derivative of N3 (number per solid angle) with respect to r.
    inline double dN3drPerOmega(double r, double n1, double n4, double gamma3);
};

/********************************************************************************************************************
 * TYPE ALIAS: ShockPair
 * DESCRIPTION: Defines a pair of Shock objects, typically representing forward and reverse shocks.
 ********************************************************************************************************************/
using ShockPair = std::pair<Shock, Shock>;

/********************************************************************************************************************
 * FUNCTION PROTOTYPES: Shock Generation Interfaces
 * DESCRIPTION: These function templates declare interfaces to generate forward shocks (2D and 3D) and
 *              forward/reverse shock pairs.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
Shock genForwardShock(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, double eps_e,
                      double eps_B);

template <typename Jet, typename Injector>
Shock genForwardShock3D(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, double eps_e,
                        double eps_B);

template <typename Jet, typename Injector>
ShockPair genFRShocks(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, double eps_e,
                      double eps_B);

template <typename Jet, typename Injector>
ShockPair genFRShocks3D(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, double eps_e,
                        double eps_B);

/********************************************************************************************************************
 * FUNCTION PROTOTYPES: Additional Shock Equation Helpers
 * DESCRIPTION: These function templates declare interfaces for finding the maximum radius and for
 *              initializing the forward shock state.
 ********************************************************************************************************************/
template <typename ShockEqn>
double find_r_max(ShockEqn& eqn, double r_min, double t_max);

template <typename ShockEqn>
void setForwardInit(ShockEqn& eqn, typename ShockEqn::State& state, double r0);

/********************************************************************************************************************
 * INLINE FUNCTIONS: Shock Utilities
 * DESCRIPTION: This section defines a set of inline functions used in shock calculations. These functions compute
 *              various physical quantities such as the comoving magnetic field (via the Weibel instability),
 *              thermal energy density, time derivatives, jet width derivative, downstream number density, fluid
 *              velocities, and update the shock state.
 ********************************************************************************************************************/

// Computes the comoving magnetic field due to the Weibel instability,
// given the magnetic energy fraction (eps_B) and the thermal energy density (e_thermal).
inline double coMovingWeibelB(double eps_B, double e_thermal) { return sqrt(8 * con::pi * eps_B * e_thermal); }

// Computes the thermal energy density in the downstream region based on the relative Lorentz factor (gamma_rel)
// and the downstream number density (n_down_str).
inline double e_ThermalDownStr(double gamma_rel, double n_down_str) {
    return n_down_str * (gamma_rel - 1) * con::mp * con::c2;
}

// Computes the derivative of the engine time with respect to radius based on the fluid speed (beta).
inline double dtdr_Engine(double beta) { return std::fabs(1 - beta) / (beta * con::c); }

// Computes the derivative of the comoving time with respect to radius,
// using the Lorentz factor (Gamma) and speed (beta).  // co-moving time
inline double dtdr_CoMoving(double Gamma, double beta) { return 1 / (Gamma * beta * con::c); };

// Computes the derivative of the jet width with respect to radius.
// Note: Does not apply to the initial non-relativistic jet.
inline double dDdr_Jet(double Gamma, double beta) {
    double constexpr cs = 0.5773502691896258 * con::c;  // sound speed approximation factor
    return cs * dtdr_CoMoving(Gamma, beta) / Gamma;
}

// Calculates the downstream number density (n4) based on energy per solid angle (dEdOmega),
// initial Lorentz factor (Gamma0), radius (r), jet width in the lab frame (D_jet_lab), and magnetization (sigma).
inline double calc_n4(double dEdOmega, double Gamma0, double r, double D_jet_lab, double sigma) {
    return dEdOmega / (Gamma0 * con::mp * con::c2 * r * r * Gamma0 * D_jet_lab) / (1 + sigma);
}

// Calculates the upstream magnetic pressure (pB4) from the downstream number density (n4) and magnetization (sigma).
inline double calc_pB4(double n4, double sigma) { return sigma * n4 * con::mp * con::c2 / 2; }

// Computes the downstream fluid velocity (u) for a given relative Lorentz factor (gamma_rel) and magnetization (sigma).
// This function uses the adiabatic index computed from gamma_rel and applies different formulas based on sigma.
inline double u_DownStr(double gamma_rel, double sigma) {
    double ad_idx = adiabaticIndex(gamma_rel);
    double gamma_m_1 = gamma_rel - 1;  // (gamma_rel - 1)
    double ad_idx_m_2 = ad_idx - 2;    // (ad_idx - 2)
    double ad_idx_m_1 = ad_idx - 1;    // (ad_idx - 1)
    if (sigma == 0) {
        return std::sqrt(gamma_m_1 * ad_idx_m_1 * ad_idx_m_1 / (-ad_idx * ad_idx_m_2 * gamma_m_1 + 2));
    } else {
        double gamma_sq = gamma_rel * gamma_rel;  // gamma_rel^2
        double gamma_p_1 = gamma_rel + 1;         // (gamma_rel + 1)
        double ad_idx_sq = ad_idx * ad_idx;       // ad_idx^2

        // Precompute common terms
        double term1 = -ad_idx * ad_idx_m_2;
        double term2 = gamma_sq - 1;
        double term3 = gamma_sq - 2;
        double term4 = gamma_p_1 * gamma_m_1;

        // Compute coefficients
        double A = term1 * gamma_m_1 + 2;
        double B = -gamma_p_1 * (-ad_idx_m_2 * (ad_idx * gamma_sq + 1) + ad_idx * ad_idx_m_1 * gamma_rel) * sigma -
                   gamma_m_1 * (term1 * term3 + 2 * gamma_rel + 3);
        double C = gamma_p_1 * (ad_idx * (1 - ad_idx / 4) * term2 + 1) * sigma * sigma +
                   term2 * (2 * gamma_rel + ad_idx_m_2 * (ad_idx * gamma_rel - 1)) * sigma +
                   term4 * gamma_m_1 * ad_idx_m_1 * ad_idx_m_1;
        double D = -gamma_m_1 * gamma_p_1 * gamma_p_1 * ad_idx_m_2 * ad_idx_m_2 * sigma * sigma / 4;

        double b = B / A;
        double c = C / A;
        double d = D / A;
        double P = c - b * b / 3;
        double Q = 2 * b * b * b / 27 - b * c / 3 + d;
        double u = std::sqrt(-P / 3);
        double uds = 2 * u * std::cos(std::acos((3 * Q / (2 * P * u)) - 2 * con::pi) / 3) - b / 3;
        return std::sqrt(uds);
    }
}

// Computes the upstream fluid velocity given the downstream velocity (u_down) and the relative Lorentz factor
// (gamma_rel).
inline double u_UpStr(double u_down, double gamma_rel) {
    return std::sqrt((1 + u_down * u_down) * (gamma_rel * gamma_rel - 1)) + u_down * gamma_rel;
}

// Computes the ratio of the upstream to downstream fluid velocities given gamma_rel and sigma.
// If the downstream velocity is zero, a fallback ratio is provided.
inline double u_UpStr2u_DownStr(double gamma_rel, double sigma) {
    double u_down_s_ = u_DownStr(gamma_rel, sigma);
    double u_up_s_ = u_UpStr(u_down_s_, gamma_rel);
    double ratio_u = u_up_s_ / u_down_s_;
    if (u_down_s_ == 0) {
        ratio_u = (7 * gamma_rel + 1) / (gamma_rel + 1);  // (g_hat+1)/(g_hat-1)
    }
    return ratio_u;
}

// Computes the downstream number density from the upstream number density (n_up_str),
// scaled by the ratio of upstream to downstream fluid velocities.
inline double n_DownStr(double n_up_str, double gamma_rel, double sigma) {
    return n_up_str * u_UpStr2u_DownStr(gamma_rel, sigma);
}

// Updates the shock state at grid cell (i, j, k) based on the provided parameters.
// If Gamma_rel > 1, the shock state is updated using computed quantities; otherwise, it is set to zero.
inline void updateShockState(Shock& shock, size_t i, size_t j, size_t k, double r, double Gamma_rel, double t_com,
                             double t_eng, double dMdOmega_up, double n_up_str, double sigma) {
    if (Gamma_rel > 1) {
        double ratio_u = u_UpStr2u_DownStr(Gamma_rel, sigma);
        double pB_up = calc_pB4(n_up_str, sigma);
        double pB_down = pB_up * ratio_u * ratio_u;
        double n_down_str = n_up_str * ratio_u;
        double co_moving_width = dMdOmega_up / (r * r * n_down_str * con::mp);
        double e_th = e_ThermalDownStr(Gamma_rel, n_down_str);
        shock.Gamma_rel[i][j][k] = Gamma_rel;
        shock.t_com[i][j][k] = t_com;
        shock.t_eng[i][j][k] = t_eng;
        shock.column_num_den[i][j][k] = n_down_str * co_moving_width;
        shock.B[i][j][k] = coMovingWeibelB(shock.eps_B, e_th) + std::sqrt(pB_down * 8 * con::pi);
    } else {
        shock.Gamma_rel[i][j][k] = 1;
        shock.t_com[i][j][k] = 0;
        shock.t_eng[i][j][k] = 0;
        shock.column_num_den[i][j][k] = 0;
        shock.B[i][j][k] = 0;
    }
}

/********************************************************************************************************************
 * FUNCTION: genForwardShock
 * DESCRIPTION: Generates a forward shock (2D) using the provided coordinates, medium, jet, injector, and energy
 *fractions. It creates a Shock object for a single phi value and iterates over theta values, solving the shock
 *              evolution for each theta slice.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
Shock genForwardShock(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, double eps_e,
                      double eps_B) {
    auto [phi_size, theta_size, r_size] = coord.shape();  // Unpack coordinate dimensions
    Shock f_shock(1, theta_size, r_size, eps_e, eps_B);   // Create Shock with 1 phi slice

    for (size_t j = 0; j < theta_size; ++j) {
        // Create a ForwardShockEqn for each theta slice (phi is set to 0)
        auto eqn = ForwardShockEqn(medium, jet, inject, 0, coord.theta[j], eps_e);
        // Solve the shock shell for this theta slice
        solveForwardShell(0, j, coord.r, f_shock, eqn, coord.t_max);
    }

    return f_shock;
}

/********************************************************************************************************************
 * FUNCTION: genForwardShock3D
 * DESCRIPTION: Generates a forward shock (3D) using the provided coordinates, medium, jet, injector, and energy
 *fractions. It creates a Shock object covering all phi and theta slices and iterates over both dimensions.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
Shock genForwardShock3D(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, double eps_e,
                        double eps_B) {
    auto [phi_size, theta_size, r_size] = coord.shape();

    Shock f_shock(phi_size, theta_size, r_size, eps_e, eps_B);  // Create Shock with full 3D dimensions
    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            // Create a ForwardShockEqn for each (phi, theta) pair
            auto eqn = ForwardShockEqn(medium, jet, inject, coord.phi[i], coord.theta[j], eps_e);
            // Solve the shock shell for this (phi, theta) slice
            solveForwardShell(i, j, coord.r, f_shock, eqn, coord.t_max);
        }
    }
    return f_shock;
}

/********************************************************************************************************************
 * FUNCTION: genFRShocks
 * DESCRIPTION: Generates a pair of forward and reverse shocks (2D) using the provided coordinates, medium, jet,
 *              injector, and energy fractions. It creates two Shock objects (one forward, one reverse) and solves
 *              the shock shells for each theta slice.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
ShockPair genFRShocks(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, double eps_e,
                      double eps_B) {
    auto [phi_size, theta_size, r_size] = coord.shape();

    Shock f_shock(1, theta_size, r_size, eps_e, eps_B);  // Forward shock for 1 phi slice
    Shock r_shock(1, theta_size, r_size, eps_e, eps_B);  // Reverse shock for 1 phi slice

    for (size_t j = 0; j < theta_size; ++j) {
        // Create equations for forward and reverse shocks for each theta slice (phi is 0)
        auto eqn_f = ForwardShockEqn(medium, jet, inject, 0, coord.theta[j], eps_e);
        auto eqn_r = FRShockEqn(medium, jet, inject, 0, coord.theta[j]);
        // Solve the forward-reverse shock shell
        solveFRShell(0, j, coord.r, f_shock, r_shock, eqn_f, eqn_r, coord.t_max);
    }

    return std::make_pair(std::move(f_shock), std::move(r_shock));
}

/********************************************************************************************************************
 * FUNCTION: genFRShocks3D
 * DESCRIPTION: Generates a pair of forward and reverse shocks (3D) using the provided coordinates, medium, jet,
 *              injector, and energy fractions. It creates two Shock objects covering all phi and theta slices and
 *              solves the shock shells for each (phi, theta) pair.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
ShockPair genFRShocks3D(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, double eps_e,
                        double eps_B) {
    auto [phi_size, theta_size, r_size] = coord.shape();

    Shock f_shock(phi_size, theta_size, r_size, eps_e, eps_B);  // Forward shock for full 3D dimensions
    Shock r_shock(phi_size, theta_size, r_size, eps_e, eps_B);  // Reverse shock for full 3D dimensions
    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            // Create equations for forward and reverse shocks for each (phi, theta) pair
            auto eqn_f = ForwardShockEqn(medium, jet, inject, coord.phi[i], coord.theta[j], eps_e);
            auto eqn_r = FRShockEqn(medium, jet, inject, coord.phi[i], coord.theta[j]);
            // Solve the forward-reverse shock shell for this (phi, theta) pair
            solveFRShell(i, j, coord.r, f_shock, r_shock, eqn_f, eqn_r, coord.t_max);
        }
    }
    return std::make_pair(std::move(f_shock), std::move(r_shock));
}

/********************************************************************************************************************
 * FUNCTION: find_r_max
 * DESCRIPTION: Finds the maximum radius (r_max) for the shock evolution by integrating the shock equations until
 *              the engine time exceeds t_max. Uses a dense output stepper.
 ********************************************************************************************************************/
template <typename ShockEqn>
double find_r_max(ShockEqn& eqn, double r_min, double t_max) {
    using namespace boost::numeric::odeint;
    double atol = 0, rtol = 1e-6, r0 = r_min;
    double dr = r_min / 100;

    typename ShockEqn::State state;
    setForwardInit(eqn, state, r0);  // Initialize the state at r0

    auto stepper = make_dense_output(atol, rtol, runge_kutta_dopri5<typename ShockEqn::State>());
    stepper.initialize(state, r0, dr);
    // Integrate until the engine time in the state exceeds t_max
    for (; state[2] <= t_max;) {
        stepper.do_step(eqn);
        state = stepper.current_state();
    }
    return stepper.current_time() + stepper.current_time_step();
}

/********************************************************************************************************************
 * FUNCTION: setReverseInit
 * DESCRIPTION: Initializes the state vector for reverse shock evolution at the given radius r0.
 ********************************************************************************************************************/
template <typename ShockEqn>
void setReverseInit(ShockEqn& eqn, typename ShockEqn::State& state, double r0) {
    double gamma4 = eqn.jet.Gamma0(eqn.phi, eqn.theta, 0);  // Obtain initial Gamma from the jet
    double u0 = (gamma4 - 1) * eqn.medium.mass(r0) / (4 * con::pi) * con::c2;
    double beta0 = gammaTobeta(gamma4);
    double t_eng0 = r0 * (1 - beta0) / beta0 / con::c;
    double t_com0 = r0 / std::sqrt(gamma4 * gamma4 - 1) / con::c;
    double D_jet0 = con::c * eqn.jet.duration;
    double dN3dOmega = 0;  // Initialize number per unit solid angle to zero
    state = {0., dN3dOmega, t_eng0, t_com0, D_jet0};
}

/********************************************************************************************************************
 * TEMPLATE HELP FUNCTIONS
 * DESCRIPTION: The following template functions are helper functions for the shock generation and evolution process.
 ********************************************************************************************************************/

// Determines the range of radii (r_min, r_max) for the shock evolution using medium and jet properties.
// The minimum radius is based on t_min and the maximum is found by solving the forward shock equation on-axis.
template <typename Jet, typename Injector>
std::tuple<double, double> findRadiusRange(Medium const& medium, Jet const& jet, Injector const& inj, double t_min,
                                           double t_max, double z = 0) {
    double gamma0 = jet.Gamma0(0, 0, 0);  // On-axis Gamma
    double beta0 = std::sqrt(1 - 1 / (gamma0 * gamma0));
    double gamma_min = (gamma0 - 1) / 100 + 1;
    double beta_min = std::sqrt(1 - 1 / (gamma_min * gamma_min));
    double r_min = beta_min * con::c * t_min / (1 + beta_min);

    // Find the on-axis r_max by solving for theta = 0, phi = 0
    auto eqn = ForwardShockEqn<Jet, Injector>(medium, jet, inj, 0, 0, 0);
    double r_max = find_r_max(eqn, r_min, t_max);
    return {r_min, r_max};
}

// Updates the forward shock state at grid index (i, j, k) using the current state vector and the medium properties.
template <typename ShockEqn>
void updateForwardShock(size_t i, size_t j, int k, double r_k, ShockEqn& eqn, const typename ShockEqn::State& state,
                        Shock& f_shock) {
    double n1 = eqn.medium.rho(r_k) / con::mp;                // Compute upstream number density
    double Gamma = state[0];                                  // Lorentz factor from state
    double t_eng = state[2];                                  // Engine time from state
    double t_com = state[3];                                  // Comoving time from state
    double dM1dOmega = eqn.medium.mass(r_k) / (4 * con::pi);  // Mass per unit solid angle

    updateShockState(f_shock, i, j, k, r_k, Gamma, t_com, t_eng, dM1dOmega, n1, eqn.jet_sigma);
}

// Initializes the forward shock state vector at radius r0.
template <typename ShockEqn>
void setForwardInit(ShockEqn& eqn, typename ShockEqn::State& state, double r0) {
    double gamma2 = eqn.jet.Gamma0(eqn.phi, eqn.theta, 0);  // Initial Lorentz factor
    double u0 = (gamma2 - 1) * eqn.medium.mass(r0) / (4 * con::pi) * con::c2;
    double beta0 = gammaTobeta(gamma2);
    double t_eng0 = r0 * (1 - beta0) / beta0 / con::c;
    double t_com0 = r0 / std::sqrt(gamma2 * gamma2 - 1) / con::c;
    double D_jet0 = con::c * eqn.jet.duration;
    state = {gamma2, u0, t_eng0, t_com0, D_jet0};
}

// Determines whether a reverse shock exists based on current shock parameters.
template <typename ShockEqn>
bool reverseShockExists(ShockEqn const& eqn, double r, double gamma, double t_eng, double D_jet) {
    double n4 = calc_n4(eqn.jet.dEdOmega(eqn.phi, eqn.theta, t_eng), gamma, r, D_jet, eqn.jet_sigma);
    double n1 = eqn.medium.rho(r) / con::mp;
    return eqn.jet_sigma < 8. / 3 * gamma * gamma * n1 / n4;
}

// Solves the forward shock evolution for a given shell (across radius values in array r) and updates the Shock object.
template <typename ShockEqn>
void solveForwardShell(size_t i, size_t j, const Array& r, Shock& f_shock, ShockEqn& eqn, double t_max) {
    using namespace boost::numeric::odeint;

    double atol = 0, rtol = 1e-6, r0 = r[0];
    double dr = (r[1] - r[0]) / 100;

    typename ShockEqn::State state;
    setForwardInit(eqn, state, r0);  // Initialize state at starting radius

    if (state[0] <= con::Gamma_cut) {  // If initial Lorentz factor is too low, exit early
        return;
    }

    // Create a dense output stepper for integrating the shock equations
    auto stepper = make_dense_output(atol, rtol, runge_kutta_dopri5<typename ShockEqn::State>());
    stepper.initialize(state, r0, dr);

    double r_back = r[r.size() - 1];  // Last radius in the array

    // Iterate over the radius array, updating the state and the Shock object as needed
    for (int k = 0; stepper.current_time() <= r_back && state[2] <= t_max;) {
        stepper.do_step(eqn);
        while (k < r.size() && stepper.current_time() > r[k]) {
            stepper.calc_state(r[k], state);
            updateForwardShock(i, j, k, r[k], eqn, state, f_shock);
            ++k;
        }
    }
}
//----------------------------------- definitions of ForwardShockEqn -----------------------------------//
/********************************************************************************************************************
 * METHOD: ForwardShockEqn::operator()(State const& y, State& dydr, double r)
 * DESCRIPTION: Computes the derivatives of the state variables with respect to radius r.
 *              The state vector components are:
 *                  y[0] - Gamma (Lorentz factor)
 *                  y[1] - u (internal energy-related variable)
 *                  y[2] - t_eng (engine time)
 *                  y[3] - t_com (co-moving time) [unused here]
 *                  y[4] - D_jet (jet shell width) [unused here]
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
void ForwardShockEqn<Jet, Injector>::operator()(State const& y, State& dydr, double r) {
    double Gamma = y[0];
    double u = y[1];
    double t_eng = y[2];  // engine time
    // double t_com = y[3];  // co-moving time (unused)
    // double D_jet = y[4];  // co-moving jet shell width (unused)

    double ad_idx = adiabaticIndex(Gamma);  // Compute adiabatic index based on Gamma
    double rho = medium.rho(r);             // Get medium density at radius r
    double beta = gammaTobeta(Gamma);       // Convert Gamma to beta (velocity/c)
    double beta4 = gammaTobeta(gamma4);     // Convert gamma4 to beta

    dydr[2] = dtdr_Engine(beta);  // Compute derivative of engine time with respect to r
    dydr[0] = dGammadr(r, Gamma, u, t_eng, ad_idx, rho, dydr[2]);  // d(Gamma)/dr
    dydr[1] = dUdr(r, Gamma, u, t_eng, ad_idx, rho, dydr[0]);      // d(u)/dr
    dydr[3] = dtdr_CoMoving(Gamma, beta);                          // d(t_com)/dr
    dydr[4] = dDdr_Jet(gamma4, beta4);                             // d(D_jet)/dr
}

/********************************************************************************************************************
 * CONSTRUCTOR: ForwardShockEqn::ForwardShockEqn
 * DESCRIPTION: Initializes a ForwardShockEqn object with references to the medium, jet, and injector,
 *              along with the angular coordinates and energy fraction.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
ForwardShockEqn<Jet, Injector>::ForwardShockEqn(Medium const& medium, Jet const& jet, Injector const& inject,
                                                double phi, double theta, double eps_e)
    : medium(medium),
      jet(jet),
      inject(inject),
      phi(phi),
      theta(theta),
      eps_e(eps_e),
      jet_sigma(jet.sigma0(phi, theta, 0)),
      inj_sigma(inject.sigma0(phi, theta, 0)),
      spreading_factor(1),
      jet_Gamma0(jet.Gamma0(phi, theta, 0)),
      inj_Gamma0(inject.Gamma0(phi, theta, 0)),
      gamma4(jet_Gamma0),
      dM0(jet.dEdOmega(phi, theta, 0) / (jet_Gamma0 * (1 + jet_sigma) * con::c2)) {
    // dM0dOmega(jet.dE0dOmega(theta) / (jet.Gamma0(theta) * con::c2)) is commented out.
}

/********************************************************************************************************************
 * METHOD: ForwardShockEqn::dGammadr
 * DESCRIPTION: Computes the derivative of Gamma with respect to radius r.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
double ForwardShockEqn<Jet, Injector>::dGammadr(double r, double Gamma, double u, double t_eng, double ad_idx,
                                                double rho, double dtdr) {
    double ad_idx_m1 = ad_idx - 1;
    double Gamma2_m1 = Gamma * Gamma - 1;
    double term1 = ad_idx * Gamma2_m1 + 1;

    double dm = medium.mass(r) / (4 * con::pi);  // Mass per unit solid angle from medium
    double dm_inj = inject.dEdOmega(phi, theta, t_eng) / (inj_Gamma0 * (1 + inj_sigma) * con::c2);  // Injected mass
    double L_inj = inject.dLdOmega(phi, theta, t_eng);  // Injected luminosity per unit solid angle

    double a1 = -Gamma2_m1 * (ad_idx * Gamma - ad_idx + 1) * r * r * rho * con::c2;
    double a2 = ad_idx_m1 * term1 * 3 * u / r;
    double a3 = Gamma * dtdr * (L_inj * (1 - Gamma / (inj_Gamma0 * (1 + inj_sigma))));

    double b1 = Gamma * (dM0 + dm + dm_inj) * con::c2;
    double b2 = (ad_idx * term1 + 2 * ad_idx_m1) / Gamma * u;

    return (a1 + a2 + a3) / (b1 + b2);
}

/********************************************************************************************************************
 * METHOD: ForwardShockEqn::dUdr
 * DESCRIPTION: Computes the derivative of u with respect to radius r.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
double ForwardShockEqn<Jet, Injector>::dUdr(double r, double Gamma, double u, double t_eng, double ad_idx, double rho,
                                            double dGdr) {
    double E = r * r * rho * con::c2;
    return (1 - eps_e) * (Gamma - 1) * E - (ad_idx - 1) * (3 / r - dGdr / Gamma) * u * spreading_factor;
}

//----------------------------------- definition of RShockEqn -----------------------------------//
/********************************************************************************************************************
 * INLINE FUNCTION: calc_gamma3
 * DESCRIPTION: Computes gamma3 for the reverse shock based on radius, upstream and downstream densities,
 *              the jet Lorentz factor (gamma4), and magnetization (sigma).
 ********************************************************************************************************************/
inline double calc_gamma3(double r, double n1, double n4, double gamma4, double sigma) {
    double C = n4 / n1 * (1 + sigma);
    double gamma3 =
        gamma4 * std::sqrt(((C - 2) - 2 * std::sqrt(C * (gamma4 * gamma4 - 1) + 1)) / (C - 4 * gamma4 * gamma4));
    return gamma3;
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::dN3drPerOmega
 * DESCRIPTION: Computes the derivative of N3 (number per unit solid angle) with respect to radius.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
double FRShockEqn<Jet, Injector>::dN3drPerOmega(double r, double n1, double n4, double gamma3) {
    double gamma34 = (gamma4 / gamma3 + gamma3 / gamma4) / 2;
    double ratio_u = u_UpStr2u_DownStr(gamma34, this->jet_sigma);

    /* The following commented-out section shows an alternative computation.
    double ad_idx2 = adiabatic_index(Gamma);
    double ad_idx3 = adiabatic_index(Gamma34);
    double u3s_ = u_down_str(Gamma34, this->sigma);
    double u4s_ = u_up_str(u3s_, Gamma34);
    double n2 = n_down_str(n1, Gamma, this->sigma);
    double e2 = e_thermal_down_str(Gamma, n2);
    double p2 = (ad_idx2 - 1) * e2;
    double pB4 = calc_pB4(n4, this->sigma);
    double pB3 = pB4 * ratio_u * ratio_u;
    double f_a = fa(Gamma34, u3s_, this->sigma);
    double f_b = ratio_u / ((ad_idx3 * Gamma34 + 1) / (ad_idx3 - 1));
    double f_c = fc(p2, pB3);
    double F = f_a * f_b * f_c;
    */
    double n3 = n4 * ratio_u;
    double dxdr = 1. / (gamma4 * std::sqrt((1 + this->jet_sigma) * n4 / n1) * std::fabs(1 - gamma4 * n4 / gamma3 / n3));
    return n3 * r * r * gamma3 * dxdr;
}

/********************************************************************************************************************
 * CONSTRUCTOR: FRShockEqn::FRShockEqn
 * DESCRIPTION: Initializes an FRShockEqn object with references to the medium, jet, and injector, and sets the
 *              angular coordinates, jet magnetization, and initial Gamma.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
FRShockEqn<Jet, Injector>::FRShockEqn(Medium const& medium, Jet const& jet, Injector const& inject, double phi,
                                      double theta)
    : medium(medium),
      jet(jet),
      inject(inject),
      phi(phi),
      theta(theta),
      jet_sigma(jet.sigma0(phi, theta, 0)),
      gamma4(jet.Gamma0(phi, theta, 0)) {}

/********************************************************************************************************************
 * METHOD: FRShockEqn::operator()(State const& y, State& dydr, double r)
 * DESCRIPTION: Computes the derivatives for the reverse shock evolution.
 *              The state vector for FRShockEqn is similar to that of ForwardShockEqn.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
void FRShockEqn<Jet, Injector>::operator()(State const& y, State& dydr, double r) {
    // y[0] is left blank.
    // double D_rs = y[1];
    double t_eng = y[2];
    // double t_com = y[3];
    double D_jet_lab = y[4];

    double n4 = calc_n4(jet.dEdOmega(phi, theta, t_eng), gamma4, r, D_jet_lab, this->jet_sigma);
    double n1 = medium.rho(r) / con::mp;

    double gamma3 = calc_gamma3(r, n1, n4, this->gamma4, this->jet_sigma);
    double beta3 = gammaTobeta(gamma3);
    double beta4 = gammaTobeta(this->gamma4);
    dydr[0] = 0;
    dydr[1] = dN3drPerOmega(r, n1, n4, gamma3);
    dydr[2] = dtdr_Engine(beta3);
    dydr[3] = dtdr_CoMoving(gamma3, beta3);
    dydr[4] = dDdr_Jet(this->gamma4, beta4);
}

/********************************************************************************************************************
 * STRUCT: CrossState
 * DESCRIPTION: Represents the state variables at the crossing between forward and reverse shock phases.
 ********************************************************************************************************************/
struct CrossState {
    double gamma_rel;
    double r;
    double column_num_den;
    double B;
};

/********************************************************************************************************************
 * INLINE FUNCTION: Blandford_McKee
 * DESCRIPTION: Updates the shock state at a grid cell using the Blandford–McKee self-similar solution.
 ********************************************************************************************************************/
inline void Blandford_McKee(size_t i, size_t j, size_t k, Shock& shock, CrossState const& state_c, double r,
                            double t_com, double t_eng) {
    double const g = 2.;
    shock.t_com[i][j][k] = t_com;
    shock.t_eng[i][j][k] = t_eng;
    // The following lines are alternative formulations (commented out):
    // shock.Gamma[j][k] = (state_c.gamma - 1) * std::pow(r / state_c.r, -g) + 1;
    // shock.n_p[j][k] = state_c.n3 * std::pow(r / state_c.r, -6 * (3 + g) / 7);
    // shock.e_th[j][k] = state_c.e3 * std::pow(r / state_c.r, -8 * (3 + g) / 7);
    // shock.width_eff[j][k] = state_c.D_eff * std::pow(r / state_c.r, (6. * (3. + g) - 14.) / 7.);
    shock.Gamma_rel[i][j][k] = (state_c.gamma_rel - 1) * std::pow(r / state_c.r, -g) + 1;
    shock.column_num_den[i][j][k] = state_c.column_num_den * std::pow(r / state_c.r, -2);
    shock.B[i][j][k] = state_c.B * std::pow(r / state_c.r, (6. * (3. + g) - 14.) / 7.);
}

/********************************************************************************************************************
 * FUNCTION: solveFRShell
 * DESCRIPTION: Solves the evolution of the forward–reverse shock shell over the radius array r, updating both
 *              the forward shock (f_shock) and reverse shock (r_shock) objects accordingly.
 ********************************************************************************************************************/
template <typename FShockEqn, typename RShockEqn>
void solveFRShell(size_t i, size_t j, Array const& r, Shock& f_shock, Shock& r_shock, FShockEqn& eqn_f,
                  RShockEqn& eqn_r, double t_max) {
    using namespace boost::numeric::odeint;
    double atol = 0, rtol = 1e-6, r0 = r[0];
    double dr = (r[1] - r[0]) / 100;
    typename FShockEqn::State state;

    auto stepper = bulirsch_stoer_dense_out<typename FShockEqn::State>{atol, rtol};
    // Alternatively:
    // auto stepper = make_dense_output(atol, rtol, runge_kutta_dopri5<std::vector<double>>());

    setForwardInit(eqn_f, state, r0);
    if (state[0] <= 1 + 1e-6) {  // If the initial Lorentz factor is too low, exit early
        return;
    }
    double gamma3 = eqn_f.gamma4;
    double t_eng = state[2];
    double t_com = state[3];
    double D_jet = state[4];

    bool RS_crossing = reverseShockExists(eqn_r, r0, gamma3, t_eng, D_jet);
    bool crossed = false;

    if (RS_crossing) {
        setReverseInit(eqn_r, state, r0);
    }

    double t_com_last = state[3];
    double n1 = 0, n4 = 0;
    double dN3dOmega = 0;
    double dN4dOmega =
        eqn_r.jet.dEdOmega(eqn_r.phi, eqn_r.theta, t_eng) / (eqn_r.gamma4 * con::mp * con::c2 * (1 + eqn_r.jet_sigma));

    CrossState state_c;
    stepper.initialize(state, r0, dr);
    // Integrate the shell over the radius array r
    double r_back = r[r.size() - 1];
    for (int k = 0; stepper.current_time() <= r_back && state[2] <= t_max;) {
        RS_crossing ? stepper.do_step(eqn_r) : stepper.do_step(eqn_f);

        for (; stepper.current_time() > r[k] && k < r.size(); k++) {
            n1 = eqn_f.medium.rho(r[k]) / con::mp;
            stepper.calc_state(r[k], state);

            t_eng = state[2];
            t_com = state[3];
            D_jet = state[4];
            if (RS_crossing) {
                n4 = calc_n4(eqn_r.jet.dEdOmega(eqn_r.phi, eqn_r.theta, t_eng), eqn_r.gamma4, r[k], D_jet,
                             eqn_r.jet_sigma);
                gamma3 = calc_gamma3(r[k], n1, n4, eqn_r.gamma4, eqn_r.jet_sigma);
            } else {
                gamma3 = state[0];
            }

            double dM2dOmega = eqn_f.medium.mass(r[k]) / (4 * con::pi);
            updateShockState(f_shock, i, j, k, r[k], gamma3, t_com, t_eng, dM2dOmega, n1, eqn_f.jet_sigma);

            if (!crossed && RS_crossing) {
                dN3dOmega = state[1];
                double gamma34 = (eqn_r.gamma4 / gamma3 + gamma3 / eqn_r.gamma4) / 2;
                double dM3dOmega = dN3dOmega * con::mp * con::c2;
                updateShockState(r_shock, i, j, k, r[k], gamma34, t_com, t_eng, dM3dOmega, n4, eqn_r.jet_sigma);

                crossed = dN3dOmega >= dN4dOmega;
                if (crossed) {
                    RS_crossing = false;
                    state_c = {r_shock.Gamma_rel[i][j][k], r[k], r_shock.column_num_den[i][j][k], r_shock.B[i][j][k]};
                    double u0 = (gamma3 - 1) * eqn_r.medium.mass(r[k]) / (4 * con::pi) * con::c2;
                    state = {gamma3, u0, t_eng, t_com, D_jet};
                    eqn_f.gamma4 = gamma3;
                    stepper.initialize(state, r[k], dr);
                }
            } else if (!crossed && !RS_crossing) {
                RS_crossing = reverseShockExists(eqn_r, r0, gamma3, t_eng, D_jet);
                if (RS_crossing) {
                    state = {0., 0., t_eng, t_com, D_jet};
                    stepper.initialize(state, r[k], dr);
                }
            } else {
                Blandford_McKee(i, j, k, r_shock, state_c, r[k], t_com, t_eng);
            }
        }
    }
}
#endif