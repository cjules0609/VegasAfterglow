//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <tuple>

#include "boost/numeric/odeint.hpp"
#include "jet.h"
#include "medium.h"
#include "mesh.h"
#include "physics.h"

/********************************************************************************************************************
 * @class Shock
 * @brief Represents a shock wave in an astrophysical environment.
 * @details The class stores physical properties of the shock across a 3D grid defined by azimuthal angle (phi),
 *          polar angle (theta), and time bins. Provides methods for shock calculations, including relativistic
 *          jump conditions, magnetic field calculations, and energy density computations.
 ********************************************************************************************************************/
class Shock {
   public:
    // Constructor: Initialize with grid dimensions and energy fractions
    Shock(size_t phi_size, size_t theta_size, size_t t_size, Real eps_e, Real eps_B);
    Shock() noexcept = default;
    Shock(Shock const& other) noexcept = default;
    Shock(Shock&& other) noexcept = default;
    Shock& operator=(Shock&& other) noexcept = default;
    Shock& operator=(Shock const& other) noexcept = default;

    MeshGrid3d t_comv;          // Comoving time
    MeshGrid3d r;               // Radius
    MeshGrid3d theta;           // Theta for jet spreading
    MeshGrid3d Gamma;           // Bulk Lorentz factor
    MeshGrid3d Gamma_rel;       // Relative Lorentz factor between downstream and upstream
    MeshGrid3d B;               // Comoving magnetic field
    MeshGrid3d column_num_den;  // Downstream proton column number density
    MeshGrid injection_idx;     // Beyond which grid index there is no electron injection
    MaskGrid required;          // Grid points actually required for final flux calculation
    Real eps_e{0};              // Electron energy fraction
    Real eps_B{0};              // Magnetic energy fraction

    // Returns grid dimensions as a tuple
    auto shape() const { return std::make_tuple(phi_size, theta_size, t_size); }

    // Resizes all grids to the given dimensions
    void resize(size_t phi_size, size_t theta_size, size_t t_size);

   private:
    size_t phi_size{0};    // Number of grid points in phi direction
    size_t theta_size{0};  // Number of grid points in theta direction
    size_t t_size{0};      // Number of grid points in time direction
};

/********************************************************************************************************************
 * @defgroup ShockUtilities Shock Utilities
 * @brief Inline functions used in shock calculations.
 * @details This section defines a set of inline functions used in shock calculations. These functions compute
 *          various physical quantities such as the comoving magnetic field (via the Weibel instability),
 *          thermal energy density, time derivatives, jet width derivative, downstream number density, fluid
 *          velocities, and update the shock state.
 ********************************************************************************************************************/

/********************************************************************************************************************
 * @brief Computes the downstream four-velocity for a given relative Lorentz factor and magnetization parameter.
 * @param gamma_rel Relative Lorentz factor between upstream and downstream regions
 * @param sigma Magnetization parameter (ratio of magnetic to rest-mass energy density)
 * @return The downstream four-velocity in the shock frame
 * @details The calculation handles both magnetized (sigma > 0) and non-magnetized (sigma = 0) cases using
 *          different equations based on jump conditions across the shock front.
 ********************************************************************************************************************/
Real compute_downstr_4vel(Real gamma_rel, Real sigma);

/********************************************************************************************************************
 * @brief Computes the ratio of upstream to downstream four-velocity across the shock front.
 * @param gamma_rel Relative Lorentz factor between upstream and downstream regions
 * @param sigma Magnetization parameter (ratio of magnetic to rest-mass energy density)
 * @return The ratio of upstream to downstream four-velocity
 * @details This ratio is a key parameter in determining various shock properties, such as compression ratio
 *          and jump conditions for density, pressure, and magnetic field.
 ********************************************************************************************************************/
Real compute_4vel_jump(Real gamma_rel, Real sigma);

/********************************************************************************************************************
 * @brief Computes the sound speed in the shocked medium based on the relative Lorentz factor.
 * @param Gamma_rel Relative Lorentz factor between upstream and downstream regions
 * @return The sound speed as a fraction of light speed
 ********************************************************************************************************************/
inline Real compute_sound_speed(Real Gamma_rel) {
    Real ad_idx = adiabatic_idx(Gamma_rel);
    return std::sqrt(ad_idx * (ad_idx - 1) * (Gamma_rel - 1) / (1 + (Gamma_rel - 1) * ad_idx));
}

/********************************************************************************************************************
 * @brief Computes the effective Lorentz factor accounting for the adiabatic index.
 * @param adx Adiabatic index of the medium
 * @param Gamma Bulk Lorentz factor
 * @return The effective Lorentz factor
 ********************************************************************************************************************/
inline Real compute_effective_Gamma(Real adx, Real Gamma) { return (adx * Gamma * Gamma - adx + 1) / Gamma; }

/********************************************************************************************************************
 * @brief Computes the comoving magnetic field using the Weibel instability mechanism.
 * @param eps_B Fraction of thermal energy in magnetic fields
 * @param e_thermal Thermal energy density
 * @return The comoving magnetic field strength
 ********************************************************************************************************************/
inline Real compute_comv_weibel_B(Real eps_B, Real e_thermal) { return std::sqrt(8 * con::pi * eps_B * e_thermal); }

/********************************************************************************************************************
 * @brief Computes the time derivative of radius (dr/dt) based on the shock velocity.
 * @param beta Shock velocity as a fraction of light speed
 * @return The rate of change of radius with respect to observer time
 ********************************************************************************************************************/
inline Real compute_dr_dt(Real beta) { return (beta * con::c) / (1 - beta); }

/********************************************************************************************************************
 * @brief Computes the time derivative of theta (dθ/dt) for jet spreading.
 * @param theta_s Typical spreading angle of the jet
 * @param theta Theta of the current grid point
 * @param drdt Time derivative of radius
 * @param r Current radius
 * @param Gamma Current bulk Lorentz factor
 * @return The rate of change of the half-opening angle
 ********************************************************************************************************************/
inline Real compute_dtheta_dt(Real theta_s, Real theta, Real drdt, Real r, Real Gamma) {
    constexpr Real Q = 2.82842712475;
    Real u2 = Gamma * Gamma - 1;
    Real ratio = u2 / (Q * Q * theta_s * theta_s);
    Real x = theta / theta_s;
    Real f = 1 / (1 + ratio) * x;
    return drdt / (2 * Gamma * r) * std::sqrt((2 * u2 + 3) / (4 * u2 + 3)) * f;
}

/********************************************************************************************************************
 * @brief Computes the time derivative of comoving time (dt_comv/dt) based on the Lorentz factor.
 * @param Gamma Bulk Lorentz factor
 * @param beta Bulk velocity as a fraction of light speed
 * @return The rate of change of comoving time with respect to observer time
 ********************************************************************************************************************/
inline Real compute_dt_dt_comv(Real Gamma, Real beta) { return 1 / (Gamma * (1 - beta)); };

/********************************************************************************************************************
 * @brief Computes the upstream magnetic pressure based on number density and magnetization.
 * @param n_up Upstream number density
 * @param sigma Magnetization parameter
 * @return The upstream magnetic pressure
 ********************************************************************************************************************/
inline Real compute_upstr_mag_p(Real n_up, Real sigma) { return sigma * n_up * con::mp * con::c2 / 2; }

/********************************************************************************************************************
 * @brief Computes the upstream four-velocity from downstream four-velocity and relative Lorentz factor.
 * @param u_down Downstream four-velocity
 * @param gamma_rel Relative Lorentz factor
 * @return The upstream four-velocity in the shock frame
 ********************************************************************************************************************/
inline Real compute_upstr_4vel(Real u_down, Real gamma_rel) {
    return std::sqrt((1 + u_down * u_down) * (gamma_rel * gamma_rel - 1)) + u_down * gamma_rel;
}

/********************************************************************************************************************
 * @brief Computes the downstream number density from upstream density, relative Lorentz factor, and magnetization.
 * @param n_up_str Upstream number density
 * @param gamma_rel Relative Lorentz factor
 * @param sigma Magnetization parameter
 * @return The downstream number density
 * @details Uses the shock jump conditions to determine the density of the downstream region.
 ********************************************************************************************************************/
inline Real compute_downstr_num_den(Real n_up_str, Real gamma_rel, Real sigma) {
    return n_up_str * compute_4vel_jump(gamma_rel, sigma);
}

/********************************************************************************************************************
 * @brief Computes the relative Lorentz factor between two frames with given Lorentz factors.
 * @param gamma1 First frame's Lorentz factor
 * @param gamma2 Second frame's Lorentz factor
 * @return The relative Lorentz factor between the two frames
 ********************************************************************************************************************/
inline Real compute_rel_Gamma(Real gamma1, Real gamma2) {
    return gamma1 * gamma2 - std::sqrt((gamma1 * gamma1 - 1) * (gamma2 * gamma2 - 1));
}

/********************************************************************************************************************
 * @brief Computes the relative Lorentz factor between two frames with given Lorentz factors and velocities.
 * @param gamma1 First frame's Lorentz factor
 * @param gamma2 Second frame's Lorentz factor
 * @param beta1 First frame's velocity as fraction of light speed
 * @param beta2 Second frame's velocity as fraction of light speed
 * @return The relative Lorentz factor between the two frames
 ********************************************************************************************************************/
inline Real compute_rel_Gamma(Real gamma1, Real gamma2, Real beta1, Real beta2) {
    return gamma1 * gamma2 * (1 - beta1 * beta2);
}

/********************************************************************************************************************
 * @brief Computes a Lorentz factor from a reference Lorentz factor and relative Lorentz factor.
 * @param gamma4 Reference Lorentz factor (typically for region 4, unshocked ejecta)
 * @param gamma_rel Relative Lorentz factor
 * @return The derived Lorentz factor
 ********************************************************************************************************************/
inline Real compute_Gamma_from_relative(Real gamma4, Real gamma_rel) {
    Real b = -2 * gamma4 * gamma_rel;
    Real c = gamma4 * gamma4 + gamma_rel * gamma_rel - 1;
    return (-b - std::sqrt(b * b - 4 * c)) / 2;
}

/********************************************************************************************************************
 * @brief Computes the downstream thermal energy density.
 * @param gamma_rel Relative Lorentz factor
 * @param n_down_str Downstream number density
 * @return The thermal energy density in the downstream region
 ********************************************************************************************************************/
inline Real compute_downstr_eth(Real gamma_rel, Real n_down_str) {
    return n_down_str * (gamma_rel - 1) * con::mp * con::c2;
}

/********************************************************************************************************************
 * @brief Computes the rate at which the shock shell spreads in the comoving frame.
 * @param Gamma_rel Relative Lorentz factor
 * @param dtdt_com Rate of change of comoving time with respect to observer time
 * @return The shell spreading rate in the comoving frame
 ********************************************************************************************************************/
inline Real compute_shell_spreading_rate(Real Gamma_rel, Real dtdt_com) {
    Real cs = compute_sound_speed(Gamma_rel);
    return cs * dtdt_com;
}
/*
inline Real dN3dt(Real r, Real n4, Real gamma3, Real gamma4, Real sigma) {
    if (gamma3 == gamma4) {
        return 0.;
    }
    Real beta3 = gammaTobeta(gamma3);
    Real beta4 = gammaTobeta(gamma4);
    Real gamma34 = relativeLorentz(gamma4, gamma3, beta4, beta3);
    Real ratio_u = u_UpStr2u_DownStr(gamma34, sigma);
    Real n3 = n4 * ratio_u;
    Real dxdt = (beta4 - beta3) * con::c / ((1 - beta3) * (1 - gamma4 / (gamma3 * ratio_u)));
    return n3 * r * r * gamma3 * dxdt;
}*/

/********************************************************************************************************************
 * @brief Computes the rate of change of shocked ejecta mass per solid angle.
 * @param width Width of the shocked region
 * @param m4 Mass per solid angle in region 4 (unshocked ejecta)
 * @param gamma3 Lorentz factor in region 3 (shocked ejecta)
 * @param gamma4 Lorentz factor in region 4 (unshocked ejecta)
 * @param sigma Magnetization parameter
 * @return The rate of change of mass in region 3
 ********************************************************************************************************************/
inline Real compute_dm3_dt(Real width, Real m4, Real gamma3, Real gamma4, Real sigma) {
    if (gamma3 == gamma4) {
        return 0.;
    }
    Real beta3 = gamma_to_beta(gamma3);
    Real beta4 = gamma_to_beta(gamma4);
    Real gamma34 = compute_rel_Gamma(gamma4, gamma3, beta4, beta3);
    Real ratio_u = compute_4vel_jump(gamma34, sigma);
    Real column_den3 = m4 * ratio_u / width;
    Real dxdt = (beta4 - beta3) * con::c / ((1 - beta3) * (1 - gamma4 / (gamma3 * ratio_u)));
    return column_den3 * gamma3 * dxdt;
}

/********************************************************************************************************************
 * @brief Computes the number density in region 4 (unshocked ejecta).
 * @param dEdOmega Energy per solid angle
 * @param Gamma0 Initial Lorentz factor
 * @param r Radius
 * @param D_jet Jet thickness
 * @param sigma Magnetization parameter
 * @return The number density in region 4
 ********************************************************************************************************************/
inline Real compute_region4_num_den(Real dEdOmega, Real Gamma0, Real r, Real D_jet, Real sigma) {
    return dEdOmega / ((Gamma0 * con::mp * con::c2 * r * r * D_jet) * (1 + sigma));
}

/********************************************************************************************************************
 * @brief Sets a stopping shock state when the Lorentz factor drops below threshold.
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param shock Reference to the Shock object to be updated
 * @param state0 Initial state to be used for some parameters
 ********************************************************************************************************************/
template <typename State>
inline void set_stopping_shock(size_t i, size_t j, Shock& shock, State const& state0) {
    xt::view(shock.t_comv, i, j, xt::all()) = state0.t_comv;
    xt::view(shock.r, i, j, xt::all()) = state0.r;
    xt::view(shock.theta, i, j, xt::all()) = state0.theta;
    xt::view(shock.Gamma, i, j, xt::all()) = 1;
    xt::view(shock.Gamma_rel, i, j, xt::all()) = 1;
    xt::view(shock.B, i, j, xt::all()) = 0;
    xt::view(shock.column_num_den, i, j, xt::all()) = 0;
}

/********************************************************************************************************************
 * @brief Saves the shock state at a given grid point (i,j,k).
 * @param shock Reference to the Shock object to store the state
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param k Grid index for time
 * @param state Current state of the system
 * @param Gamma_downstr Downstream Lorentz factor
 * @param Gamma_upstr Upstream Lorentz factor
 * @param N_downstr Downstream total number of shocked particles
 * @param n_upstr Upstream number density
 * @param sigma_upstr Upstream magnetization
 ********************************************************************************************************************/
template <typename State>
void save_shock_state(Shock& shock, size_t i, size_t j, size_t k, State const& state, Real Gamma_downstr,
                      Real Gamma_upstr, Real N_downstr, Real n_upstr, Real sigma_upstr) {
    Real Gamma_rel = compute_rel_Gamma(Gamma_upstr, Gamma_downstr);
    Real ratio_u = compute_4vel_jump(Gamma_rel, sigma_upstr);
    Real pB_upstr = compute_upstr_mag_p(n_upstr, sigma_upstr);
    Real pB_downstr = pB_upstr * ratio_u * ratio_u;
    Real n_downstr = n_upstr * ratio_u;
    Real e_th = compute_downstr_eth(Gamma_rel, n_downstr);
    shock.t_comv(i, j, k) = state.t_comv;
    shock.r(i, j, k) = state.r;
    shock.theta(i, j, k) = state.theta;
    shock.Gamma(i, j, k) = Gamma_downstr;
    shock.Gamma_rel(i, j, k) = Gamma_rel;
    shock.B(i, j, k) = compute_comv_weibel_B(shock.eps_B, e_th) + std::sqrt(pB_downstr * 8 * con::pi);
    shock.column_num_den(i, j, k) = N_downstr / (state.r * state.r);
}

/********************************************************************************************************************
 * @brief Computes the swept-up mass for a shock based on the equation system and current state.
 * @tparam Eqn Type of the equation system
 * @param eqn The equation system containing medium properties and other parameters
 * @param state The current state of the system
 * @return The swept-up mass per solid angle
 * @details Handles both cases where mass profile is provided or calculated.
 ********************************************************************************************************************/
template <typename Eqn>
Real compute_swept_mass(Eqn const& eqn, typename Eqn::State const& state) {
    return eqn.medium.mass(eqn.phi, state.theta, state.r);
}

/********************************************************************************************************************
 * @brief Computes the deceleration time for the shock based on the equation system and time bounds.
 * @tparam Eqn Type of the equation system
 * @param eqn The equation system containing ejecta and medium properties
 * @param t_max Maximum time to consider
 * @return The estimated deceleration time
 * @details The deceleration time marks when the shock begins to significantly slow down due to mass sweeping.
 ********************************************************************************************************************/
template <typename Eqn>
Real compute_dec_time(Eqn const& eqn, Real t_max) {
    Real e_k = eqn.ejecta.eps_k(eqn.phi, eqn.theta0);
    Real gamma = eqn.ejecta.Gamma0(eqn.phi, eqn.theta0);
    Real beta = gamma_to_beta(gamma);

    Real m_shell = e_k / (gamma * con::c2);

    if constexpr (HasSigma<decltype(eqn.ejecta)>) {
        m_shell /= 1 + eqn.ejecta.sigma0(eqn.phi, eqn.theta0);
    }

    Real r_dec = beta * con::c * t_max / (1 - beta);
    for (size_t i = 0; i < 30; i++) {
        Real m_swept = eqn.medium.mass(eqn.phi, eqn.theta0, r_dec);
        if (m_swept < m_shell / gamma) {
            break;
        }
        r_dec /= 3;
    }
    Real t_dec = r_dec * (1 - beta) / (beta * con::c);
    return t_dec;
}