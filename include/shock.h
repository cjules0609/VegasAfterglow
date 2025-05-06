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
 * CLASS: Shock
 * DESCRIPTION: Represents a shock wave in an astrophysical environment. The class stores physical properties of the
 *              shock across a 3D grid defined by azimuthal angle (phi), polar angle (theta), and time bins.
 *              Provides methods for shock calculations, including relativistic jump conditions, magnetic field
 *              calculations, and energy density computations.
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
 * INLINE FUNCTIONS: Shock Utilities
 * DESCRIPTION: This section defines a set of inline functions used in shock calculations. These functions compute
 *              various physical quantities such as the comoving magnetic field (via the Weibel instability),
 *              thermal energy density, time derivatives, jet width derivative, downstream number density, fluid
 *              velocities, and update the shock state.
 ********************************************************************************************************************/

/********************************************************************************************************************
 * FUNCTION: compute_downstr_4vel(gamma_rel, sigma)
 * DESCRIPTION: Computes the downstream four-velocity for a given relative Lorentz factor and magnetization parameter.
 *              The calculation handles both magnetized (sigma > 0) and non-magnetized (sigma = 0) cases using
 *              different equations based on jump conditions across the shock front.
 * PARAMETERS:
 *   - gamma_rel: Relative Lorentz factor between upstream and downstream regions
 *   - sigma: Magnetization parameter (ratio of magnetic to rest-mass energy density)
 * RETURNS: The downstream four-velocity in the shock frame
 ********************************************************************************************************************/
Real compute_downstr_4vel(Real gamma_rel, Real sigma);

/********************************************************************************************************************
 * FUNCTION: compute_4vel_jump(gamma_rel, sigma)
 * DESCRIPTION: Computes the ratio of upstream to downstream four-velocity across the shock front.
 *              This ratio is a key parameter in determining various shock properties, such as compression ratio
 *              and jump conditions for density, pressure, and magnetic field.
 * PARAMETERS:
 *   - gamma_rel: Relative Lorentz factor between upstream and downstream regions
 *   - sigma: Magnetization parameter (ratio of magnetic to rest-mass energy density)
 * RETURNS: The ratio of upstream to downstream four-velocity
 ********************************************************************************************************************/
Real compute_4vel_jump(Real gamma_rel, Real sigma);

/********************************************************************************************************************
 * FUNCTION: compute_sound_speed(Gamma_rel)
 * DESCRIPTION: Computes the sound speed in the shocked medium based on the relative Lorentz factor.
 * PARAMETERS:
 *   - Gamma_rel: Relative Lorentz factor between upstream and downstream regions
 * RETURNS: The sound speed as a fraction of light speed
 ********************************************************************************************************************/
inline Real compute_sound_speed(Real Gamma_rel) {
    Real ad_idx = adiabatic_idx(Gamma_rel);
    return std::sqrt(ad_idx * (ad_idx - 1) * (Gamma_rel - 1) / (1 + (Gamma_rel - 1) * ad_idx));
}

/********************************************************************************************************************
 * FUNCTION: compute_effective_Gamma(adx, Gamma)
 * DESCRIPTION: Computes the effective Lorentz factor accounting for the adiabatic index.
 * PARAMETERS:
 *   - adx: Adiabatic index of the medium
 *   - Gamma: Bulk Lorentz factor
 * RETURNS: The effective Lorentz factor
 ********************************************************************************************************************/
inline Real compute_effective_Gamma(Real adx, Real Gamma) { return (adx * Gamma * Gamma - adx + 1) / Gamma; }

/********************************************************************************************************************
 * FUNCTION: compute_comv_weibel_B(eps_B, e_thermal)
 * DESCRIPTION: Computes the comoving magnetic field using the Weibel instability mechanism.
 * PARAMETERS:
 *   - eps_B: Fraction of thermal energy in magnetic fields
 *   - e_thermal: Thermal energy density
 * RETURNS: The comoving magnetic field strength
 ********************************************************************************************************************/
inline Real compute_comv_weibel_B(Real eps_B, Real e_thermal) { return std::sqrt(8 * con::pi * eps_B * e_thermal); }

/********************************************************************************************************************
 * FUNCTION: compute_dr_dt(beta)
 * DESCRIPTION: Computes the time derivative of radius (dr/dt) based on the shock velocity.
 * PARAMETERS:
 *   - beta: Shock velocity as a fraction of light speed
 * RETURNS: The rate of change of radius with respect to observer time
 ********************************************************************************************************************/
inline Real compute_dr_dt(Real beta) { return (beta * con::c) / (1 - beta); }

/********************************************************************************************************************
 * FUNCTION: compute_dtheta_dt(theta_s, theta, drdt, r, Gamma)
 * DESCRIPTION: Computes the time derivative of theta (dθ/dt) for jet spreading.
 * PARAMETERS:
 *   - theta_s: typical spreading angle of the jet
 *   - theta: theta of the current grid point
 *   - drdt: Time derivative of radius
 *   - r: Current radius
 *   - Gamma: Current bulk Lorentz factor
 * RETURNS: The rate of change of the half-opening angle
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
 * FUNCTION: compute_dt_dt_comv(Gamma)
 * DESCRIPTION: Computes the time derivative of comoving time (dt_comv/dt) based on the Lorentz factor.
 * PARAMETERS:
 *   - Gamma: Bulk Lorentz factor
 *   - beta: Bulk velocity as a fraction of light speed
 * RETURNS: The rate of change of comoving time with respect to observer time
 ********************************************************************************************************************/
inline Real compute_dt_dt_comv(Real Gamma, Real beta) { return 1 / (Gamma * (1 - beta)); };

/********************************************************************************************************************
 * FUNCTION: compute_upstr_mag_p(n_up, sigma)
 * DESCRIPTION: Computes the upstream magnetic pressure based on number density and magnetization.
 * PARAMETERS:
 *   - n_up: Upstream number density
 *   - sigma: Magnetization parameter
 * RETURNS: The upstream magnetic pressure
 ********************************************************************************************************************/
inline Real compute_upstr_mag_p(Real n_up, Real sigma) { return sigma * n_up * con::mp * con::c2 / 2; }

/********************************************************************************************************************
 * FUNCTION: compute_upstr_4vel(u_down, gamma_rel)
 * DESCRIPTION: Computes the upstream four-velocity from downstream four-velocity and relative Lorentz factor.
 * PARAMETERS:
 *   - u_down: Downstream four-velocity
 *   - gamma_rel: Relative Lorentz factor
 * RETURNS: The upstream four-velocity in the shock frame
 ********************************************************************************************************************/
inline Real compute_upstr_4vel(Real u_down, Real gamma_rel) {
    return std::sqrt((1 + u_down * u_down) * (gamma_rel * gamma_rel - 1)) + u_down * gamma_rel;
}

/********************************************************************************************************************
 * FUNCTION: compute_downstr_num_den(n_up_str, gamma_rel, sigma)
 * DESCRIPTION: Computes the downstream number density from upstream density, relative Lorentz factor,
 *              and magnetization using the shock jump conditions.
 * PARAMETERS:
 *   - n_up_str: Upstream number density
 *   - gamma_rel: Relative Lorentz factor
 *   - sigma: Magnetization parameter
 * RETURNS: The downstream number density
 ********************************************************************************************************************/
inline Real compute_downstr_num_den(Real n_up_str, Real gamma_rel, Real sigma) {
    return n_up_str * compute_4vel_jump(gamma_rel, sigma);
}

/********************************************************************************************************************
 * FUNCTION: compute_rel_Gamma(gamma1, gamma2)
 * DESCRIPTION: Computes the relative Lorentz factor between two frames with given Lorentz factors.
 * PARAMETERS:
 *   - gamma1: First frame's Lorentz factor
 *   - gamma2: Second frame's Lorentz factor
 * RETURNS: The relative Lorentz factor between the two frames
 ********************************************************************************************************************/
inline Real compute_rel_Gamma(Real gamma1, Real gamma2) {
    return gamma1 * gamma2 - std::sqrt((gamma1 * gamma1 - 1) * (gamma2 * gamma2 - 1));
}

/********************************************************************************************************************
 * FUNCTION: compute_rel_Gamma(gamma1, gamma2, beta1, beta2)
 * DESCRIPTION: Computes the relative Lorentz factor between two frames with given Lorentz factors and velocities.
 * PARAMETERS:
 *   - gamma1: First frame's Lorentz factor
 *   - gamma2: Second frame's Lorentz factor
 *   - beta1: First frame's velocity as fraction of light speed
 *   - beta2: Second frame's velocity as fraction of light speed
 * RETURNS: The relative Lorentz factor between the two frames
 ********************************************************************************************************************/
inline Real compute_rel_Gamma(Real gamma1, Real gamma2, Real beta1, Real beta2) {
    return gamma1 * gamma2 * (1 - beta1 * beta2);
}

/********************************************************************************************************************
 * FUNCTION: compute_Gamma_from_relative(gamma4, gamma_rel)
 * DESCRIPTION: Computes a Lorentz factor from a reference Lorentz factor and relative Lorentz factor.
 * PARAMETERS:
 *   - gamma4: Reference Lorentz factor (typically for region 4, unshocked ejecta)
 *   - gamma_rel: Relative Lorentz factor
 * RETURNS: The derived Lorentz factor
 ********************************************************************************************************************/
inline Real compute_Gamma_from_relative(Real gamma4, Real gamma_rel) {
    Real b = -2 * gamma4 * gamma_rel;
    Real c = gamma4 * gamma4 + gamma_rel * gamma_rel - 1;
    return (-b - std::sqrt(b * b - 4 * c)) / 2;
}

/********************************************************************************************************************
 * FUNCTION: compute_downstr_eth(gamma_rel, n_down_str)
 * DESCRIPTION: Computes the downstream thermal energy density.
 * PARAMETERS:
 *   - gamma_rel: Relative Lorentz factor
 *   - n_down_str: Downstream number density
 * RETURNS: The thermal energy density in the downstream region
 ********************************************************************************************************************/
inline Real compute_downstr_eth(Real gamma_rel, Real n_down_str) {
    return n_down_str * (gamma_rel - 1) * con::mp * con::c2;
}

/********************************************************************************************************************
 * FUNCTION: compute_shell_spreading_rate(Gamma_rel, dtdt_com)
 * DESCRIPTION: Computes the rate at which the shock shell spreads in the comoving frame.
 * PARAMETERS:
 *   - Gamma_rel: Relative Lorentz factor
 *   - dtdt_com: Rate of change of comoving time with respect to observer time
 * RETURNS: The shell spreading rate in the comoving frame
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
 * FUNCTION: compute_dm3_dt(width, m4, gamma3, gamma4, sigma)
 * DESCRIPTION: Computes the rate of change of shocked ejecta mass per solid angle.
 * PARAMETERS:
 *   - width: Width of the shocked region
 *   - m4: Mass per solid angle in region 4 (unshocked ejecta)
 *   - gamma3: Lorentz factor in region 3 (shocked ejecta)
 *   - gamma4: Lorentz factor in region 4 (unshocked ejecta)
 *   - sigma: Magnetization parameter
 * RETURNS: The rate of change of mass in region 3
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
 * FUNCTION: compute_region4_num_den(dEdOmega, Gamma0, r, D_jet, sigma)
 * DESCRIPTION: Computes the number density in region 4 (unshocked ejecta).
 * PARAMETERS:
 *   - dEdOmega: Energy per solid angle
 *   - Gamma0: Initial Lorentz factor
 *   - r: Radius
 *   - D_jet: Jet thickness
 *   - sigma: Magnetization parameter
 * RETURNS: The number density in region 4
 ********************************************************************************************************************/
inline Real compute_region4_num_den(Real dEdOmega, Real Gamma0, Real r, Real D_jet, Real sigma) {
    return dEdOmega / ((Gamma0 * con::mp * con::c2 * r * r * D_jet) * (1 + sigma));
}

/********************************************************************************************************************
 * FUNCTION: set_stopping_shock(i, j, shock, state0)
 * DESCRIPTION: Sets a stopping shock state when the Lorentz factor drops below threshold.
 * PARAMETERS:
 *   - i, j: Grid indices for phi and theta
 *   - shock: Reference to the Shock object to be updated
 *   - state0: Initial state to be used for some parameters
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
 * FUNCTION: save_shock_state(shock, i, j, k, state, Gamma_downstr, Gamma_upstr, N_downstr, n_upstr, sigma_upstr)
 * DESCRIPTION: Saves the shock state at a given grid point (i,j,k).
 * PARAMETERS:
 *   - shock: Reference to the Shock object to store the state
 *   - i, j, k: Grid indices for phi, theta, and time
 *   - state: Current state of the system
 *   - Gamma_downstr: Downstream Lorentz factor
 *   - Gamma_upstr: Upstream Lorentz factor
 *   - N_downstr: Downstream total number of shocked particles
 *   - n_upstr: Upstream number density
 *   - sigma_upstr: Upstream magnetization
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
 * FUNCTION: compute_swept_mass(eqn, state)
 * DESCRIPTION: Computes the swept-up mass for a shock based on the equation system and current state.
 *              Handles both cases where mass profile is provided or calculated.
 * PARAMETERS:
 *   - eqn: The equation system containing medium properties and other parameters
 *   - state: The current state of the system
 * RETURNS: The swept-up mass per solid angle
 ********************************************************************************************************************/
template <typename Eqn>
Real compute_swept_mass(Eqn const& eqn, typename Eqn::State const& state) {
    return eqn.medium.mass(eqn.phi, state.theta, state.r);
}

/********************************************************************************************************************
 * FUNCTION: compute_dec_time(eqn, t0, t_max)
 * DESCRIPTION: Computes the deceleration time for the shock based on the equation system and time bounds.
 *              The deceleration time marks when the shock begins to significantly slow down due to mass sweeping.
 * PARAMETERS:
 *   - eqn: The equation system containing ejecta and medium properties
 *   - t0: Initial time
 *   - t_max: Maximum time to consider
 * RETURNS: The estimated deceleration time
 ********************************************************************************************************************/
template <typename Eqn>
Real compute_dec_time(Eqn const& eqn, Real t_max) {
    Real e_k = eqn.ejecta.eps_k(eqn.phi, eqn.theta0);
    Real gamma = eqn.ejecta.Gamma0(eqn.phi, eqn.theta0);
    Real beta = gamma_to_beta(gamma);

    Real m_shell = e_k / (gamma * con::c2);

    if constexpr (HasSigma<decltype(eqn.ejecta)>::value) {
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