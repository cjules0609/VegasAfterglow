//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "reverse-shock.hpp"
#include "shock.h"
/********************************************************************************************************************
 * @brief Constructor for the FRShockEqn class.
 * @details Initializes the forward-reverse shock equation with the given medium, ejecta, and parameters.
 * @param medium The medium through which the shock propagates
 * @param ejecta The ejecta driving the shock
 * @param phi Azimuthal angle
 * @param theta Polar angle
 * @param eps_e Electron energy fraction
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
FRShockEqn<Ejecta, Medium>::FRShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, Real eps_e)
    : medium(medium),
      ejecta(ejecta),
      phi(phi),
      theta0(theta),
      eps_e(eps_e),
      Gamma4(ejecta.Gamma0(phi, theta)),
      deps0_dt(ejecta.eps_k(phi, theta) / ejecta.T0),
      dm0_dt(deps0_dt / (Gamma4 * con::c2)),
      u4(std::sqrt(Gamma4 * Gamma4 - 1) * con::c) {
    if constexpr (HasSigma<Ejecta>::value) {
        dm0_dt /= 1 + ejecta.sigma0(phi, theta);
    }
}

/********************************************************************************************************************
 * @brief Implements the reverse shock ODE system.
 * @details Computes the derivatives of state variables with respect to time.
 * @param state Current state of the system
 * @param diff Output derivatives to be populated
 * @param t Current time
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
void FRShockEqn<Ejecta, Medium>::operator()(State const& state, State& diff, Real t) {
    Real Gamma3 = 1;

    // Gamma_rel = 1;

    auto [deps_shell_dt, dm_shell_dt] = get_injection_rate(t);
    diff.eps_shell = deps_shell_dt;
    diff.m_shell = dm_shell_dt;

    bool is_injecting = diff.eps_shell > 0 || diff.m_shell > 0;

    if (!crossed) {
        Gamma3 = compute_crossing_Gamma3(state);
        // Gamma_rel = compute_rel_Gamma(Gamma4, Gamma3);

        Real sigma4 = compute_shell_sigma(state);
        diff.m3 = compute_dm3_dt(state.width_shell, state.m_shell, Gamma3, Gamma4, sigma4);
        if (state.m3 >= state.m_shell) {
            diff.m3 = std::min(diff.m3, diff.m_shell);
        }
    } else {
        Real Gamma_rel = compute_crossed_Gamma_rel(state);
        Gamma3 = compute_crossed_Gamma3(Gamma_rel, state.r);
        diff.m3 = 0;
    }

    Real beta3 = gamma_to_beta(Gamma3);
    diff.r = compute_dr_dt(beta3);
    diff.t_comv = compute_dt_dt_comv(Gamma3, beta3);

    diff.width_shell = is_injecting ? u4 : compute_shell_spreading_rate(Gamma3, diff.t_comv);
    diff.theta = 0;  // no lateral spreading
}
inline Real compute_init_comv_shell_width(Real Gamma4, Real t0, Real T);
/********************************************************************************************************************
 * @brief Set the initial conditions for the reverse shock ODE.
 * @details Sets up initial state values and determines if the shock has already crossed.
 * @param state State vector to initialize
 * @param t0 Initial time
 * @return True if the shock has already crossed at the initial time, false otherwise
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
bool FRShockEqn<Ejecta, Medium>::set_init_state(State& state, Real t0) const noexcept {
    Real beta4 = gamma_to_beta(Gamma4);
    state.theta = theta0;
    state.width_shell = compute_init_comv_shell_width(Gamma4, t0, ejecta.T0);

    Real dt = std::min(t0, ejecta.T0);
    state.eps_shell = deps0_dt * dt;
    state.m_shell = dm0_dt * dt;

    state.r = beta4 * con::c * t0 / (1 - beta4);
    state.t_comv = state.r / std::sqrt(Gamma4 * Gamma4 - 1) / con::c;

    Real rho = medium.rho(phi, theta0, state.r);

    if (t0 < ejecta.T0) {
        // thick shell deceleration radius where gamma starts to drop propto t^{-1/4}
        Real r_dec = 0.35 * std::sqrt(3 * dm0_dt * (1 - beta4) / (beta4 * con::c * rho * Gamma4));

        // r0 is larger than the decelertation radius, so r=beta4 * con::c * t0 / (1 - beta4) is not appropriate
        if (state.r > r_dec) {
            Real t_dec = r_dec * (1 - beta4) / (beta4 * con::c);
            state.r = std::sqrt(4 * Gamma4 * Gamma4 * r_dec * (t0 - t_dec) + r_dec * r_dec);
            rho = medium.rho(phi, theta0, state.r);

            Real t_dec_com = r_dec / std::sqrt(Gamma4 * Gamma4 - 1) / con::c;
            state.t_comv = (std::sqrt(state.r * state.r * state.r) - std::sqrt(r_dec * r_dec * r_dec)) /
                               (1.5 * Gamma4 * std::sqrt(r_dec)) +
                           t_dec_com;
        }
    }

    auto [crossed, m3] = solve_init_m3(*this, state, Gamma4, t0);
    state.m3 = m3;
    return crossed;
}

/********************************************************************************************************************
 * @brief Calculates the energy and mass injection rates at time t.
 * @details Returns a pair containing energy injection rate and mass injection rate.
 * @param t Current time
 * @return A pair {energy_rate, mass_rate} with injection rates at time t
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
std::pair<Real, Real> FRShockEqn<Ejecta, Medium>::get_injection_rate(Real t) const {
    Real deps_shell_dt = 0;
    Real dm_shell_dt = 0;

    if (t < ejecta.T0) {
        deps_shell_dt = deps0_dt;
        dm_shell_dt = dm0_dt;
    }

    if constexpr (State::energy_inject) {
        deps_shell_dt += ejecta.deps_dt(phi, theta0, t);
    }

    if constexpr (State::mass_inject) {
        dm_shell_dt += ejecta.dm_dt(phi, theta0, t);
    }

    return {deps_shell_dt, dm_shell_dt};
}

/********************************************************************************************************************
 * @brief Checks if the ejecta is still injecting mass or energy at time t.
 * @details Returns true if either energy or mass is being injected.
 * @param t Current time
 * @return Boolean indicating if injection is happening at time t
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
bool FRShockEqn<Ejecta, Medium>::is_injecting(Real t) const {
    auto [deps_shell_dt, dm_shell_dt] = get_injection_rate(t);
    return deps_shell_dt > 0 || dm_shell_dt > 0;
}

/********************************************************************************************************************
 * @brief Computes the Lorentz factor (Gamma3) during shock crossing.
 * @details Uses energy conservation to determine the appropriate Gamma3 value.
 * @param state Current state of the system
 * @return The computed Lorentz factor for region 3 (shocked ejecta)
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_crossing_Gamma3(State const& state) const {
    Real m2 = compute_swept_mass(*this, state);
    Real m3 = state.m3;
    Real sigma4 = compute_shell_sigma(state);

    auto func = [=, this](Real Gamma3) -> Real {
        Real Gamma34 = compute_rel_Gamma(Gamma4, Gamma3);
        Real adx3 = adiabatic_idx(Gamma34);
        Real adx2 = adiabatic_idx(Gamma3);
        Real g_eff2 = compute_effective_Gamma(adx2, Gamma3);
        Real g_eff3 = compute_effective_Gamma(adx3, Gamma3);

        Real eps2 = m2 * (Gamma3 - 1 + g_eff2 * (Gamma3 - 1));
        Real eps3 = m3 * (Gamma3 - Gamma4 + g_eff3 * (Gamma34 - 1)) * (1 + sigma4);

        return eps2 + eps3;
    };
    constexpr Real r_tol = 1e-3;
    return root_bisect(func, 1, Gamma4, r_tol);
}

/********************************************************************************************************************
 * @brief Sets constants at the shock crossing point for later calculations.
 * @details Stores key parameters and switches the ODE from crossing to crossed state.
 * @param state State at shock crossing
 * @param B Magnetic field at shock crossing
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
void FRShockEqn<Ejecta, Medium>::set_cross_state(State const& state, Real B) {
    this->r_x = state.r;
    constexpr Real n3_norm = 1;  // normalized n3/n3_x = 1
    this->N_electron = n3_norm * state.width_shell * state.r * state.r;

    Real Gamma3 = compute_crossing_Gamma3(state);
    this->u_x = std::sqrt(Gamma3 * Gamma3 - 1);

    Real Gamma_rel = compute_rel_Gamma(Gamma4, Gamma3);
    Real gamma_electron = (Gamma_rel - 1) * con::mp / con::me * eps_e + 1;  // electron Lorentz factor

    // use electron adiabatic index see section 3.2 https://arxiv.org/pdf/astro-ph/9910241
    this->gamma_hat_x = adiabatic_idx(gamma_electron);
    Real p_norm = (gamma_hat_x - 1) * (Gamma_rel - 1) * n3_norm;      // p = (ad-1)(gamma-1)n m_p c^2
    this->adiabatic_const = fast_pow(n3_norm, gamma_hat_x) / p_norm;  // p \propto n^ad
    this->e_mag_const = p_norm / (B * B);
    this->crossed = true;
}

/********************************************************************************************************************
 * @brief Calculates the power-law index for post-crossing four-velocity evolution.
 * @details The index transitions from g_low=1.5 for low relative Lorentz factors to g_high=3.5
 *          for high relative Lorentz factors (Blandford-McKee limit).
 * @param gamma_rel Relative Lorentz factor
 * @param k Medium power law index (default: 0)
 * @return The power-law index for velocity evolution
 ********************************************************************************************************************/
inline Real get_post_cross_g(Real gamma_rel, Real k = 0) {
    constexpr Real g_low = 1.5;   // k is the medium power law index
    constexpr Real g_high = 3.5;  // Blandford-McKee limit// TODO: need to be modified for non ISM medium
    Real p = std::sqrt(std::sqrt(gamma_rel - 1));
    return g_low + (g_high - g_low) * p / (1 + p);
}

/********************************************************************************************************************
 * @brief Computes the Lorentz factor for region 3 after shock crossing.
 * @details Uses a power-law profile that applies both in relativistic and Newtonian regimes.
 * @param gamma_rel Relative Lorentz factor
 * @param r Current radius
 * @return The Lorentz factor for region 3 post-crossing
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_crossed_Gamma3(Real gamma_rel, Real r) const {
    Real g = get_post_cross_g(gamma_rel);
    Real u = u_x * fast_pow(r / r_x, -g);
    return std::sqrt(u * u + 1);
}

/********************************************************************************************************************
 * @brief Computes the relative Lorentz factor after shock crossing.
 * @details Uses adiabatic expansion to determine the relative Lorentz factor.
 * @param state Current state of the system
 * @return The relative Lorentz factor post-crossing
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_crossed_Gamma_rel(State const& state) const {
    Real n3_norm = N_electron / (state.width_shell * state.r * state.r);  // proton number conservation
    Real p3_norm = std::pow(n3_norm, gamma_hat_x) / adiabatic_const;      // adiabatic expansion
    return p3_norm / ((gamma_hat_x - 1) * n3_norm) + 1;
}

/********************************************************************************************************************
 * @brief Computes the magnetic field after shock crossing.
 * @details Calculates magnetic field based on pressure and energy conservation.
 * @param state Current state of the system
 * @return Magnetic field strength in the shocked region
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_crossed_B(State const& state) const {
    Real n3_norm = N_electron / (state.width_shell * state.r * state.r);  // proton number conservation
    Real p3_norm = std::pow(n3_norm, gamma_hat_x) / adiabatic_const;      // adiabatic expansion
    return std::sqrt(p3_norm / e_mag_const);
}

/********************************************************************************************************************
 * @brief Calculates the magnetization parameter of the shell.
 * @details Sigma is defined as (ε/Γmc²) - 1, where ε is the energy per solid angle.
 * @param state Current state of the system
 * @return The magnetization parameter of the shell
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_shell_sigma(State const& state) const {
    return std::max(0.0, state.eps_shell / (Gamma4 * state.m_shell * con::c2) - 1);
}

//---------------------------------------------------------------------------------------------------------------------
// Helper functions
//---------------------------------------------------------------------------------------------------------------------

/********************************************************************************************************************
 * @brief Calculates the initial shocked mass m3 at time t0.
 * @details Ensures energy conservation at the initial time.
 * @param eqn The equation system containing medium and ejecta parameters
 * @param state Current state of the system
 * @param Gamma4 Lorentz factor of the unshocked ejecta
 * @param t0 Initial time
 * @return A tuple containing {crossed, initial_m3} where crossed indicates if shock has crossed
 ********************************************************************************************************************/
template <typename Eqn, typename State>
inline auto solve_init_m3(Eqn const& eqn, State const& state, Real Gamma4, Real t0) {
    Real E_iso = eqn.ejecta.eps_k(eqn.phi, eqn.theta0) * 4 * con::pi;
    Real n1 = eqn.medium.rho(eqn.phi, eqn.theta0, state.r) / con::mp;
    Real l = sedov_length(E_iso, n1);
    Real sigma0 = eqn.compute_shell_sigma(state);
    Real beta4 = gamma_to_beta(Gamma4);
    Real Delta0 = con::c * eqn.ejecta.T0;  // lab frame shell width
    Real r_x = std::sqrt(std::sqrt(Delta0 * l * l * l / 3) / (1 + sigma0) * 3. / 4);
    Real Rs = con::c * eqn.ejecta.T0 * beta4 / (1 - beta4);

    if (r_x > Rs) {
        constexpr Real rtol = 1e-3;
        Real r_x_new = (1 + 2 * rtol) * r_x;
        for (; std::fabs((r_x - r_x_new) / r_x) > rtol;) {
            Delta0 = con::c * eqn.ejecta.T0 + (r_x_new - Rs) / (std::sqrt(3) * Gamma4 * Gamma4);
            r_x = r_x_new;
            r_x_new = std::sqrt(std::sqrt(Delta0 * l * l * l / 3) / (1 + sigma0) * 3. / 4);
        }
    }

    Real m4_tot = E_iso / (4 * con::pi * con::c2 * Gamma4 * (1 + sigma0));

    if (!eqn.is_injecting(t0)) {
        Real m3 = m4_tot * std::pow(state.r / r_x, 1.5);
        if (m3 >= state.m_shell) {
            return std::make_tuple(true, state.m_shell);
        } else {
            return std::make_tuple(false, m3);
        }
    } else {
        Real m3 = m4_tot * state.r / r_x;
        return std::make_tuple(false, std::min(m3, state.m_shell));
    }
}

/********************************************************************************************************************
 * @brief Calculates the comoving shell width at initial radius.
 * @details Accounts for both pure injection phase and shell spreading phase.
 * @param Gamma4 Lorentz factor of the unshocked ejecta
 * @param t0 Initial time
 * @param T Engine duration
 * @return The comoving shell width
 ********************************************************************************************************************/
inline Real compute_init_comv_shell_width(Real Gamma4, Real t0, Real T) {
    Real beta4 = gamma_to_beta(Gamma4);
    if (t0 < T) {  // pure injection
        return Gamma4 * t0 * beta4 * con::c;
    } else {  // injection+shell spreading
        Real cs = con::c / std::sqrt(3);
        return Gamma4 * T * beta4 * con::c + cs * (t0 - T) * Gamma4;
    }
}
// inline Real initShellWidth(Real gamma, Real T, Real r) { return r / gamma / std::sqrt(3); }

/********************************************************************************************************************
 * @brief Updates the post-shock crossing state at a grid point.
 * @details Stores computed values in the shock object for the given grid indices.
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param k Grid index for time
 * @param eqn The equation system containing parameters
 * @param state Current state of the system
 * @param shock Shock object to update
 ********************************************************************************************************************/
template <typename Eqn, typename State>
inline void save_crossed_RS_state(size_t i, size_t j, size_t k, Eqn const& eqn, State const& state, Shock& shock) {
    size_t k0 = shock.injection_idx(i, j);
    Real r0 = shock.r(i, j, k0);
    shock.t_comv(i, j, k) = state.t_comv;
    shock.r(i, j, k) = state.r;
    shock.theta(i, j, k) = state.theta;
    shock.Gamma_rel(i, j, k) = eqn.compute_crossed_Gamma_rel(state);
    shock.Gamma(i, j, k) = eqn.compute_crossed_Gamma3(shock.Gamma_rel(i, j, k), state.r);
    shock.column_num_den(i, j, k) = shock.column_num_den(i, j, k0) * (r0 * r0) / (state.r * state.r);
    shock.B(i, j, k) = eqn.compute_crossed_B(state);
}

/********************************************************************************************************************
 * @brief Sets the initial condition for forward shock based on the reverse shock state.
 * @details Used at the crossing point to maintain consistency between forward and reverse shocks.
 * @param eqn_rvs The reverse shock equation system
 * @param state_fwd Forward shock state to initialize
 * @param state_rvs Reverse shock state at crossing
 * @param gamma2 Lorentz factor of region 2 (shocked medium)
 ********************************************************************************************************************/
template <typename Eqn, typename FState, typename RState>
void set_fwd_state_from_rvs_state(Eqn const& eqn_rvs, FState& state_fwd, RState const& state_rvs, Real gamma2) {
    state_fwd.r = state_rvs.r;
    state_fwd.t_comv = state_rvs.t_comv;
    state_fwd.theta = state_rvs.theta;
    state_fwd.Gamma = gamma2;

    Real m_swept = compute_swept_mass(eqn_rvs, state_rvs);

    if constexpr (HasU<FState>::value) {
        state_fwd.u = (gamma2 - 1) * m_swept * con::c2;
    }

    if constexpr (FState::energy_inject) {
        state_fwd.eps_shell = state_rvs.eps_shell;
    }

    if constexpr (FState::mass_inject) {
        state_fwd.m_shell = state_rvs.m_shell;
    }
}

/********************************************************************************************************************
 * @brief Saves the state of both forward and reverse shocks at a grid point.
 * @details Updates shock properties for both shocks and checks if crossing is complete.
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param k Grid index for time
 * @param eqn_rvs The reverse shock equation system
 * @param state Current state of the system
 * @param t Current time
 * @param shock_fwd Forward shock object to update
 * @param shock_rvs Reverse shock object to update
 * @return True if the shock has crossed and there's no more injection
 ********************************************************************************************************************/
template <typename Eqn, typename State>
bool save_shock_pair_state(size_t i, size_t j, int k, Eqn const& eqn_rvs, State const& state, Real t, Shock& shock_fwd,
                           Shock& shock_rvs) {
    Real n4 = state.m_shell / (state.r * state.r * state.width_shell * con::mp);
    Real sigma4 = eqn_rvs.compute_shell_sigma(state);

    Real m2 = compute_swept_mass(eqn_rvs, state);
    Real n1 = eqn_rvs.medium.rho(eqn_rvs.phi, state.theta, state.r) / con::mp;

    constexpr Real Gamma1 = 1;
    constexpr Real sigma1 = 0;

    Real Gamma3 = eqn_rvs.compute_crossing_Gamma3(state);

    save_shock_state(shock_fwd, i, j, k, state, Gamma3, Gamma1, m2 / con::mp, n1, sigma1);
    save_shock_state(shock_rvs, i, j, k, state, Gamma3, eqn_rvs.Gamma4, state.m3 / con::mp, n4, sigma4);
    return state.m3 >= state.m_shell && !eqn_rvs.is_injecting(t);
}

/********************************************************************************************************************
 * @brief Handles low Gamma scenario by stopping both shocks.
 * @details Sets appropriate values when Lorentz factor drops below threshold.
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param shock_fwd Forward shock object to update
 * @param shock_rvs Reverse shock object to update
 * @param state_fwd Forward shock state
 * @param state_rvs Reverse shock state
 * @return True if stopping condition was applied
 ********************************************************************************************************************/
template <typename FState, typename RState>
bool set_stopping_shocks(size_t i, size_t j, Shock& shock_fwd, Shock& shock_rvs, FState const& state_fwd,
                         RState const& state_rvs) {
    if (state_fwd.Gamma < con::Gamma_cut) {
        set_stopping_shock(i, j, shock_fwd, state_fwd);
        set_stopping_shock(i, j, shock_rvs, state_rvs);
        return true;
    }
    return false;
}

/********************************************************************************************************************
 * @brief Solves forward/reverse shock until crossing occurs.
 * @details Advances the ODE solver and saves shock states at each time step.
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param t View of time points
 * @param stepper_rvs ODE stepper for reverse shock
 * @param eqn_rvs Reverse shock equation system
 * @param state_rvs Reverse shock state
 * @param shock_fwd Forward shock object
 * @param shock_rvs Reverse shock object
 * @return Index where shock crossing occurs
 ********************************************************************************************************************/
template <typename View, typename RState, typename REqn, typename Stepper>
size_t solve_until_cross(size_t i, size_t j, View const& t, Stepper& stepper_rvs, REqn& eqn_rvs, RState& state_rvs,
                         Shock& shock_fwd, Shock& shock_rvs) {
    size_t k0 = 0;
    Real t_back = t.back();
    bool crossed = false;
    size_t k = 0;

    while (!crossed && stepper_rvs.current_time() <= t_back) {
        stepper_rvs.do_step(eqn_rvs);
        while (k < t.size() && stepper_rvs.current_time() > t(k)) {
            stepper_rvs.calc_state(t(k), state_rvs);
            crossed = save_shock_pair_state(i, j, k, eqn_rvs, state_rvs, t(k), shock_fwd, shock_rvs);
            if (crossed) {
                k0 = k;
                shock_rvs.injection_idx(i, j) = k0;
                eqn_rvs.set_cross_state(state_rvs, shock_rvs.B(i, j, k));
                stepper_rvs.initialize(state_rvs, t(k0), stepper_rvs.current_time_step());
                break;
            }
            k++;
        }
    }

    return k0;
}

/********************************************************************************************************************
 * @brief Solves shock after crossing using the provided update function.
 * @details Advances the ODE solver and saves shock states at each time step.
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param t View of time points
 * @param k0 Index at which crossing occurred
 * @param stepper ODE stepper
 * @param eqn Equation system
 * @param state Current state
 * @param shock Shock object to update
 * @param update Function to update shock state
 ********************************************************************************************************************/
template <typename UpdateFunc, typename Stepper, typename Eqn, typename State, typename View>
void solve_post_cross(size_t i, size_t j, View const& t, size_t k0, Stepper& stepper, Eqn& eqn, State& state,
                      Shock& shock, UpdateFunc update) {
    size_t k = k0 + 1;
    Real t_back = t.back();

    while (stepper.current_time() <= t_back) {
        stepper.do_step(eqn);
        while (k < t.size() && stepper.current_time() > t(k)) {
            stepper.calc_state(t(k), state);
            update(i, j, k, eqn, state, shock);
            k++;
        }
    }
}

/********************************************************************************************************************
 * @brief Solves the reverse/forward shock ODE at a grid point.
 * @details Manages the evolution of both shocks before and after crossing.
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param t View of time points
 * @param shock_fwd Forward shock object
 * @param shock_rvs Reverse shock object
 * @param eqn_fwd Forward shock equation system
 * @param eqn_rvs Reverse shock equation system
 * @param rtol Relative tolerance for ODE solver
 ********************************************************************************************************************/
template <typename FwdEqn, typename RvsEqn, typename View>
void grid_solve_shock_pair(size_t i, size_t j, View const& t, Shock& shock_fwd, Shock& shock_rvs, FwdEqn const& eqn_fwd,
                           RvsEqn& eqn_rvs, Real rtol = 1e-6) {
    using namespace boost::numeric::odeint;
    auto stepper_fwd = make_dense_output(rtol, rtol, runge_kutta_dopri5<typename FwdEqn::State>());

    typename FwdEqn::State state_fwd;
    typename RvsEqn::State state_rvs;

    Real t_dec = compute_dec_time(eqn_rvs, t.back());
    Real t0 = std::min(t.front(), t_dec / 100);

    eqn_fwd.set_init_state(state_fwd, t0);

    bool crossed = eqn_rvs.set_init_state(state_rvs, t0);

    if (set_stopping_shocks(i, j, shock_fwd, shock_rvs, state_fwd, state_rvs)) return;

    if (crossed) {
        grid_solve_fwd_shock(i, j, t, shock_fwd, eqn_fwd, rtol);
        return;
    }

    auto stepper_rvs = make_dense_output(rtol, rtol, runge_kutta_dopri5<typename RvsEqn::State>());
    stepper_rvs.initialize(state_rvs, t0, 0.1 * t0);

    size_t k0 = solve_until_cross(i, j, t, stepper_rvs, eqn_rvs, state_rvs, shock_fwd, shock_rvs);

    set_fwd_state_from_rvs_state(eqn_rvs, state_fwd, state_rvs, shock_rvs.Gamma(i, j, k0));
    stepper_fwd.initialize(state_fwd, t(k0), stepper_rvs.current_time_step());

    solve_post_cross(i, j, t, k0, stepper_fwd, eqn_fwd, state_fwd, shock_fwd,
                     save_fwd_shock_state<FwdEqn, typename FwdEqn::State>);

    solve_post_cross(i, j, t, k0, stepper_rvs, eqn_rvs, state_rvs, shock_rvs,
                     save_crossed_RS_state<RvsEqn, typename RvsEqn::State>);
}

/********************************************************************************************************************
 * @brief Generates a pair of forward and reverse shocks.
 * @details Creates two Shock objects and solves the shock evolution for each grid point.
 * @param coord Coordinate system definition
 * @param medium The medium through which the shock propagates
 * @param jet The jet (ejecta) driving the shock
 * @param eps_e_f Electron energy fraction for forward shock
 * @param eps_B_f Magnetic energy fraction for forward shock
 * @param eps_e_r Electron energy fraction for reverse shock
 * @param eps_B_r Magnetic energy fraction for reverse shock
 * @param rtol Relative tolerance for ODE solver
 * @return A pair of Shock objects {forward_shock, reverse_shock}
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
ShockPair generate_shock_pair(Coord const& coord, Medium const& medium, Ejecta const& jet, Real eps_e_f, Real eps_B_f,
                              Real eps_e_r, Real eps_B_r, Real rtol = 1e-5) {
    auto [phi_size, theta_size, t_size] = coord.shape();
    size_t phi_size_needed = coord.t.shape()[0];
    Shock f_shock(phi_size_needed, theta_size, t_size, eps_e_f, eps_B_f);
    Shock r_shock(phi_size_needed, theta_size, t_size, eps_e_r, eps_B_r);
    for (size_t i = 0; i < phi_size_needed; ++i) {
        Real theta_s =
            jet_spreading_edge(jet, medium, coord.phi(i), coord.theta.front(), coord.theta.back(), coord.t.front());
        for (size_t j = 0; j < theta_size; ++j) {
            // auto eqn_f = ForwardShockEqn(medium, jet, inject, coord.phi(i), coord.theta(j), eps_e);
            auto eqn_f = SimpleShockEqn(medium, jet, coord.phi(i), coord.theta(j), eps_e_f, theta_s);
            auto eqn_r = FRShockEqn(medium, jet, coord.phi(i), coord.theta(j), eps_e_r);
            // Solve the forward-reverse shock shell
            grid_solve_shock_pair(i, j, xt::view(coord.t, i, j, xt::all()), f_shock, r_shock, eqn_f, eqn_r, rtol);
        }
    }

    return std::make_pair(std::move(f_shock), std::move(r_shock));
}
