//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "reverse-shock.hpp"
#include "shock.h"

template <typename State>
bool is_injecting(State const& diff) {
    if (diff.eps4 > 0 || diff.m4 > 0) {
        return true;
    } else {
        return false;
    }
}

template <typename State>
bool is_crossing(State const& state, State const& diff) {
    return (state.m3 < state.m4) || is_injecting(diff);
}

template <typename Eqn, typename State>
bool is_crossing(Eqn const& eqn, State const& state, Real t) {
    Real dedt = eqn.ejecta.deps_dt(eqn.phi, state.theta, t);
    Real dmdt = eqn.ejecta.dm_dt(eqn.phi, state.theta, t);

    return (state.m3 < state.m4) || dedt > 0 || dmdt > 0;
}

template <typename Ejecta, typename Medium>
FRShockEqn<Ejecta, Medium>::FRShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta,
                                       RadParams const& rad_fwd, RadParams const& rad_rvs)
    : medium(medium),
      ejecta(ejecta),
      rad_fwd(rad_fwd),
      rad_rvs(rad_rvs),
      phi(phi),
      theta0(theta),
      Gamma4(ejecta.Gamma0(phi, theta)),
      deps0_dt(ejecta.eps_k(phi, theta) / ejecta.T0),
      dm0_dt(deps0_dt / (Gamma4 * con::c2)),
      u4(std::sqrt(Gamma4 * Gamma4 - 1) * con::c) {
    if constexpr (HasSigma<Ejecta>) {
        dm0_dt /= 1 + ejecta.sigma0(phi, theta);
    }
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_dGamma_dt(State const& state, State const& diff, Real t) const noexcept {
    /*Real Gamma34 = compute_rel_Gamma(Gamma4, state.Gamma);
    Real ad_idx2 = adiabatic_idx(state.Gamma);
    Real ad_idx3 = adiabatic_idx(Gamma34);

    Real Gamma_eff2 = compute_effective_Gamma(ad_idx2, state.Gamma);
    Real Gamma_eff3 = compute_effective_Gamma(ad_idx3, state.Gamma);

    Real dGamma_eff2_dGamma = compute_effective_Gamma_dGamma(ad_idx2, state.Gamma);
    Real dGamma_eff3_dGamma = compute_effective_Gamma_dGamma(ad_idx3, state.Gamma);

    Real rho = medium.rho(phi, state.theta, state.r);
    Real m2 = medium.mass(phi, state.theta, state.r);

    Real dm2_dt = state.r * state.r * rho * diff.r;

    Real term2 = Gamma_eff2 * (ad_idx2 - 1);
    Real term3 = Gamma_eff3 * (ad_idx3 - 1);

    Real a = (Gamma_eff2 + 1) * (state.Gamma - 1) * con::c2 * dm2_dt +
             (state.Gamma - Gamma4 + Gamma_eff3 * (Gamma34 - 1)) * con::c2 * diff.m3 -
             3 / state.r * (term2 * state.U_th2 + term3 * state.U_th3) * diff.r;
    Real b = (m2 + state.m3) * con::c2 + (dGamma_eff2_dGamma + term2 / state.Gamma) * state.U_th2 +
             (dGamma_eff3_dGamma + term3 / state.Gamma) * state.U_th3;

    return -a / b;*/

    Real Gamma34 = compute_rel_Gamma(Gamma4, state.Gamma);
    Real ad_idx2 = adiabatic_idx(state.Gamma);
    Real ad_idx3 = adiabatic_idx(Gamma34);

    Real Gamma_eff2 = compute_effective_Gamma(ad_idx2, state.Gamma);
    Real Gamma_eff3 = compute_effective_Gamma(ad_idx3, state.Gamma);

    Real dGamma_eff2_dGamma = compute_effective_Gamma_dGamma(ad_idx2, state.Gamma);
    Real dGamma_eff3_dGamma = compute_effective_Gamma_dGamma(ad_idx3, state.Gamma);

    Real rho = medium.rho(phi, state.theta, state.r);
    Real m2 = medium.mass(phi, state.theta, state.r);

    Real dm2_dt = state.r * state.r * rho * diff.r;

    Real deps_dt = 0;

    if constexpr (State::energy_inject) {
        deps_dt = ejecta.deps_dt(phi, state.theta, t);
    }

    Real a = (state.Gamma - 1) * con::c2 * dm2_dt + (state.Gamma - Gamma4) * con::c2 * diff.m3 +
             Gamma_eff2 * diff.U2_th + Gamma_eff3 * diff.U3_th - deps_dt;
    Real b = (m2 + state.m3) * con::c2 + dGamma_eff2_dGamma * state.U2_th + dGamma_eff3_dGamma * state.U3_th;

    return -a / b;
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_dU2_dt(State const& state, State const& diff, Real t) const noexcept {
    Real e_th = (state.Gamma - 1) * 4 * state.Gamma * medium.rho(phi, state.theta, state.r) * con::c2;
    Real eps_rad = compute_radiative_efficiency(state.t_comv, state.Gamma, e_th, rad_fwd);

    Real ad_idx = adiabatic_idx(state.Gamma);

    Real rho = medium.rho(phi, state.theta, state.r);
    Real dm2_dt = rho * state.r * state.r * diff.r;

    Real shock_heating = compute_shock_heating_rate(state.Gamma, dm2_dt);

    Real adiabatic_cooling = compute_adiabatic_cooling_rate2(ad_idx, state.r, state.x4, state.U2_th, diff.r, diff.x4);

    return (1 - eps_rad) * shock_heating + adiabatic_cooling;
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_dU3_dt(State const& state, State const& diff, Real t) const noexcept {
    Real Gamma34 = compute_rel_Gamma(this->Gamma4, state.Gamma);
    Real ad_idx = adiabatic_idx(Gamma34);
    Real adiabatic_cooling = compute_adiabatic_cooling_rate2(ad_idx, state.r, state.x3, state.U3_th, diff.r, diff.x3);

    Real shock_heating = compute_shock_heating_rate(Gamma34, diff.m3);
    if (is_crossing(state, diff)) {
        Real eps_rad = 0;  ////compute_radiative_efficiency(state.t_comv, state.Gamma, e_th, rad_rvs);
        return (1 - eps_rad) * shock_heating + adiabatic_cooling;
    } else {
        return adiabatic_cooling;
    }
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_dx3_dt(State const& state, State const& diff, Real t) const noexcept {
    if (is_crossing(state, diff)) {
        if (state.Gamma == this->Gamma4) {
            return 0.;
        } else {
            Real sigma = compute_shell_sigma(state);
            Real Gamma34 = compute_rel_Gamma(this->Gamma4, state.Gamma);
            Real beta3 = gamma_to_beta(state.Gamma);
            Real beta4 = gamma_to_beta(this->Gamma4);
            Real comp_ratio = compute_4vel_jump(Gamma34, sigma);
            Real dx3dt = (beta4 - beta3) * con::c / ((1 - beta3) * (state.Gamma * comp_ratio / this->Gamma4 - 1));

            return dx3dt * state.Gamma;
        }
    } else {
        Real Gamma34 = compute_rel_Gamma(this->Gamma4, state.Gamma);
        return compute_shell_spreading_rate(Gamma34, diff.t_comv);
    }
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_dx4_dt(State const& state, State const& diff, Real t) const noexcept {
    if (is_injecting(diff)) {
        return this->u4;
    } else {
        return compute_shell_spreading_rate(this->Gamma4, diff.t_comv);
    }
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_dm3_dt(State const& state, State const& diff, Real t) const noexcept {
    if (is_crossing(state, diff)) {
        if (state.Gamma == this->Gamma4) {
            return 0.;
        } else {
            Real sigma = compute_shell_sigma(state);
            Real Gamma34 = compute_rel_Gamma(this->Gamma4, state.Gamma);
            Real comp_ratio = compute_4vel_jump(Gamma34, sigma);
            Real column_den3 = state.m4 * comp_ratio / state.x4;
            Real dm3dt = column_den3 * diff.x3;

            if (state.m3 >= state.m4) {
                return std::min(dm3dt, diff.m4);
            } else {
                return dm3dt;
            }
        }
    } else {
        return 0.;
    }
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_deps4_dt(State const& state, State const& diff, Real t) const noexcept {
    Real deps4_dt = 0;

    if (t < ejecta.T0) {
        deps4_dt = deps0_dt;
    }

    if constexpr (State::energy_inject) {
        deps4_dt += ejecta.deps_dt(phi, theta0, t);
    }

    return deps4_dt;
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_dm4_dt(State const& state, State const& diff, Real t) const noexcept {
    Real dm4_dt = 0;

    if (t < ejecta.T0) {
        dm4_dt = dm0_dt;
    }

    if constexpr (State::mass_inject) {
        dm4_dt += ejecta.dm_dt(phi, theta0, t);
    }

    return dm4_dt;
}

template <typename Ejecta, typename Medium>
void FRShockEqn<Ejecta, Medium>::operator()(State const& state, State& diff, Real t) {
    Real beta3 = gamma_to_beta(state.Gamma);
    diff.r = compute_dr_dt(beta3);
    diff.t_comv = compute_dt_dt_comv(state.Gamma, beta3);

    diff.eps4 = compute_deps4_dt(state, diff, t);
    diff.m4 = compute_dm4_dt(state, diff, t);

    diff.x4 = compute_dx4_dt(state, diff, t);
    diff.x3 = compute_dx3_dt(state, diff, t);

    diff.m3 = compute_dm3_dt(state, diff, t);

    diff.U2_th = compute_dU2_dt(state, diff, t);
    diff.U3_th = compute_dU3_dt(state, diff, t);

    diff.Gamma = compute_dGamma_dt(state, diff, t);

    diff.theta = 0;
}

inline Real compute_init_comv_shell_width(Real Gamma4, Real t0, Real T);

template <typename Ejecta, typename Medium>
void FRShockEqn<Ejecta, Medium>::save_cross_state(State const& state) {
    r_x = state.r;
    u_x = std::sqrt(state.Gamma * state.Gamma - 1);

    V3_comv_x = r_x * r_x * state.x3;

    Real sigma4 = compute_shell_sigma(state);
    Real comp_ratio34 = compute_compression(Gamma4, state.Gamma, sigma4);
    Real rho4 = state.m4 / (state.r * state.r * state.x4);
    rho3_x = rho4 * comp_ratio34;

    Real B4 = compute_upstr_B(rho4, sigma4);
    B3_ordered_x = B4 * comp_ratio34;
}

template <typename Ejecta, typename Medium>
void FRShockEqn<Ejecta, Medium>::set_init_state(State& state, Real t0) const noexcept {
    Real beta4 = gamma_to_beta(Gamma4);

    state.r = beta4 * con::c * t0 / (1 - beta4);
    state.t_comv = state.r / std::sqrt(Gamma4 * Gamma4 - 1) / con::c;
    state.theta = theta0;

    state.x4 = compute_init_comv_shell_width(Gamma4, t0, ejecta.T0);
    state.x3 = 0;

    Real dt = std::min(t0, ejecta.T0);
    state.eps4 = deps0_dt * dt;
    state.m4 = dm0_dt * dt;
    state.m3 = 0;
    state.Gamma = Gamma4;
    state.U2_th = (state.Gamma - 1) * medium.mass(phi, theta0, state.r) * con::c2;
    state.U3_th = 0;
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the power-law index for post-crossing four-velocity evolution.
 * @details The index transitions from g_low=1.5 for low relative Lorentz factors to g_high=3.5
 *          for high relative Lorentz factors (Blandford-McKee limit).
 * @param gamma_rel Relative Lorentz factor
 * @param k Medium power law index (default: 0)
 * @return The power-law index for velocity evolution
 * <!-- ************************************************************************************** -->
 */
inline Real get_post_cross_g(Real gamma_rel, Real k = 0) {
    constexpr Real g_low = 1.5;   // k is the medium power law index
    constexpr Real g_high = 3.5;  // Blandford-McKee limit// TODO: need to be modified for non ISM medium
    Real p = std::sqrt(std::sqrt(gamma_rel - 1));
    return g_low + (g_high - g_low) * p / (1 + p);
}

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::compute_shell_sigma(State const& state) const {
    return std::max(0.0, state.eps4 / (Gamma4 * state.m4 * con::c2) - 1);
}

//---------------------------------------------------------------------------------------------------------------------
// Helper functions
//---------------------------------------------------------------------------------------------------------------------

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the comoving shell width at initial radius.
 * @details Accounts for both pure injection phase and shell spreading phase.
 * @param Gamma4 Lorentz factor of the unshocked ejecta
 * @param t0 Initial time
 * @param T Engine duration
 * @return The comoving shell width
 * <!-- ************************************************************************************** -->
 */
inline Real compute_init_comv_shell_width(Real Gamma4, Real t0, Real T) {
    Real beta4 = gamma_to_beta(Gamma4);
    if (t0 < T) {  // pure injection
        return Gamma4 * t0 * beta4 * con::c;
    } else {  // injection+shell spreading
        Real cs = compute_sound_speed(Gamma4);
        return Gamma4 * T * beta4 * con::c + cs * (t0 - T) * Gamma4;
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Saves the state of both forward and reverse shocks at a grid point.
 * @details Updates shock properties for both shocks and checks if crossing is complete.
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param k Grid index for time
 * @param eqn The reverse shock equation system
 * @param state Current state of the system
 * @param shock_fwd Forward shock object to update
 * @param shock_rvs Reverse shock object to update
 * @return True if the shock has crossed and there's no more injection
 * <!-- ************************************************************************************** -->
 */
template <typename Eqn, typename State>
bool save_shock_pair_state(size_t i, size_t j, int k, Eqn const& eqn, State const& state, Shock& shock_fwd,
                           Shock& shock_rvs) {
    constexpr Real gamma1 = 1;  // Lorentz factor of unshocked medium (at rest)
    constexpr Real sigma1 = 0;  // Magnetization of unshocked medium

    Real comp_ratio12 = compute_compression(gamma1, state.Gamma, sigma1);
    Real rho1 = eqn.medium.rho(eqn.phi, state.theta, state.r);

    Real m2 = compute_swept_mass(eqn, state);
    Real Gamma2_th = state.U2_th / (m2 * con::c2) + 1;

    Real B1 = compute_upstr_B(rho1, sigma1);
    Real B2 = compute_downstr_B(shock_fwd.rad.eps_B, rho1, B1, Gamma2_th, comp_ratio12);

    save_shock_state(shock_fwd, i, j, k, state.t_comv, state.r, state.theta, state.Gamma, Gamma2_th, B2, m2);

    //------------------------------------------------------------------------------------------------//

    if (k <= shock_rvs.injection_idx(i, j)) {
        Real Gamma4 = eqn.Gamma4;
        Real sigma4 = eqn.compute_shell_sigma(state);

        Real comp_ratio34 = compute_compression(Gamma4, state.Gamma, sigma4);
        Real rho4 = state.m4 / (state.r * state.r * state.x4);

        Real Gamma3_th = state.U3_th / (state.m3 * con::c2) + 1;

        Real B4 = compute_upstr_B(rho4, sigma4);
        Real B3 = compute_downstr_B(shock_rvs.rad.eps_B, rho4, B4, Gamma3_th, comp_ratio34);

        save_shock_state(shock_rvs, i, j, k, state.t_comv, state.r, state.theta, state.Gamma, Gamma3_th, B3, state.m3);
    } else {
        Real V3_comv = state.r * state.r * state.x3;
        Real comp_ratio = eqn.V3_comv_x / V3_comv;

        Real Gamma3_th = state.U3_th / (state.m3 * con::c2) + 1;

        Real B3 = compute_downstr_B(shock_rvs.rad.eps_B, eqn.rho3_x, eqn.B3_ordered_x, Gamma3_th, comp_ratio);

        save_shock_state(shock_rvs, i, j, k, state.t_comv, state.r, state.theta, state.Gamma, Gamma3_th, B3, state.m3);
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Solves the reverse/forward shock ODE at a grid point.
 * @details Manages the evolution of both shocks before and after crossing.
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param t View of time points
 * @param shock_fwd Forward shock object
 * @param shock_rvs Reverse shock object
 * @param eqn Reverse shock equation system
 * @param rtol Relative tolerance for ODE solver
 * <!-- ************************************************************************************** -->
 */
template <typename Eqn, typename View>
void grid_solve_shock_pair(size_t i, size_t j, View const& t, Shock& shock_fwd, Shock& shock_rvs, Eqn& eqn,
                           Real rtol = 1e-6) {
    using namespace boost::numeric::odeint;

    typename Eqn::State state;
    Real t_dec = compute_dec_time(eqn, t.back());
    // Real t0 = min(t_min, t_dec / 100, 1 * unit::sec);
    Real t0 = std::min(1e-1 * unit::sec, t_dec / 100);

    eqn.set_init_state(state, t0);

    if (state.Gamma <= con::Gamma_cut) {
        set_stopping_shock(i, j, shock_fwd, state);
        set_stopping_shock(i, j, shock_rvs, state);
        return;
    }

    auto stepper = make_dense_output(rtol, rtol, runge_kutta_dopri5<typename Eqn::State>());
    stepper.initialize(state, t0, 1e-9 * t0);

    size_t k = 0;
    for (; t(k) < t0; k++) {
        eqn.set_init_state(state, t(k));
        save_shock_pair_state(i, j, k, eqn, state, shock_fwd, shock_rvs);
    }

    bool reverse_shock_crossing = true;

    for (; stepper.current_time() <= t.back();) {
        stepper.do_step(eqn);
        while (k < t.size() && stepper.current_time() > t(k)) {
            stepper.calc_state(t(k), state);
            if (reverse_shock_crossing && !is_crossing(eqn, state, t(k))) {
                shock_rvs.injection_idx(i, j) = k;
                reverse_shock_crossing = false;
                eqn.save_cross_state(state);
            }
            save_shock_pair_state(i, j, k, eqn, state, shock_fwd, shock_rvs);
            ++k;
        }
    }
}

template <typename Ejecta, typename Medium>
ShockPair generate_shock_pair(Coord const& coord, Medium const& medium, Ejecta const& jet, RadParams const& rad_fwd,
                              RadParams const& rad_rvs, Real rtol) {
    auto [phi_size, theta_size, t_size] = coord.shape();
    size_t phi_size_needed = coord.t.shape()[0];
    Shock f_shock(phi_size_needed, theta_size, t_size, rad_fwd);
    Shock r_shock(phi_size_needed, theta_size, t_size, rad_rvs);
    for (size_t i = 0; i < phi_size_needed; ++i) {
        Real theta_s =
            jet_spreading_edge(jet, medium, coord.phi(i), coord.theta.front(), coord.theta.back(), coord.t.front());
        for (size_t j = 0; j < theta_size; ++j) {
            auto eqn_r = FRShockEqn(medium, jet, coord.phi(i), coord.theta(j), rad_fwd, rad_rvs);
            // Solve the forward-reverse shock shell
            grid_solve_shock_pair(i, j, xt::view(coord.t, i, j, xt::all()), f_shock, r_shock, eqn_r, rtol);
        }
    }
    return std::make_pair(std::move(f_shock), std::move(r_shock));
}
