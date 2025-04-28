//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "forward-shock.hpp"
#include "simple-shock.hpp"

/********************************************************************************************************************
 * CONSTRUCTOR: ForwardShockEqn::ForwardShockEqn
 * DESCRIPTION: ForwardShockEqn constructor
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
ForwardShockEqn<Ejecta, Medium>::ForwardShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta,
                                                 Real eps_e, Real theta_s)
    : medium(medium),
      ejecta(ejecta),
      phi(phi),
      theta0(theta),
      eps_e(eps_e),
      dOmega0(1 - std::cos(theta)),
      theta_s(theta_s),
      m_shell(0) {
    m_shell = ejecta.eps_k(phi, theta0) / ejecta.Gamma0(phi, theta0) / con::c2;
    if constexpr (HasSigma<Ejecta>) {
        m_shell /= 1 + ejecta.sigma0(phi, theta0);
    }
}

/********************************************************************************************************************
 * METHOD: ForwardShockEqn::operator()(State const& y, State& dydr, Real t)
 * DESCRIPTION: Computes the derivatives of the state variables with respect to engine t.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
void ForwardShockEqn<Ejecta, Medium>::operator()(State const& state, State& diff, Real t) const noexcept {
    Real beta = gamma_to_beta(state.Gamma);

    diff.r = compute_dr_dt(beta);
    diff.t_comv = compute_dt_dt_comv(state.Gamma);

    if (ejecta.spreading && state.theta < 0.5 * con::pi) {
        diff.theta = compute_dtheta_dt(theta_s, state.theta, diff.r, state.r, state.Gamma);
    } else {
        diff.theta = 0;
    }

    if constexpr (State::mass_inject) {
        diff.m_shell = ejecta.dm_dt(phi, theta0, t);
    }

    if constexpr (State::energy_inject) {
        diff.eps_shell = ejecta.deps_dt(phi, theta0, t);
    }

    Real dm_dt_swpt = state.r * state.r * medium.rho(phi, state.theta, state.r) * diff.r;
    Real m_swept = compute_swept_mass(*this, state);

    if constexpr (!State::mass_profile) {
        diff.m_swept = dm_dt_swpt;
    }

    Real ad_idx = adiabatic_idx(state.Gamma);

    diff.Gamma = dGamma_dt(m_swept, dm_dt_swpt, state, diff, ad_idx);
    diff.u = dU_dt(m_swept, dm_dt_swpt, state, diff, ad_idx);
}

/********************************************************************************************************************
 * METHOD: ForwardShockEqn::dGamma_dt
 * DESCRIPTION: dGamma_dt with respect to engine time t.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
Real ForwardShockEqn<Ejecta, Medium>::dGamma_dt(Real m_swept, Real dm_dt_swept, State const& state, State const& diff,
                                                Real ad_idx) const noexcept {
    Real Gamma2 = state.Gamma * state.Gamma;
    Real Gamma_eff = (ad_idx * (Gamma2 - 1) + 1) / state.Gamma;
    Real dGamma_eff = (ad_idx * (Gamma2 + 1) - 1) / Gamma2;
    Real dlnVdt = 3 / state.r * diff.r;  // only r term

    Real m_shell = this->m_shell;
    Real u = state.u;  // Internal energy per unit solid angle

    if (ejecta.spreading) {
        Real cos_theta = std::cos(state.theta);
        Real sin_theta = std::sin(state.theta);
        Real f_spread = (1 - cos_theta) / dOmega0;
        dm_dt_swept = dm_dt_swept * f_spread + m_swept / dOmega0 * sin_theta * diff.theta;
        m_swept *= f_spread;
        dlnVdt += sin_theta / (1 - cos_theta) * diff.theta;
        u *= f_spread;
    }

    Real a1 = -(state.Gamma - 1) * (Gamma_eff + 1) * con::c2 * dm_dt_swept;
    Real a2 = (ad_idx - 1) * Gamma_eff * u * dlnVdt;

    if constexpr (State::energy_inject) {
        a1 += diff.eps_shell;
    }

    if constexpr (State::mass_inject) {
        a1 -= state.Gamma * diff.m_shell * con::c2;
        m_shell = state.m_shell;
    }

    Real b1 = (m_shell + m_swept) * con::c2;
    Real b2 = (dGamma_eff + Gamma_eff * (ad_idx - 1) / state.Gamma) * u;

    return (a1 + a2) / (b1 + b2);
}

/********************************************************************************************************************
 * METHOD: ForwardShockEqn::dUdt
 * DESCRIPTION: Computes the derivative of u with respect to time t.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
Real ForwardShockEqn<Ejecta, Medium>::dU_dt(Real m_swept, Real dm_dt_swept, State const& state, State const& diff,
                                            Real ad_idx) const noexcept {
    Real dlnVdt = 3 / state.r * diff.r - diff.Gamma / state.Gamma;
    if (ejecta.spreading) {
        Real factor = std::sin(state.theta) / (1 - std::cos(state.theta)) * diff.theta;
        dm_dt_swept = dm_dt_swept + m_swept * factor;
        dlnVdt += factor;
        dlnVdt += factor / (ad_idx - 1);
    }

    return (1 - eps_e) * (state.Gamma - 1) * con::c2 * dm_dt_swept - (ad_idx - 1) * dlnVdt * state.u;
}

/********************************************************************************************************************
 * FUNCTION: setForwardInit
 * DESCRIPTION: Set the initial conditions for the forward shock ODE solver.
 *              This function computes the initial state of the shock based on the ejecta properties
 *              and the ambient medium.
 * RETURNS: The deceleration time, which helps determine an appropriate time step.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
void ForwardShockEqn<Ejecta, Medium>::set_init_state(State& state, Real t0) const noexcept {
    state.Gamma = ejecta.Gamma0(phi, theta0);

    Real beta0 = gamma_to_beta(state.Gamma);
    state.r = beta0 * con::c * t0 / (1 - beta0);

    state.t_comv = state.r / std::sqrt(state.Gamma * state.Gamma - 1) / con::c;

    state.theta = theta0;

    if constexpr (State::energy_inject) {
        state.eps_shell = ejecta.eps_k(phi, theta0);
    }

    if constexpr (State::mass_inject) {
        state.m_shell = m_shell;
    }

    if constexpr (!State::mass_profile) {
        state.m_swept = medium.rho(phi, theta0, state.r) * state.r * state.r * state.r / 3;
        state.u = (state.Gamma - 1) * state.m_swept * con::c2;
    } else {
        state.u = (state.Gamma - 1) * medium.mass(phi, theta0, state.r) * con::c2;
    }
}

/********************************************************************************************************************
 * FUNCTION: updateForwardShock
 * DESCRIPTION: Updates the forward shock state at grid index (i, j, k) using the current ODE solution at t.
 *              This function computes the physical quantities needed to describe the shock state and updates
 *              the shock object accordingly.
 ********************************************************************************************************************/
template <typename Eqn, typename State>
void save_fwd_shock_state(size_t i, size_t j, size_t k, Eqn const& eqn, State const& state, Shock& shock) {
    // Calculate number density of the ambient medium at current position
    Real n1 = eqn.medium.rho(eqn.phi, state.theta, state.r) / con::mp;
    Real m2 = compute_swept_mass(eqn, state);

    // Set constant parameters for the unshocked medium
    constexpr Real gamma1 = 1;  // Lorentz factor of unshocked medium (at rest)
    constexpr Real sigma1 = 0;  // Magnetization of unshocked medium

    // Update the shock state with calculated values
    save_shock_state(shock, i, j, k, state, state.Gamma, gamma1, m2 / con::mp, n1, sigma1);
}

/********************************************************************************************************************
 * FUNCTION: solveForwardShell
 * DESCRIPTION: Solve the forward shock ODE at grid (phi[i], theta[j]) as a function of t (on-axis observation time).
 *              This function uses an adaptive step size ODE solver to evolve the shock state through time.
 ********************************************************************************************************************/
template <typename FwdEqn, typename View>
void grid_solve_fwd_shock(size_t i, size_t j, View const& t, Shock& shock, FwdEqn const& eqn, double rtol) {
    using namespace boost::numeric::odeint;

    // Initialize state array
    typename FwdEqn::State state;

    // Get initial time and set up initial conditions
    Real t_dec = compute_dec_time(eqn, t.front(), t.back());
    Real t0 = std::min(t.front(), std::max(t_dec / 2, 10 * con::sec));
    eqn.set_init_state(state, t0);

    // Early exit if initial Lorentz factor is below cutoff
    if (state.Gamma <= con::Gamma_cut) {
        set_stopping_shock(i, j, shock, state);
        return;
    }

    // Set up ODE solver with adaptive step size control
    auto stepper = make_dense_output(rtol, rtol, runge_kutta_dopri5<typename FwdEqn::State>());

    stepper.initialize(state, t0, 0.1 * t0);

    // Solve ODE and update shock state at each requested time point
    for (size_t k = 0; stepper.current_time() <= t.back();) {
        // Advance solution by one adaptive step
        stepper.do_step(eqn);

        // Update shock state for all time points that have been passed in this step
        while (k < t.size() && stepper.current_time() > t(k)) {
            stepper.calc_state(t(k), state);
            save_fwd_shock_state(i, j, k, eqn, state, shock);
            ++k;
        }
    }
}

/********************************************************************************************************************
 * FUNCTION PROTOTYPES: Shock Generation Interfaces
 * DESCRIPTION: These function templates declare interfaces to generate forward shocks (2D and 3D) and
 *              forward/reverse shock pairs.
 ********************************************************************************************************************/
using ShockPair = std::pair<Shock, Shock>;

template <typename Ejecta, typename Medium>
Shock generate_fwd_shock(Coord const& coord, Medium const& medium, Ejecta const& jet, Real eps_e, Real eps_B,
                         Real rtol = 1e-5) {
    auto [phi_size, theta_size, t_size] = coord.shape();  // Unpack coordinate dimensions
    size_t phi_size_needed = coord.t.shape()[0];
    Shock shock(phi_size_needed, theta_size, t_size, eps_e, eps_B);

    for (size_t i = 0; i < phi_size_needed; ++i) {
        Real theta_s =
            jet_spreading_edge(jet, medium, coord.phi(i), coord.theta.front(), coord.theta.back(), coord.t.front());
        for (size_t j = 0; j < theta_size; ++j) {
            // auto eqn = ForwardShockEqn(medium, jet, coord.phi[i], coord.theta[j], eps_e, theta_s);
            auto eqn = SimpleShockEqn(medium, jet, coord.phi(i), coord.theta(j), eps_e, theta_s);
            //  Solve the shock shell for this theta slice
            grid_solve_fwd_shock(i, j, xt::view(coord.t, i, j, xt::all()), shock, eqn, rtol);
        }
    }

    return shock;
}
