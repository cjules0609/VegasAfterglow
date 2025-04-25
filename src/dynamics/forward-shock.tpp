//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "forward-shock.hpp"
#include "simple-shock.hpp"
/********************************************************************************************************************
 * METHOD: ForwardShockEqn::operator()(State const& y, State& dydr, Real t)
 * DESCRIPTION: Computes the derivatives of the state variables with respect to radius t.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
void ForwardShockEqn<Ejecta, Medium>::operator()(StateArray const& y, StateArray& dydt, Real t) const noexcept {
    FState const state(y);
    FState const diff(dydt);

    Real ad_idx = adiabaticIndex(state.Gamma);
    Real rho = medium.rho(phi, state.theta, state.r);
    Real beta = gammaTobeta(state.Gamma);

    diff.r = drdt(beta);

    if (ejecta.spreading && state.theta < 0.5 * con::pi) {
        diff.theta = dtheta_dt(theta_s, state.theta, diff.r, state.r, state.Gamma);
    } else {
        diff.theta = 0;
    }

    diff.t_com = dtdt_CoMoving(state.Gamma);
    diff.M_sw = state.r * state.r * rho * diff.r;
    diff.M_ej = ejecta.dMdtdOmega(phi, theta0, t);
    diff.E_ej = ejecta.dEdtdOmega(phi, theta0, t);
    diff.Gamma = dGammadt(t, state, diff, ad_idx);
    diff.u = dUdt(state, diff, ad_idx);
}

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
      theta_s(theta_s) {}

/********************************************************************************************************************
 * METHOD: ForwardShockEqn::dGammadt
 * DESCRIPTION: dGammadt with respect to on-axis observe t.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
Real ForwardShockEqn<Ejecta, Medium>::dGammadt(Real t, constState const& state, State const& diff,
                                               Real ad_idx) const noexcept {
    Real Gamma2 = state.Gamma * state.Gamma;
    Real Gamma_eff = (ad_idx * (Gamma2 - 1) + 1) / state.Gamma;
    Real dGamma_eff = (ad_idx * (Gamma2 + 1) - 1) / Gamma2;
    Real dmdt = diff.M_sw;
    Real dlnVdt = 3 / state.r * diff.r;  // only r term
    Real M_sw = state.M_sw;              // Mass per unit solid angle from medium
    Real u = state.u;                    // Internal energy per unit solid angle

    if (ejecta.spreading) {
        Real cos_theta = std::cos(state.theta);
        Real sin_theta = std::sin(state.theta);
        Real f_spread = (1 - cos_theta) / dOmega0;
        dmdt = dmdt * f_spread + M_sw / dOmega0 * sin_theta * diff.theta;
        M_sw *= f_spread;
        dlnVdt += sin_theta / (1 - cos_theta) * diff.theta;
        u *= f_spread;
    }

    Real a1 = -(state.Gamma - 1) * (Gamma_eff + 1) * con::c2 * dmdt;
    Real a2 = (ad_idx - 1) * Gamma_eff * u * dlnVdt;
    Real a3 = diff.E_ej - state.Gamma * diff.M_ej * con::c2;

    Real b1 = (state.M_ej + M_sw) * con::c2;
    Real b2 = (dGamma_eff + Gamma_eff * (ad_idx - 1) / state.Gamma) * u;

    return (a1 + a2 + a3) / (b1 + b2);
}

/********************************************************************************************************************
 * METHOD: ForwardShockEqn::dUdt
 * DESCRIPTION: Computes the derivative of u with respect to time t.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
Real ForwardShockEqn<Ejecta, Medium>::dUdt(constState const& state, State const& diff, Real ad_idx) const noexcept {
    Real dmdt = diff.M_sw;
    Real dlnVdt = 3 / state.r * diff.r - diff.Gamma / state.Gamma;
    if (ejecta.spreading) {
        Real factor = std::sin(state.theta) / (1 - std::cos(state.theta)) * diff.theta;
        dmdt = dmdt + state.M_sw * factor;
        dlnVdt += factor;
        dlnVdt += factor / (ad_idx - 1);
    }

    return (1 - eps_e) * (state.Gamma - 1) * con::c2 * dmdt - (ad_idx - 1) * dlnVdt * state.u;
}

/********************************************************************************************************************
 * FUNCTION: updateForwardShock
 * DESCRIPTION: Updates the forward shock state at grid index (i, j, k) using the current ODE solution at t.
 *              This function computes the physical quantities needed to describe the shock state and updates
 *              the shock object accordingly.
 ********************************************************************************************************************/
template <typename Eqn, typename State>
void updateForwardShock(size_t i, size_t j, size_t k, Eqn const& eqn, State const& state, Shock& shock) {
    // Calculate number density of the ambient medium at current position
    Real n1 = eqn.medium.rho(eqn.phi, state.theta, state.r) / con::mp;

    // Calculate number of protons per unit solid angle in the swept-up material
    Real N2 = state.M_sw / con::mp;

    // Set constant parameters for the unshocked medium
    constexpr Real gamma1 = 1;  // Lorentz factor of unshocked medium (at rest)
    constexpr Real sigma1 = 0;  // Magnetization of unshocked medium

    // Update the shock state with calculated values
    updateShockState(shock, i, j, k, state, state.Gamma, gamma1, N2, n1, sigma1);
}

/********************************************************************************************************************
 * FUNCTION: setForwardInit
 * DESCRIPTION: Set the initial conditions for the forward shock ODE solver.
 *              This function computes the initial state of the shock based on the ejecta properties
 *              and the ambient medium.
 * RETURNS: The deceleration time, which helps determine an appropriate time step.
 ********************************************************************************************************************/
template <typename Eqn, typename State>
void setForwardInit(Eqn const& eqn, State& state, Real t0) {
    // Set initial Lorentz factor from ejecta model
    state.Gamma = eqn.ejecta.Gamma0(eqn.phi, eqn.theta0);

    // Calculate initial radius based on observer time t0
    Real beta0 = gammaTobeta(state.Gamma);
    state.r = beta0 * con::c * t0 / (1 - beta0);

    // Get ambient medium density at initial position
    Real rho = eqn.medium.rho(eqn.phi, eqn.theta0, state.r);

    // Calculate initial swept-up mass per solid angle (assuming spherical expansion)
    state.M_sw = 1. / 3 * rho * state.r * state.r * state.r;

    // Set initial ejecta energy per solid angle from model
    state.E_ej = eqn.ejecta.dE0dOmega(eqn.phi, eqn.theta0);

    // Calculate initial ejecta mass per solid angle
    state.M_ej = state.E_ej / (state.Gamma * (1 + eqn.ejecta.sigma0(eqn.phi, eqn.theta0)) * con::c2);

    // Set initial internal energy
    state.u = (state.Gamma - 1) * state.M_sw * con::c2;

    // Calculate initial comoving time
    state.t_com = state.r / std::sqrt(state.Gamma * state.Gamma - 1) / con::c;

    // Set initial theta to the input value
    state.theta = eqn.theta0;
}
/********************************************************************************************************************
 * FUNCTION: solveForwardShell
 * DESCRIPTION: Solve the forward shock ODE at grid (phi[i], theta[j]) as a function of t (on-axis observation time).
 *              This function uses an adaptive step size ODE solver to evolve the shock state through time.
 ********************************************************************************************************************/
template <typename FwdEqn, typename View>
void solveForwardShell(size_t i, size_t j, View const& t, Shock& shock, FwdEqn const& eqn, double rtol) {
    using namespace boost::numeric::odeint;

    // Initialize state array and wrapper
    typename FwdEqn::StateArray y{};
    FState const state(y);

    // Get initial time and set up initial conditions
    Real t0 = std::min(t.front(), 1 * con::sec);

    setForwardInit(eqn, state, t0);
    Real dt = t0 / 10;  // Initial time step based on deceleration time

    // Early exit if initial Lorentz factor is below cutoff
    if (state.Gamma <= con::Gamma_cut) {
        setStoppingShock(i, j, shock, state);
        return;
    }

    // Set up ODE solver with adaptive step size control
    auto stepper = make_dense_output(rtol, rtol, runge_kutta_dopri5<typename FwdEqn::StateArray>());
    stepper.initialize(y, t0, dt);

    // Solve ODE and update shock state at each requested time point
    for (size_t k = 0; stepper.current_time() <= t.back();) {
        // Advance solution by one adaptive step
        stepper.do_step(eqn);

        // Update shock state for all time points that have been passed in this step
        while (k < t.size() && stepper.current_time() > t(k)) {
            stepper.calc_state(t(k), y);
            updateForwardShock(i, j, k, eqn, state, shock);
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
Shock genForwardShock(Coord const& coord, Medium const& medium, Ejecta const& jet, Real eps_e, Real eps_B,
                      Real rtol = 1e-6) {
    auto [phi_size, theta_size, t_size] = coord.shape();  // Unpack coordinate dimensions
    size_t phi_size_needed = coord.t.shape()[0];
    Shock shock(phi_size_needed, theta_size, t_size, eps_e, eps_B);

    for (size_t i = 0; i < phi_size_needed; ++i) {
        Real theta_s =
            jetSpreadingEdge(jet, medium, coord.phi(i), coord.theta.front(), coord.theta.back(), coord.t.front());
        for (size_t j = 0; j < theta_size; ++j) {
            // Create a ForwardShockEqn for each theta slice
            // auto eqn = ForwardShockEqn(medium, jet, coord.phi[i], coord.theta[j], eps_e, theta_s);
            auto eqn = SimpleShockEqn(medium, jet, coord.phi(i), 0, coord.theta(j), eps_e, theta_s);
            //   Solve the shock shell for this theta slice
            solveForwardShell(i, j, xt::view(coord.t, i, j, xt::all()), shock, eqn, rtol);
        }
    }

    return shock;
}

template <typename Ejecta, typename Medium>
void genForwardShock(Shock& shock, Coord const& coord, Medium const& medium, Ejecta const& jet, Real eps_e, Real eps_B,
                     Real rtol = 1e-6) {
    auto [phi_size, theta_size, t_size] = coord.shape();  // Unpack coordinate dimensions
    size_t phi_size_needed = coord.t.shape()[0];
    shock.resize(phi_size_needed, theta_size, t_size);
    shock.eps_B = eps_B;
    shock.eps_e = eps_e;
    for (size_t i = 0; i < phi_size_needed; ++i) {
        Real theta_s =
            jetSpreadingEdge(jet, medium, coord.phi(i), coord.theta.front(), coord.theta.back(), coord.t.front());
        for (size_t j = 0; j < theta_size; ++j) {
            // Create a ForwardShockEqn for each theta slice
            // auto eqn = ForwardShockEqn(medium, jet, coord.phi(i), coord.theta(j), eps_e, theta_s);
            auto eqn = SimpleShockEqn(medium, jet, coord.phi(i), 0, coord.theta(j), eps_e, theta_s);
            //   Solve the shock shell for this theta slice
            solveForwardShell(i, j, xt::view(coord.t, i, j, xt::all()), shock, eqn, rtol);
        }
    }
}