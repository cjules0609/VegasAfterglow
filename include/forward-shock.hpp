//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _FORWARDSHOCK_
#define _FORWARDSHOCK_
#include "jet.h"
#include "shock.h"
/********************************************************************************************************************
 * CLASS: FState
 * DESCRIPTION: This is a helper class to provide named access to the forward ODE state array (required by ODE solver)
 *              components.
 ********************************************************************************************************************/
template <typename StateArray>
struct FState {
    using T = decltype(std::declval<StateArray>()[0]);

    FState() = delete;

    FState(StateArray& y)
        : Gamma(y[0]), u(y[1]), r(y[2]), t_com(y[3]), theta(y[4]), M_sw(y[5]), M_ej(y[6]), E_ej(y[7]) {}

    T& Gamma;
    T& u;
    T& r;
    T& t_com;
    T& theta;
    T& M_sw;
    T& M_ej;
    T& E_ej;
};

/********************************************************************************************************************
 * CLASS: ForwardShockEqn
 * DESCRIPTION: Represents the forward shock equation for a given jet and medium. It defines a state vector
 *              (an array of 8 Reals) and overloads operator() to compute the derivatives of the state with
 *              respect to t. It also declares helper functions for the derivatives. Comprehensive models
 ********************************************************************************************************************/
class ForwardShockEqn {
   public:
    // State vector: [Gamma, u, r, t_com, theta, m, M_ej, E_ej]
    // - Gamma: Bulk Lorentz factor
    // - u: internal energy per solid angle
    // - r: radius
    // - t_com: comoving time
    // - theta: jet opening angle
    // - M_sw: swept mass per solid angle
    // - M_ej: ejecta mass per solid angle
    // - E_ej: ejecta energy per solid angle
    using StateArray = std::array<Real, 8>;
    using State = FState<StateArray>;
    using constState = FState<const StateArray>;

    ForwardShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, Real eps_e, Real theta_s);

    Medium const& medium;  // Reference to the medium properties
    Ejecta const& ejecta;  // Reference to the ejecta properties
    Real const phi{0};     // Angular coordinate phi
    Real const theta0{0};  // Angular coordinate theta
    Real const eps_e{0};   // Electron energy fraction

    // forward shock ODE equation
    void operator()(StateArray const& y, StateArray& dydt, Real t);

   private:
    // dGammadt with respect to on-axis observe t.
    inline Real dGammadt(Real t, constState const& state, State const& diff, Real ad_idx);
    // dUdt with respect to on-axis observe t.
    inline Real dUdt(constState const& state, State const& diff, Real ad_idx);

    Real const dOmega0{0};  // Initial solid angle
    Real const theta_s{0};
};

/********************************************************************************************************************
 * FUNCTION: updateForwardShock
 * DESCRIPTION: Updates the forward shock state at grid index (i, j, k) using the current ODE solution at t.
 ********************************************************************************************************************/
template <typename Eqn, typename State>
void updateForwardShock(size_t i, size_t j, int k, Eqn const& eqn, State const& state, Shock& shock) {
    Real n1 = eqn.medium.rho(eqn.phi, state.theta, state.r) / con::mp;
    Real N2 = state.M_sw / con::mp;  // number of proton per unit solid angle
    constexpr Real gamma1 = 1;
    constexpr Real sigma1 = 0;
    updateShockState(shock, i, j, k, state, state.Gamma, gamma1, N2, n1, sigma1);
}

/********************************************************************************************************************
 * FUNCTION: setForwardInit
 * DESCRIPTION: Set the initial conditions for the forward shock ODE solver.
 ********************************************************************************************************************/
template <typename Eqn, typename State>
Real setForwardInit(Eqn const& eqn, State& state, Real t0) {
    state.Gamma = eqn.ejecta.Gamma0(eqn.phi, eqn.theta0);  // Initial Lorentz factor
    Real beta0 = gammaTobeta(state.Gamma);
    state.r = beta0 * con::c * t0 / (1 - beta0);
    Real rho = eqn.medium.rho(eqn.phi, eqn.theta0, state.r);
    state.M_sw = 1. / 3 * rho * state.r * state.r * state.r;  // swept mass per solid angle
    state.E_ej = eqn.ejecta.dE0dOmega(eqn.phi, eqn.theta0);   //
    state.M_ej = state.E_ej / (state.Gamma * (1 + eqn.ejecta.sigma0(eqn.phi, eqn.theta0)) * con::c2);
    state.u = (state.Gamma - 1) * state.M_sw * con::c2;
    state.t_com = state.r / std::sqrt(state.Gamma * state.Gamma - 1) / con::c;
    state.theta = eqn.theta0;
    Real n_ism = eqn.medium.rho(eqn.phi, eqn.theta0, state.r) / con::mp;

    Real r_dec = thinShellDecRadius(state.E_ej * 4 * con::pi, n_ism, state.Gamma);
    Real t_dec = r_dec / (2 * state.Gamma * state.Gamma * con::c);
    return t_dec;
}

/********************************************************************************************************************
 * FUNCTION: solveForwardShell
 * DESCRIPTION: Solve the forward shock ODE at grid (phi[i], theta[j]) as a function of t (on-axis observation time).
 ********************************************************************************************************************/
// Solves the forward shock evolution for a given shell (across radius values in array r) and updates the Shock object.
template <typename FwdEqn>
void solveForwardShell(size_t i, size_t j, const Array& t, Shock& shock, FwdEqn const& eqn, double rtol) {
    using namespace boost::numeric::odeint;

    typename FwdEqn::StateArray y;
    FState state(y);

    Real t0 = t[0];
    Real t_dec = setForwardInit(eqn, state, t0);  // Initialize state at starting radius
    Real dt = t_dec / 100;

    if (state.Gamma <= con::Gamma_cut) {                         // If initial Lorentz factor is too low, exit early
        setStoppingShock(i, j, shock, t, state.r, state.theta);  // Set the shock state to zero
        return;
    }

    // ODE solver with adaptive step size, and relative tolerance rtol
    // auto stepper = bulirsch_stoer_dense_out<typename FwdEqn::StateArray>{0, rtol};
    auto stepper = make_dense_output(0, rtol, runge_kutta_dopri5<typename FwdEqn::StateArray>());
    stepper.initialize(y, t0, dt);
    Real t_back = t[t.size() - 1];  // Last time in the array

    // Iterate over the radius array, updating the state and the Shock object as needed
    for (int k = 0; stepper.current_time() <= t_back;) {
        stepper.do_step(eqn);
        while (k < t.size() && stepper.current_time() > t[k]) {
            stepper.calc_state(t[k], y);
            updateForwardShock(i, j, k, eqn, state, shock);
            ++k;
        }
    }
}
#endif