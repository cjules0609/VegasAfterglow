//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include "jet.h"
#include "shock.h"
#include <array>
#include <cmath>

/********************************************************************************************************************
 * CLASS: FState
 * DESCRIPTION: This is a helper class to provide named access to the forward ODE state array (required by ODE solver)
 *              components.
 ********************************************************************************************************************/
template <typename StateArray>
struct FState {
    using T = decltype(std::declval<StateArray>()[0]);

    // Prevent default construction to ensure proper initialization
    FState() = delete;

    // Constructor that maps array elements to named references
    FState(StateArray& y)
        : Gamma(y[0]), u(y[1]), r(y[2]), t_com(y[3]), theta(y[4]), M_sw(y[5]), M_ej(y[6]), E_ej(y[7]) {}

    // References to state variables
    T& Gamma;   // Bulk Lorentz factor
    T& u;       // Internal energy per solid angle
    T& r;       // Radius
    T& t_com;   // Comoving time
    T& theta;   // Jet opening angle
    T& M_sw;    // Swept mass per solid angle
    T& M_ej;    // Ejecta mass per solid angle
    T& E_ej;    // Ejecta energy per solid angle
};

/********************************************************************************************************************
 * CLASS: ForwardShockEqn
 * DESCRIPTION: Represents the forward shock equation for a given jet and medium. It defines a state vector
 *              (an array of 8 Reals) and overloads operator() to compute the derivatives of the state with
 *              respect to t. It also declares helper functions for the derivatives.
 ********************************************************************************************************************/
class ForwardShockEqn {
   public:
    // State vector: [Gamma, u, r, t_com, theta, M_sw, M_ej, E_ej]
    using StateArray = std::array<Real, 8>;
    using State = FState<StateArray>;
    using constState = FState<const StateArray>;

    // Constructor
    ForwardShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, Real eps_e, Real theta_s);

    // References to model components
    Medium const& medium;  // Reference to the medium properties
    Ejecta const& ejecta;  // Reference to the ejecta properties
    
    // Model parameters
    Real const phi{0};     // Angular coordinate phi
    Real const theta0{0};  // Angular coordinate theta
    Real const eps_e{0};   // Electron energy fraction

    // Forward shock ODE equation - callable interface for ODE solver
    void operator()(StateArray const& y, StateArray& dydt, Real t);

   private:
    
    // dGammadt with respect to on-axis observer time
    inline Real dGammadt(Real t, constState const& state, State const& diff, Real ad_idx);
    
    // dUdt with respect to on-axis observer time
    inline Real dUdt(constState const& state, State const& diff, Real ad_idx);

    // Private member variables
    Real const dOmega0{0};  // Initial solid angle
    Real const theta_s{0};  // Jet structure parameter
};

/********************************************************************************************************************
 * FUNCTION: updateForwardShock
 * DESCRIPTION: Updates the forward shock state at grid index (i, j, k) using the current ODE solution at t.
 ********************************************************************************************************************/
template <typename Eqn, typename State>
void updateForwardShock(size_t i, size_t j, int k, Eqn const& eqn, State const& state, Shock& shock) {
    // Calculate number density of the ambient medium
    Real n1 = eqn.medium.rho(eqn.phi, state.theta, state.r) / con::mp;
    
    // Calculate number of protons per unit solid angle in the swept-up material
    Real N2 = state.M_sw / con::mp;
    
    // Set constant parameters for the unshocked medium
    constexpr Real gamma1 = 1;     // Lorentz factor of unshocked medium (at rest)
    constexpr Real sigma1 = 0;     // Magnetization of unshocked medium
    
    // Update the shock state with calculated values
    updateShockState(shock, i, j, k, state, state.Gamma, gamma1, N2, n1, sigma1);
}

/********************************************************************************************************************
 * FUNCTION: setForwardInit
 * DESCRIPTION: Set the initial conditions for the forward shock ODE solver.
 * RETURNS: The deceleration time, which helps determine an appropriate time step.
 ********************************************************************************************************************/
template <typename Eqn, typename State>
Real setForwardInit(Eqn const& eqn, State& state, Real t0) {
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
    
    // Calculate number density for deceleration radius calculation
    Real n_ism = eqn.medium.rho(eqn.phi, eqn.theta0, state.r) / con::mp;

    // Calculate deceleration radius and time
    Real r_dec = thinShellDecRadius(state.E_ej * 4 * con::pi, n_ism, state.Gamma);
    Real t_dec = r_dec / (2 * state.Gamma * state.Gamma * con::c);
    
    return t_dec;
}

/********************************************************************************************************************
 * FUNCTION: solveForwardShell
 * DESCRIPTION: Solve the forward shock ODE at grid (phi[i], theta[j]) as a function of t (on-axis observation time).
 ********************************************************************************************************************/
template <typename FwdEqn>
void solveForwardShell(size_t i, size_t j, const Array& t, Shock& shock, FwdEqn const& eqn, double rtol) {
    using namespace boost::numeric::odeint;

    // Initialize state array and wrapper
    typename FwdEqn::StateArray y;
    FState state(y);

    // Get initial time and set up initial conditions
    Real t0 = t.front();
    Real t_dec = setForwardInit(eqn, state, t0);
    Real dt = t_dec / 100;  // Initial time step based on deceleration time
    
    // Set default injection index (will be updated if shock propagates)
    shock.injection_idx(i, j) = t.size();
    
    // Early exit if initial Lorentz factor is below cutoff
    if (state.Gamma <= con::Gamma_cut) {
        setStoppingShock(i, j, shock, t, state.r, state.theta);
        return;
    }

    // Set up ODE solver with adaptive step size control
    auto stepper = make_dense_output(0, rtol, runge_kutta_dopri5<typename FwdEqn::StateArray>());
    stepper.initialize(y, t0, dt);
    Real t_back = t.back();  // Last time in the array

    // Solve ODE and update shock state at each requested time point
    for (int k = 0; stepper.current_time() <= t_back;) {
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
