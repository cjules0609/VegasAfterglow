//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _FORWARDSHOCK_
#define _FORWARDSHOCK_
#include "shock.h"

/********************************************************************************************************************
 * CLASS: ForwardShockEqn
 * DESCRIPTION: Represents the forward shock equation for a given Jet and Injector. It defines a state vector
 *              (an array of 5 Reals) and overloads operator() to compute the derivatives of the state with
 *              respect to radius r. It also declares helper functions for the derivatives. Comprehensive models
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
class ForwardShockEqn {
   public:
    using StateArray = std::array<Real, 5>;  // State vector: [Gamma, u, r, t_com, theta]

    ForwardShockEqn(Medium const& medium, Jet const& jet, Injector const& inject, Real phi, Real theta, Real eps_e);

    Medium const& medium;     // Reference to the medium properties
    Jet const& jet;           // Reference to the jet properties
    Injector const& inject;   // Reference to the injector properties
    Real const phi{0};        // Angular coordinate phi
    Real const theta0{0};     // Angular coordinate theta
    Real const eps_e{0};      // Electron energy fraction
    Real const jet_sigma{0};  // Jet magnetization parameter
    Real gamma0{1};           // Initial Lorentz factor (or a related parameter)

    // Overloaded operator() to compute the derivatives of the state vector with respect to radius t.
    void operator()(StateArray const& y, StateArray& dydt, Real t);

   private:
    // Helper function: computes the derivative of Gamma with respect to t.
    inline Real dGammadt(Real t, Real Gamma, Real u, Real r, Real theta, Real drdt, Real dthetadt, Real ad_idx,
                         Real rho);
    // Helper function: computes the derivative of u with respect to t.
    inline Real dUdt(Real Gamma, Real u, Real r, Real theta, Real dGdt, Real drdt, Real dthetadt, Real ad_idx,
                     Real rho);

    Real const inj_Gamma0{0};  // Initial Gamma from the injector
    Real const inj_sigma{0};   // Injector magnetization parameter
    Real const dM0{0};         // Initial mass per unit solid angle
    Real const dOmega0{0};     // Initial solid angle
};

// forward shock state with named member access
struct FState {
    FState() = delete;

    template <typename StateArray>
    FState(StateArray& y) : Gamma(y[0]), u(y[1]), r(y[2]), t_com(y[3]), theta(y[4]) {}
    Real& Gamma;
    Real& u;
    Real& r;
    Real& t_com;
    Real& theta;
};

struct constFState {
    constFState() = delete;

    template <typename StateArray>
    constFState(StateArray const& y) : Gamma(y[0]), u(y[1]), r(y[2]), t_com(y[3]), theta(y[4]) {}
    Real const& Gamma;
    Real const& u;
    Real const& r;
    Real const& t_com;
    Real const& theta;
};

/********************************************************************************************************************
 * METHOD: ForwardShockEqn::operator()(State const& y, State& dydr, Real t)
 * DESCRIPTION: Computes the derivatives of the state variables with respect to radius t.
 *              The state vector components are:
 *                  y[0] - Gamma (Lorentz factor)
 *                  y[1] - u (internal energy per solid angle)
 *                  y[2] - r (radius)
 *                  y[3] - t_com (co-moving time) [unused here]
 *                  y[4] - theta_j (jet opening angle)
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
void ForwardShockEqn<Jet, Injector>::operator()(StateArray const& y, StateArray& dydt, Real t) {
    constFState state(y);

    Real ad_idx = adiabaticIndex(state.Gamma);  // Compute adiabatic index based on Gamma
    Real rho = medium.rho(state.r);             // Get medium density at radius r
    Real beta = gammaTobeta(state.Gamma);       // Convert Gamma to beta (velocity/c)
    Real uv = state.Gamma * beta;

    dydt[2] = drdt(beta);  // Compute derivative of engine time with respect to t

    if (jet.spreading && state.theta < 0.5 * con::pi && uv * state.theta < 0.5) {
        dydt[4] = dtheta_dt(uv, dydt[2], state.r, state.Gamma);  // d(theta_jet)/dt
    } else {
        dydt[4] = 0;
    }

    dydt[0] = dGammadt(t, state.Gamma, state.u, state.r, state.theta, dydt[2], dydt[4], ad_idx, rho);    // d(Gamma)/dt
    dydt[1] = dUdt(state.Gamma, state.u, state.r, state.theta, dydt[0], dydt[2], dydt[4], ad_idx, rho);  // d(u)/dt
    dydt[3] = dtdt_CoMoving(state.Gamma);                                                                // d(t_com)/dt
}

/********************************************************************************************************************
 * CONSTRUCTOR: ForwardShockEqn::ForwardShockEqn
 * DESCRIPTION: Initializes a ForwardShockEqn object with references to the medium, jet, and injector,
 *              along with the angular coordinates and energy fraction.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
ForwardShockEqn<Jet, Injector>::ForwardShockEqn(Medium const& medium, Jet const& jet, Injector const& inject, Real phi,
                                                Real theta, Real eps_e)
    : medium(medium),
      jet(jet),
      inject(inject),
      phi(phi),
      theta0(theta),
      eps_e(eps_e),
      jet_sigma(jet.sigma0(phi, theta, 0)),
      gamma0(jet.Gamma0(phi, theta, 0)),
      inj_Gamma0(inject.Gamma0(phi, theta, 0)),
      inj_sigma(inject.sigma0(phi, theta, 0)),
      dM0(jet.dEdOmega(phi, theta, 0) / (gamma0 * (1 + jet_sigma) * con::c2)),
      dOmega0(1 - std::cos(theta)) {}

/********************************************************************************************************************
 * METHOD: ForwardShockEqn::dGammadt
 * DESCRIPTION: Computes the derivative of Gamma with respect to  t.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
Real ForwardShockEqn<Jet, Injector>::dGammadt(Real t, Real Gamma, Real u, Real r, Real theta, Real drdt, Real dthetadt,
                                              Real ad_idx, Real rho) {
    Real Gamma2 = Gamma * Gamma;
    Real Gamma_eff = (ad_idx * (Gamma2 - 1) + 1) / Gamma;
    Real dGamma_eff = (ad_idx * Gamma2 + ad_idx - 1) / Gamma2;
    Real dmdt = r * r * rho * drdt;
    Real dlnVdt = 3 / r * drdt;  // only r term

    Real dm = medium.mass(r) / (4 * con::pi);  // Mass per unit solid angle from medium
    Real dm_inj = inject.dEdOmega(phi, theta0, t) / (inj_Gamma0 * (1 + inj_sigma) * con::c2);  // Injected mass
    Real L_inj = inject.dLdOmega(phi, theta0, t);  // Injected luminosity per unit solid angle

    if (jet.spreading) {
        Real f_spread = (1 - std::cos(theta)) / dOmega0;
        dmdt = dmdt * f_spread + dm / dOmega0 * std::sin(theta) * dthetadt;
        dm *= f_spread;
        dlnVdt += std::sin(theta) / (1 - std::cos(theta)) * dthetadt;
        u *= f_spread;
    }

    Real a1 = -(Gamma - 1) * (Gamma_eff + 1) * con::c2 * dmdt;
    Real a2 = (ad_idx - 1) * Gamma_eff * u * dlnVdt;
    Real a3 = L_inj * (1 - Gamma / (inj_Gamma0 * (1 + inj_sigma)));

    Real b1 = (dM0 + dm + dm_inj) * con::c2;
    Real b2 = (dGamma_eff + Gamma_eff * (ad_idx - 1) / Gamma) * u;

    return (a1 + a2 + a3) / (b1 + b2);
}

/********************************************************************************************************************
 * METHOD: ForwardShockEqn::dUdt
 * DESCRIPTION: Computes the derivative of u with respect to time t.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
Real ForwardShockEqn<Jet, Injector>::dUdt(Real Gamma, Real u, Real r, Real theta, Real dGdt, Real drdt, Real dthetadt,
                                          Real ad_idx, Real rho) {
    Real dmdt = r * r * rho * drdt;
    Real dlnVdt = (3 / r * drdt - dGdt / Gamma);
    if (jet.spreading) {
        Real dm = medium.mass(r) / (4 * con::pi);
        Real factor = std::sin(theta) / (1 - std::cos(theta)) * dthetadt;
        dmdt = dmdt + dm * factor;
        dlnVdt += factor;
        dlnVdt += factor / (ad_idx - 1);
    }

    return (1 - eps_e) * (Gamma - 1) * con::c2 * dmdt - (ad_idx - 1) * dlnVdt * u;
}

// Updates the forward shock state at grid index (i, j, k) using the current state vector and the medium properties.
template <typename Eqn>
void updateForwardShock(size_t i, size_t j, int k, Eqn& eqn, const typename Eqn::StateArray& y, Shock& shock) {
    constFState state(y);

    Real n1 = eqn.medium.rho(state.r) / con::mp;                          // Compute upstream number density
    Real dN1dOmega = eqn.medium.mass(state.r) / (4 * con::pi * con::mp);  // number of proton per unit solid angle

    updateShockState(shock, i, j, k, state.r, state.theta, state.Gamma, state.Gamma, state.t_com, dN1dOmega, n1,
                     eqn.jet_sigma);
}

// Initializes the forward shock state vector at radius r0.
template <typename Eqn>
void setForwardInit(Eqn& eqn, typename Eqn::StateArray& y, Real t0) {
    FState state(y);
    state.Gamma = eqn.jet.Gamma0(eqn.phi, eqn.theta0, t0);  // Initial Lorentz factor
    Real beta0 = gammaTobeta(state.Gamma);
    state.r = beta0 * con::c * t0 / (1 - beta0);
    state.u = (state.Gamma - 1) * eqn.medium.mass(state.r) / (4 * con::pi) * con::c2;
    state.t_com = state.r / std::sqrt(state.Gamma * state.Gamma - 1) / con::c;
    state.theta = eqn.theta0;
    eqn.gamma4 = state.Gamma;
}

// Solves the forward shock evolution for a given shell (across radius values in array r) and updates the Shock object.
template <typename Eqn>
void solveForwardShell(size_t i, size_t j, const Array& t, Shock& shock, Eqn& eqn, double rtol) {
    using namespace boost::numeric::odeint;

    Real t0 = t[0];
    Real dt = (t[1] - t[0]) / 100;

    typename Eqn::StateArray y;
    FState state(y);
    setForwardInit(eqn, y, t0);  // Initialize state at starting radius

    if (y[0] < con::Gamma_cut) {  // If initial Lorentz factor is too low, exit early
        setStoppingShock(i, j, shock, t, state.r, state.theta);
        return;
    }

    // Create a dense output stepper for integrating the shock equations
    // auto stepper = bulirsch_stoer_dense_out<typename ShockEqn::StateArray>{0, rtol};
    auto stepper = make_dense_output(0, rtol, runge_kutta_dopri5<typename Eqn::StateArray>());

    stepper.initialize(y, t0, dt);

    Real t_back = t[t.size() - 1];  // Last radius in the array

    // Iterate over the radius array, updating the state and the Shock object as needed
    for (int k = 0; stepper.current_time() <= t_back;) {
        stepper.do_step(eqn);
        while (k < t.size() && stepper.current_time() > t[k]) {
            stepper.calc_state(t[k], y);
            updateForwardShock(i, j, k, eqn, y, shock);
            ++k;
        }
    }
}
#endif