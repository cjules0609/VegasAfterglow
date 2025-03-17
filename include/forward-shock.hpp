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
 * CLASS: FState
 * DESCRIPTION: This is a helper class to provide named access to the forward ODE state array (required by ODE solver)
 *              components.
 ********************************************************************************************************************/
struct FState {
    FState() = delete;

    template <typename StateArray>
    FState(StateArray& y)
        : Gamma(y[0]), u(y[1]), r(y[2]), t_com(y[3]), theta(y[4]), M_sw(y[5]), M_ej(y[6]), E_ej(y[7]) {}
    Real& Gamma;
    Real& u;
    Real& r;
    Real& t_com;
    Real& theta;
    Real& M_sw;
    Real& M_ej;
    Real& E_ej;
};

/********************************************************************************************************************
 * CLASS: constFState
 * DESCRIPTION: This is a helper class to provide named access to the forward ODE state array (required by ODE solver)
 *              components.
 ********************************************************************************************************************/
struct constFState {
    constFState() = delete;

    template <typename StateArray>
    constFState(StateArray const& y)
        : Gamma(y[0]), u(y[1]), r(y[2]), t_com(y[3]), theta(y[4]), M_sw(y[5]), M_ej(y[6]), E_ej(y[7]) {}
    Real const& Gamma;
    Real const& u;
    Real const& r;
    Real const& t_com;
    Real const& theta;
    Real const& M_sw;
    Real const& M_ej;
    Real const& E_ej;
};

/********************************************************************************************************************
 * CLASS: ForwardShockEqn
 * DESCRIPTION: Represents the forward shock equation for a given jet and medium. It defines a state vector
 *              (an array of 8 Reals) and overloads operator() to compute the derivatives of the state with
 *              respect to radius r. It also declares helper functions for the derivatives. Comprehensive models
 ********************************************************************************************************************/
template <typename Ejecta>
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

    ForwardShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, Real eps_e);

    Medium const& medium;  // Reference to the medium properties
    Ejecta const& ejecta;  // Reference to the ejecta properties
    Real const phi{0};     // Angular coordinate phi
    Real const theta0{0};  // Angular coordinate theta
    Real const eps_e{0};   // Electron energy fraction

    // forward shock ODE equation
    void operator()(StateArray const& y, StateArray& dydt, Real t);

   private:
    // dGammadt with respect to on-axis observe t.
    inline Real dGammadt(Real t, constFState const& state, FState const& diff, Real ad_idx);
    // dUdt with respect to on-axis observe t.
    inline Real dUdt(constFState const& state, FState const& diff, Real ad_idx);

    Real const dOmega0{0};  // Initial solid angle
};

/********************************************************************************************************************
 * METHOD: ForwardShockEqn::operator()(State const& y, State& dydr, Real t)
 * DESCRIPTION: Computes the derivatives of the state variables with respect to radius t.
 ********************************************************************************************************************/
template <typename Ejecta>
void ForwardShockEqn<Ejecta>::operator()(StateArray const& y, StateArray& dydt, Real t) {
    constFState state(y);
    FState diff(dydt);

    Real ad_idx = adiabaticIndex(state.Gamma);
    Real rho = medium.rho(phi, state.theta, state.r);
    Real beta = gammaTobeta(state.Gamma);
    Real uv = state.Gamma * beta;

    diff.r = drdt(beta);

    if (ejecta.spreading && state.theta < 0.5 * con::pi && uv * state.theta < 0.5) {
        diff.theta = dtheta_dt(uv, diff.r, state.r, state.Gamma);
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
template <typename Ejecta>
ForwardShockEqn<Ejecta>::ForwardShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, Real eps_e)
    : medium(medium), ejecta(ejecta), phi(phi), theta0(theta), eps_e(eps_e), dOmega0(1 - std::cos(theta)) {}

/********************************************************************************************************************
 * METHOD: ForwardShockEqn::dGammadt
 * DESCRIPTION: dGammadt with respect to on-axis observe t.
 ********************************************************************************************************************/
template <typename Ejecta>
Real ForwardShockEqn<Ejecta>::dGammadt(Real t, constFState const& state, FState const& diff, Real ad_idx) {
    Real Gamma2 = state.Gamma * state.Gamma;
    Real Gamma_eff = (ad_idx * (Gamma2 - 1) + 1) / state.Gamma;
    Real dGamma_eff = (ad_idx * (Gamma2 + 1) - 1) / Gamma2;
    Real dmdt = diff.M_sw;
    Real dlnVdt = 3 / state.r * diff.r;  // only r term

    Real M_sw = state.M_sw;  // Mass per unit solid angle from medium
    Real u = state.u;        // Internal energy per unit solid angle

    if (ejecta.spreading) {
        Real f_spread = (1 - std::cos(state.theta)) / dOmega0;
        dmdt = dmdt * f_spread + M_sw / dOmega0 * std::sin(state.theta) * diff.theta;
        M_sw *= f_spread;
        dlnVdt += std::sin(state.theta) / (1 - std::cos(state.theta)) * diff.theta;
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
template <typename Ejecta>
Real ForwardShockEqn<Ejecta>::dUdt(constFState const& state, FState const& diff, Real ad_idx) {
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
 ********************************************************************************************************************/
template <typename Eqn>
void updateForwardShock(size_t i, size_t j, int k, Eqn& eqn, const typename Eqn::StateArray& y, Shock& shock) {
    constFState state(y);

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
template <typename Eqn>
void setForwardInit(Eqn& eqn, typename Eqn::StateArray& y, Real t0) {
    FState state(y);
    state.Gamma = eqn.ejecta.Gamma(eqn.phi, eqn.theta0, 0);  // Initial Lorentz factor
    Real beta0 = gammaTobeta(state.Gamma);
    state.r = beta0 * con::c * t0 / (1 - beta0);
    state.M_sw = 0;                                          // swept mass per solid angle
    state.E_ej = eqn.ejecta.dE0dOmega(eqn.phi, eqn.theta0);  //
    state.M_ej = state.E_ej / (state.Gamma * (1 + eqn.ejecta.sigma0(eqn.phi, eqn.theta0)) * con::c2);
    state.u = (state.Gamma - 1) * state.M_sw * con::c2;
    state.t_com = state.r / std::sqrt(state.Gamma * state.Gamma - 1) / con::c;
    state.theta = eqn.theta0;
}

/********************************************************************************************************************
 * FUNCTION: solveForwardShell
 * DESCRIPTION: Solve the forward shock ODE at grid (phi[i], theta[j]) as a function of t (on-axis observation time).
 ********************************************************************************************************************/
// Solves the forward shock evolution for a given shell (across radius values in array r) and updates the Shock object.
template <typename Eqn>
void solveForwardShell(size_t i, size_t j, const Array& t, Shock& shock, Eqn& eqn, double rtol) {
    using namespace boost::numeric::odeint;

    Real t0 = t[0];
    Real dt = (t[1] - t[0]) / 100;

    typename Eqn::StateArray y;
    FState state(y);
    setForwardInit(eqn, y, t0);  // Initialize state at starting radius

    if (state.Gamma <= con::Gamma_cut) {                         // If initial Lorentz factor is too low, exit early
        setStoppingShock(i, j, shock, t, state.r, state.theta);  // Set the shock state to zero
        return;
    }

    // ODE solver with adaptive step size, and relative tolerance rtol
    // auto stepper = bulirsch_stoer_dense_out<typename ShockEqn::StateArray>{0, rtol};
    auto stepper = make_dense_output(0, rtol, runge_kutta_dopri5<typename Eqn::StateArray>());

    stepper.initialize(y, t0, dt);

    Real t_back = t[t.size() - 1];  // Last time in the array

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