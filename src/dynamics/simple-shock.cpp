//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "simple-shock.h"

/********************************************************************************************************************
 * METHOD: SimpleShockEqn::operator()(State const& y, State& dydr, Real t)
 * DESCRIPTION: Computes the derivatives of the state variables with respect to t.
 *              The state vector components are:
 *                  y[0] - Gamma (Lorentz factor)
 *                  y[1] - u (internal energy per solid angle)
 *                  y[2] - r (radius)
 *                  y[3] - t_com (co-moving time) [unused here]
 *                  y[4] - theta_jet (jet opening angle)
 *                  y[5] - m (swept mass per solid angle)
 *                  y[6] - M (ejecta mass per solid angle)
 *                  y[7] - E_ejecta (ejecta energy per solid angle)
 ********************************************************************************************************************/
void SimpleShockEqn::operator()(StateArray const& y, StateArray& dydt, Real t) {
    FState state(y);
    FState diff(dydt);

    Real rho = medium.rho(phi, state.theta, state.r);  // Get medium density
    Real beta = gammaTobeta(state.Gamma);              // Convert Gamma to beta (velocity/c)

    diff.r = drdt(beta);
    if (ejecta.spreading && state.theta < 0.5 * con::pi) {
        diff.theta = dtheta_dt(theta_s, state.theta, diff.r, state.r, state.Gamma);
    } else {
        diff.theta = 0;
    }

    diff.u = 0;
    diff.t_com = dtdt_CoMoving(state.Gamma);
    diff.M_sw = state.r * state.r * rho * diff.r;
    // d(M)/dt, use theta0 instead of theta as the mass injected from the bottom
    diff.M_ej = ejecta.dMdtdOmega(phi, theta0, t);
    diff.E_ej = ejecta.dEdtdOmega(phi, theta0, t);
    diff.Gamma = dGammadt(t, state, diff);
}

/********************************************************************************************************************
 * CONSTRUCTOR: SimpleShockEqn::SimpleShockEqn
 * DESCRIPTION: Initializes a SimpleShockEqn object with references to the medium and ejecta
 *              along with the angular coordinates and energy fraction.
 ********************************************************************************************************************/
SimpleShockEqn::SimpleShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta_lo, Real theta,
                               Real eps_e, Real theta_s)
    : medium(medium),
      ejecta(ejecta),
      phi(phi),
      theta_lo(theta_lo),
      theta0(theta),
      eps_e(eps_e),
      dOmega0(std::cos(theta_lo) - std::cos(theta0)),
      theta_s(theta_s) {}

/********************************************************************************************************************
 * METHOD: SimpleShockEqn::dGammadt
 * DESCRIPTION: Computes the derivative of Gamma with respect to radius t.
 ********************************************************************************************************************/
Real SimpleShockEqn::dGammadt(Real t, FState<const StateArray> const& state, FState<StateArray> const& diff) {
    Real M_sw = state.M_sw;
    Real dMdt = diff.M_sw;

    if (ejecta.spreading) {
        Real f_spread = (std::cos(theta_lo) - std::cos(state.theta)) / dOmega0;
        dMdt = dMdt * f_spread + M_sw / dOmega0 * std::sin(state.theta) * diff.theta;
        M_sw *= f_spread;
    }

    double a1 = (1 - state.Gamma * state.Gamma) * dMdt;
    double a2 = diff.E_ej / con::c2 - state.Gamma * diff.M_ej;
    return (a1 + a2) / (state.M_ej + eps_e * M_sw + 2 * (1 - eps_e) * state.Gamma * M_sw);
}
