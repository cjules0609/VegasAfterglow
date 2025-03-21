//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "forward-shock.hpp"

/********************************************************************************************************************
 * METHOD: ForwardShockEqn::operator()(State const& y, State& dydr, Real t)
 * DESCRIPTION: Computes the derivatives of the state variables with respect to radius t.
 ********************************************************************************************************************/
void ForwardShockEqn::operator()(StateArray const& y, StateArray& dydt, Real t) {
    FState state(y);
    FState diff(dydt);

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
ForwardShockEqn::ForwardShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, Real eps_e,
                                 Real theta_s)
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
Real ForwardShockEqn::dGammadt(Real t, constState const& state, State const& diff, Real ad_idx) {
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
Real ForwardShockEqn::dUdt(constState const& state, State const& diff, Real ad_idx) {
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
