//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "simple-shock.hpp"
/********************************************************************************************************************
 * CONSTRUCTOR: SimpleShockEqn::SimpleShockEqn
 * DESCRIPTION: Initializes a SimpleShockEqn object with references to the medium and ejecta
 *              along with the angular coordinates and energy fraction.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
SimpleShockEqn<Ejecta, Medium>::SimpleShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta,
                                               Real eps_e, Real theta_s)
    : medium(medium),
      ejecta(ejecta),
      phi(phi),
      theta0(theta),
      eps_e(eps_e),
      dOmega0(1 - std::cos(theta0)),
      theta_s(theta_s),
      M_ej(0) {
    M_ej = ejecta.dE0dOmega(phi, theta0) / ejecta.Gamma0(phi, theta0) / con::c2;
    if constexpr (HasSigma<Ejecta>) {
        M_ej /= 1 + ejecta.sigma0(phi, theta0);
    }
}

/********************************************************************************************************************
 * METHOD: SimpleShockEqn::operator()(State const& y, State& dydr, Real t)
 * DESCRIPTION: Computes the derivatives of the state variables with respect to t.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
void SimpleShockEqn<Ejecta, Medium>::operator()(State const& state, State& diff, Real t) const noexcept {
    Real beta = gammaTobeta(state.Gamma);

    diff.r = drdt(beta);
    diff.t_com = dtdt_CoMoving(state.Gamma);

    if (ejecta.spreading && state.theta < 0.5 * con::pi) {
        diff.theta = dtheta_dt(theta_s, state.theta, diff.r, state.r, state.Gamma);
    } else {
        diff.theta = 0;
    }

    if constexpr (State::mass_inject) {
        diff.M_ej = ejecta.dMdtdOmega(phi, theta0, t);
    }

    if constexpr (State::energy_inject) {
        diff.E_ej = ejecta.dEdtdOmega(phi, theta0, t);
    }

    Real rho = medium.rho(phi, state.theta, state.r);
    Real dMdt_sw = state.r * state.r * rho * diff.r;

    if constexpr (!State::mass_profile) {
        diff.M_sw = dMdt_sw;
    }

    diff.Gamma = dGammadt(dMdt_sw, state, diff);
}

/********************************************************************************************************************
 * METHOD: SimpleShockEqn::dGammadt
 * DESCRIPTION: Computes the derivative of Gamma with respect to radius t.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
Real SimpleShockEqn<Ejecta, Medium>::dGammadt(Real dMdt, State const& state, State const& diff) const noexcept {
    Real M_sw = 0;
    Real M_ej = this->M_ej;

    if constexpr (!State::mass_profile) {
        M_sw = state.M_sw;
    } else {
        M_sw = medium.mass(phi, state.theta, state.r);
    }

    if (ejecta.spreading) {
        Real f_spread = (1 - std::cos(state.theta)) / dOmega0;
        dMdt = dMdt * f_spread + M_sw / dOmega0 * std::sin(state.theta) * diff.theta;
        M_sw *= f_spread;
    }

    double a1 = (1 - state.Gamma * state.Gamma) * dMdt;

    if constexpr (State::energy_inject) {
        a1 += diff.E_ej / con::c2;
    }

    if constexpr (State::mass_inject) {
        a1 -= state.Gamma * diff.M_ej;
        M_ej = state.M_ej;
    }

    return a1 / (M_ej + eps_e * M_sw + 2 * (1 - eps_e) * state.Gamma * M_sw);
}

template <typename Ejecta, typename Medium>
void SimpleShockEqn<Ejecta, Medium>::setInitState(State& state, Real t0) const noexcept {
    state.Gamma = ejecta.Gamma0(phi, theta0);

    Real beta0 = gammaTobeta(state.Gamma);
    state.r = beta0 * con::c * t0 / (1 - beta0);

    state.t_com = state.r / std::sqrt(state.Gamma * state.Gamma - 1) / con::c;

    state.theta = theta0;

    if constexpr (State::energy_inject) {
        state.E_ej = ejecta.dE0dOmega(phi, theta0);
    }

    if constexpr (State::mass_inject) {
        state.M_ej = M_ej;
    }

    if constexpr (!State::mass_profile) {
        state.M_sw = medium.rho(phi, theta0, state.r) * state.r * state.r * state.r / 3;
    }
}