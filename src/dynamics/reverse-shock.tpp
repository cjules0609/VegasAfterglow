//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "reverse-shock.hpp"

/********************************************************************************************************************
 * CONSTRUCTOR: FRShockEqn::FRShockEqn
 * DESCRIPTION: constructor for the FRShockEqn class.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
FRShockEqn<Ejecta, Medium>::FRShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, Real eps_e)
    : medium(medium),
      ejecta(ejecta),
      phi(phi),
      theta0(theta),
      eps_e(eps_e),
      gamma4(ejecta.Gamma0(phi, theta)),
      dE0dt(ejecta.dE0dOmega(phi, theta) / ejecta.T0),
      dM0dt(dE0dt / (gamma4 * con::c2)),
      u4(std::sqrt(gamma4 * gamma4 - 1) * con::c) {
    if constexpr (HasSigma<Ejecta>) {
        dM0dt /= 1 + ejecta.sigma0(phi, theta);
    }
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::operator()(State const& y, State& dydt, Real t)
 * DESCRIPTION: Reverse shock ODE.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
void FRShockEqn<Ejecta, Medium>::operator()(State const& state, State& diff, Real t) {
    Real gamma3 = 1, gamma_rel = 1;

    auto [dE_dt, dM_dt] = getInjection(t);
    diff.E_ej = dE_dt;
    diff.M_ej = dM_dt;
    bool is_injecting = diff.E_ej > 0 || diff.M_ej > 0;

    if (!crossed) {
        gamma3 = crossingGamma3(state);
        gamma_rel = relativeLorentz(gamma4, gamma3);

        Real sigma4 = shellMagnetization(state);
        diff.M3 = dM3dt(state.width, state.M_ej, gamma3, gamma4, sigma4);
        if (state.M3 >= state.M_ej) {
            diff.M3 = std::min(diff.M3, diff.M_ej);
        }
    } else {
        gamma_rel = crossedGamma_rel(state);
        gamma3 = crossedGamma3(gamma_rel, state.r);
        diff.M3 = 0;
    }

    Real beta3 = gammaTobeta(gamma3);
    diff.r = drdt(beta3);
    diff.t_com = dtdt_CoMoving(gamma3);

    if constexpr (!State::mass_profile) {
        diff.M_sw = medium.rho(phi, state.theta, state.r) * state.r * state.r * diff.r;
    }

    diff.width = is_injecting ? u4 : ShellWidthSpredingRate(gamma3, diff.t_com);
    diff.theta = 0;  // no lateral spreading
}

/********************************************************************************************************************
 * FUNCTION: setReverseInit
 * DESCRIPTION: Set the initial conditions for the reverse shock ODE.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
bool FRShockEqn<Ejecta, Medium>::setInitState(State& state, Real t0) const noexcept {
    Real beta4 = gammaTobeta(gamma4);
    state.theta = theta0;
    state.width = initShellWidth(gamma4, t0, ejecta.T0);

    Real dt = std::min(t0, ejecta.T0);
    state.E_ej = dE0dt * dt;
    state.M_ej = dM0dt * dt;

    state.r = beta4 * con::c * t0 / (1 - beta4);
    state.t_com = state.r / std::sqrt(gamma4 * gamma4 - 1) / con::c;

    Real rho = medium.rho(phi, theta0, state.r);

    if (t0 < ejecta.T0) {
        // thick shell deceleration radius where gamma starts to drop propto t^{-1/4}
        Real r_dec = 0.35 * std::sqrt(3 * dM0dt * (1 - beta4) / (beta4 * con::c * rho * gamma4));

        // r0 is larger than the decelertation radius, so r=beta4 * con::c * t0 / (1 - beta4) is not appropriate
        if (state.r > r_dec) {
            Real t_dec = r_dec * (1 - beta4) / (beta4 * con::c);
            state.r = std::sqrt(4 * gamma4 * gamma4 * r_dec * (t0 - t_dec) + r_dec * r_dec);
            rho = medium.rho(phi, theta0, state.r);

            Real t_dec_com = r_dec / std::sqrt(gamma4 * gamma4 - 1) / con::c;
            state.t_com =
                (std::pow(state.r, 1.5) - std::pow(r_dec, 1.5)) / (1.5 * gamma4 * std::sqrt(r_dec)) + t_dec_com;
        }
    }

    if constexpr (!State::mass_profile) {
        state.M_sw = state.r * state.r * state.r * rho / 3;
    } else {
        state.M_sw = medium.mass(phi, theta0, state.r);
    }

    auto [crossed, M3] = calcInitN3(*this, state, gamma4, t0);
    state.M3 = M3;
    return crossed;
}

template <typename Ejecta, typename Medium>
std::pair<Real, Real> FRShockEqn<Ejecta, Medium>::getInjection(Real t) const {
    Real dE_dt = 0;
    Real dM_dt = 0;

    if (t < ejecta.T0) {
        dE_dt = dE0dt;
        dM_dt = dM0dt;
    }

    if constexpr (State::energy_inject) {
        dE_dt += ejecta.dEdtdOmega(phi, theta0, t);
    }

    if constexpr (State::mass_inject) {
        dM_dt += ejecta.dMdtdOmega(phi, theta0, t);
    }

    return {dE_dt, dM_dt};
}

template <typename Ejecta, typename Medium>
bool FRShockEqn<Ejecta, Medium>::isInjecting(Real t) const {
    auto [dE_dt, dM_dt] = getInjection(t);
    return dE_dt > 0 || dM_dt > 0;
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::crossingGamma3(State const& state, Real t)
 * DESCRIPTION: Shock crossing gamma3 calculation.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::crossingGamma3(State const& state) const {
    Real M2 = sweptUpMass(*this, state);

    Real sigma4 = shellMagnetization(state);
    return calc_gamma3(gamma4, M2, state.M3, sigma4);
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::setCrossState(State const& state, Real B, Real t)
 * DESCRIPTION: Set some constants at the shock crossed point for later scaling calculation and switch the ODE from
 *              shock crossing state to crossed state.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
void FRShockEqn<Ejecta, Medium>::setCrossState(State const& state, Real B, Real t) {
    this->r_x = state.r;
    constexpr Real n3_norm = 1;  // normalized n3/n3_x = 1
    this->N0 = n3_norm * state.width * state.r * state.r;

    Real gamma3 = crossingGamma3(state);
    this->u_x = std::sqrt(gamma3 * gamma3 - 1);

    Real gamma_rel = relativeLorentz(gamma4, gamma3);
    Real gamma_electron = (gamma_rel - 1) * con::mp / con::me * eps_e + 1;  // electron Lorentz factor

    // use electron adiabatic index see section 3.2 https://arxiv.org/pdf/astro-ph/9910241
    this->ad_idx0 = adiabaticIndex(gamma_electron);
    Real p_norm = (ad_idx0 - 1) * (gamma_rel - 1) * n3_norm;      // p = (ad-1)(gamma-1)n m_p c^2
    this->adiabatic_const = std::pow(n3_norm, ad_idx0) / p_norm;  // p \propto n^ad
    this->Emag_const = p_norm / (B * B);
    this->crossed = true;
}

inline Real gParameter(Real gamma_rel, Real k = 0);

template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::crossedGamma3(Real gamma_rel, Real r) const {
    Real g = gParameter(gamma_rel);
    Real u = u_x * std::pow(r / r_x, -g);
    return std::sqrt(u * u + 1);
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::crossedGamma_rel(State const& state)
 * DESCRIPTION: Shock crossing gamma_43 calculation.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::crossedGamma_rel(State const& state) const {
    Real n3_norm = N0 / (state.width * state.r * state.r);        // proton number conservation
    Real p3_norm = std::pow(n3_norm, ad_idx0) / adiabatic_const;  // adiabatic expansion
    return p3_norm / ((ad_idx0 - 1) * n3_norm) + 1;
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::crossedB(State const& state)
 * DESCRIPTION: Post shock crossing magnetic field calculation.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::crossedB(State const& state) const {
    Real n3_norm = N0 / (state.width * state.r * state.r);        // proton number conservation
    Real p3_norm = std::pow(n3_norm, ad_idx0) / adiabatic_const;  // adiabatic expansion
    return std::sqrt(p3_norm / Emag_const);
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::shellMagnetization(State const& state)
 * DESCRIPTION: shell magnetization calculation.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
Real FRShockEqn<Ejecta, Medium>::shellMagnetization(State const& state) const {
    return std::max(0.0, state.E_ej / (gamma4 * state.M_ej * con::c2) - 1);
}
//---------------------------------------------------------------------------------------------------------------------
// Helper functions
//---------------------------------------------------------------------------------------------------------------------
/********************************************************************************************************************
 * METHOD: FRShockEqn::crossedGamma3(Real gamma_rel, Real r)
 * DESCRIPTION: Post crossing gamma3*beta3 profile that is \propto r^-g. Using gamma3*beta3 for Newtonian regime as
 *              well.
 ********************************************************************************************************************/
inline Real gParameter(Real gamma_rel, Real k) {
    constexpr Real g_low = 1.5;   // k is the medium power law index
    constexpr Real g_high = 3.5;  // Blandford-McKee limit// TODO: need to be modified for non ISM medium
    Real p = std::sqrt(std::sqrt(gamma_rel - 1));
    return g_low + (g_high - g_low) * p / (1 + p);
}

/********************************************************************************************************************
 * FUNCTION: calc_gamma3
 * DESCRIPTION: Calculate the reverse shock crossing gamma3
 ********************************************************************************************************************/
inline Real calc_gamma3(Real gamma4, Real M2, Real M3, Real sigma) {
    auto func = [=](Real gamma3) -> Real {
        Real gamma34 = relativeLorentz(gamma4, gamma3);
        Real adx3 = adiabaticIndex(gamma34);
        Real adx2 = adiabaticIndex(gamma3);
        Real g_eff2 = Gamma_eff(adx2, gamma3);
        Real g_eff3 = Gamma_eff(adx3, gamma3);

        Real E2 = M2 * (gamma3 - 1 + g_eff2 * (gamma3 - 1));
        Real E3 = M3 * (gamma3 - gamma4 + g_eff3 * (gamma34 - 1)) * (1 + sigma);

        return E2 + E3;
    };
    constexpr Real r_tol = 1e-3;
    return rootBisection(func, 1, gamma4, r_tol);
}

/********************************************************************************************************************
 * FUNCTION: calcInitN3
 * DESCRIPTION: calculate the initial N3 at t0 (to ensure the energy conservation)
 ********************************************************************************************************************/
template <typename Eqn, typename State>
inline auto calcInitN3(Eqn const& eqn, State const& state, Real gamma4, Real t0) {
    Real E_iso = eqn.ejecta.dE0dOmega(eqn.phi, eqn.theta0) * 4 * con::pi;
    Real n1 = eqn.medium.rho(eqn.phi, eqn.theta0, state.r) / con::mp;
    Real l = SedovLength(E_iso, n1);
    Real sigma0 = eqn.shellMagnetization(state);
    Real beta4 = gammaTobeta(gamma4);
    Real Delta0 = con::c * eqn.ejecta.T0;  // lab frame shell width
    Real r_x = std::sqrt(std::sqrt(Delta0 * l * l * l / 3) / (1 + sigma0) * 3. / 4);
    Real Rs = con::c * eqn.ejecta.T0 * beta4 / (1 - beta4);

    if (r_x > Rs) {
        constexpr Real rtol = 1e-3;
        Real r_x_new = (1 + 2 * rtol) * r_x;
        for (; std::abs((r_x - r_x_new) / r_x) > rtol;) {
            Delta0 = con::c * eqn.ejecta.T0 + (r_x_new - Rs) / (std::sqrt(3) * gamma4 * gamma4);
            r_x = r_x_new;
            r_x_new = std::sqrt(std::sqrt(Delta0 * l * l * l / 3) / (1 + sigma0) * 3. / 4);
        }
    }

    Real M4_tot = E_iso / (4 * con::pi * gamma4 * (1 + sigma0) * con::c2);

    if (!eqn.isInjecting(t0)) {
        Real M3 = M4_tot * std::pow(state.r / r_x, 1.5);
        if (M3 >= state.M_ej) {
            return std::make_tuple(true, state.M_ej);
        } else {
            return std::make_tuple(false, M3);
        }
    } else {
        Real M3 = M4_tot * state.r / r_x;
        return std::make_tuple(false, std::min(M3, state.M_ej));
    }
}

// comoving shell width at r (including the spreading from r=0 to r=r0)
inline Real initShellWidth(Real gamma, Real t0, Real T) {
    Real beta = gammaTobeta(gamma);
    if (t0 < T) {  // pure injection
        return gamma * t0 * beta * con::c;
    } else {  // injection+shell spreading
        Real cs = con::c / std::sqrt(3);
        return gamma * T * beta * con::c + cs * (t0 - T) * gamma;
    }
}
// inline Real initShellWidth(Real gamma, Real T, Real r) { return r / gamma / std::sqrt(3); }

/********************************************************************************************************************
 * INLINE FUNCTION: updateCrossedReverseShock
 * DESCRIPTION: Updates the post shock crossing state based on current ODE solution at t.
 ********************************************************************************************************************/
template <typename Eqn, typename State>
inline void updateCrossedReverseShock(size_t i, size_t j, size_t k, Eqn const& eqn, State const& state, Shock& shock) {
    size_t k0 = shock.injection_idx(i, j);
    Real r0 = shock.r(i, j, k0);
    shock.t_com(i, j, k) = state.t_com;
    shock.r(i, j, k) = state.r;
    shock.theta(i, j, k) = state.theta;
    shock.Gamma_rel(i, j, k) = eqn.crossedGamma_rel(state);
    shock.Gamma(i, j, k) = eqn.crossedGamma3(shock.Gamma_rel(i, j, k), state.r);
    shock.column_num_den(i, j, k) = shock.column_num_den(i, j, k0) * (r0 * r0) / (state.r * state.r);
    shock.B(i, j, k) = eqn.crossedB(state);
}

/********************************************************************************************************************
 * FUNCTION: set_f_state
 * DESCRIPTION: Set the initial condition for forward shock ODE solver based on the reverse shock state at crossed
 *              point.
 ********************************************************************************************************************/
template <typename Eqn, typename FState, typename RState>
void set_f_state(Eqn const& eqn_rvs, FState& state_fwd, RState const& state_rvs, Real gamma2) {
    state_fwd.r = state_rvs.r;
    state_fwd.t_com = state_rvs.t_com;
    state_fwd.theta = state_rvs.theta;
    state_fwd.Gamma = gamma2;

    Real M_sw = sweptUpMass(eqn_rvs, state_rvs);
    state_fwd.u = (gamma2 - 1) * M_sw * con::c2;

    if constexpr (!FState::mass_profile) {
        state_fwd.M_sw = M_sw;
    }

    if constexpr (FState::energy_inject) {
        state_fwd.E_ej = state_rvs.E_ej;
    }

    if constexpr (FState::mass_inject) {
        state_fwd.M_ej = state_rvs.M_ej;
    }
}

template <typename Eqn, typename State>
bool updateShockPair(size_t i, size_t j, int k, Eqn const& eqn_rvs, State const& state, Real t, Shock& shock_fwd,
                     Shock& shock_rvs) {
    Real n4 = state.M_ej / (state.r * state.r * state.width * con::mp);
    Real sigma4 = eqn_rvs.shellMagnetization(state);

    Real M2 = sweptUpMass(eqn_rvs, state);
    Real n1 = eqn_rvs.medium.rho(eqn_rvs.phi, state.theta, state.r) / con::mp;

    constexpr Real gamma1 = 1;
    constexpr Real sigma1 = 0;

    Real gamma3 = calc_gamma3(eqn_rvs.gamma4, M2, state.M3, sigma4);

    updateShockState(shock_fwd, i, j, k, state, gamma3, gamma1, M2 / con::mp, n1, sigma1);
    updateShockState(shock_rvs, i, j, k, state, gamma3, eqn_rvs.gamma4, state.M3 / con::mp, n4, sigma4);
    return state.M3 >= state.M_ej && !eqn_rvs.isInjecting(t);
}

// Handle low Gamma scenario
bool handleLowGamma(size_t i, size_t j, Shock& shock_fwd, Shock& shock_rvs, auto const& state_fwd,
                    auto const& state_rvs, auto const& eqn_rvs) {
    if (state_fwd.Gamma < con::Gamma_cut) {
        setStoppingShock(i, j, shock_fwd, state_fwd);
        setStoppingShock(i, j, shock_rvs, state_rvs);
        return true;
    }
    return false;
}

// Solve reverse shock until it crosses
template <typename View>
size_t solveShockPairUntilCross(size_t i, size_t j, View const& t, auto& stepper_rvs, auto& eqn_rvs, auto& state_rvs,
                                Shock& shock_fwd, Shock& shock_rvs) {
    size_t k0 = 0;
    Real t_back = t.back();
    bool crossed = false;
    size_t k = 0;

    while (!crossed && stepper_rvs.current_time() <= t_back) {
        stepper_rvs.do_step(eqn_rvs);
        while (k < t.size() && stepper_rvs.current_time() > t(k)) {
            stepper_rvs.calc_state(t(k), state_rvs);
            crossed = updateShockPair(i, j, k, eqn_rvs, state_rvs, t(k), shock_fwd, shock_rvs);
            if (crossed) {
                k0 = k;
                shock_rvs.injection_idx(i, j) = k0;
                eqn_rvs.setCrossState(state_rvs, shock_rvs.B(i, j, k), t(k));
                stepper_rvs.initialize(state_rvs, t(k0), stepper_rvs.current_time_step());
                break;
            }
            k++;
        }
    }

    return k0;
}

// Solve shock after crossing
template <typename UpdateFunc>
void solveShockAfterCross(size_t i, size_t j, auto const& t, size_t k0, auto& stepper, const auto& eqn, auto& state,
                          Shock& shock, UpdateFunc update) {
    size_t k = k0 + 1;
    Real t_back = t.back();

    while (stepper.current_time() <= t_back) {
        stepper.do_step(eqn);
        while (k < t.size() && stepper.current_time() > t(k)) {
            stepper.calc_state(t(k), state);
            update(i, j, k, eqn, state, shock);
            k++;
        }
    }
}

/********************************************************************************************************************
 * FUNCTION: solveFRShell
 * DESCRIPTION: Solve the reverse/forward shock ODE at grid (phi[i], theta[j]) as a function of t (on-axis
 *observation time).
 ********************************************************************************************************************/
template <typename FwdEqn, typename RvsEqn, typename View>
void solveFRShell(size_t i, size_t j, View const& t, Shock& shock_fwd, Shock& shock_rvs, FwdEqn const& eqn_fwd,
                  RvsEqn& eqn_rvs, Real rtol = 1e-6) {
    using namespace boost::numeric::odeint;
    auto stepper_fwd = make_dense_output(rtol, rtol, runge_kutta_dopri5<typename FwdEqn::State>());

    typename FwdEqn::State state_fwd;
    typename RvsEqn::State state_rvs;

    Real t0 = std::min(t.front(), 1 * con::sec);

    eqn_fwd.setInitState(state_fwd, t0);

    bool crossed = eqn_rvs.setInitState(state_rvs, t0);

    if (handleLowGamma(i, j, shock_fwd, shock_rvs, state_fwd, state_rvs, eqn_rvs)) return;

    if (crossed) {
        solveForwardShell(i, j, t, shock_fwd, eqn_fwd, rtol);
        return;
    }

    auto stepper_rvs = make_dense_output(rtol, rtol, runge_kutta_dopri5<typename RvsEqn::State>());
    stepper_rvs.initialize(state_rvs, t0, 0.01 * t0);

    size_t k0 = solveShockPairUntilCross(i, j, t, stepper_rvs, eqn_rvs, state_rvs, shock_fwd, shock_rvs);

    set_f_state(eqn_rvs, state_fwd, state_rvs, shock_rvs.Gamma(i, j, k0));
    stepper_fwd.initialize(state_fwd, t(k0), stepper_rvs.current_time_step());

    solveShockAfterCross(i, j, t, k0, stepper_fwd, eqn_fwd, state_fwd, shock_fwd,
                         updateForwardShock<FwdEqn, typename FwdEqn::State>);

    solveShockAfterCross(i, j, t, k0, stepper_rvs, eqn_rvs, state_rvs, shock_rvs,
                         updateCrossedReverseShock<RvsEqn, typename RvsEqn::State>);
}

/********************************************************************************************************************
 * FUNCTION: genFRShocks
 * DESCRIPTION: Generates a pair of forward and reverse shocks using the provided coordinates, medium, jet,
 *              and energy fractions. It creates two Shock objects (one forward, one reverse) and solves
 *              the shock shells for each phi, theta slice.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
ShockPair genFRShocks(Coord const& coord, Medium const& medium, Ejecta const& jet, Real eps_e_f, Real eps_B_f,
                      Real eps_e_r, Real eps_B_r, Real rtol = 1e-6) {
    auto [phi_size, theta_size, t_size] = coord.shape();
    size_t phi_size_needed = coord.t.shape()[0];
    Shock f_shock(phi_size_needed, theta_size, t_size, eps_e_f, eps_B_f);
    Shock r_shock(phi_size_needed, theta_size, t_size, eps_e_r, eps_B_r);
    for (size_t i = 0; i < phi_size_needed; ++i) {
        Real theta_s =
            jetSpreadingEdge(jet, medium, coord.phi(i), coord.theta.front(), coord.theta.back(), coord.t.front());
        for (size_t j = 0; j < theta_size; ++j) {
            // auto eqn_f = ForwardShockEqn(medium, jet, inject, coord.phi(i), coord.theta(j), eps_e);
            auto eqn_f = SimpleShockEqn(medium, jet, coord.phi(i), coord.theta(j), eps_e_f, theta_s);
            auto eqn_r = FRShockEqn(medium, jet, coord.phi(i), coord.theta(j), eps_e_r);
            // Solve the forward-reverse shock shell
            solveFRShell(i, j, xt::view(coord.t, i, j, xt::all()), f_shock, r_shock, eqn_f, eqn_r, rtol);
        }
    }

    return std::make_pair(std::move(f_shock), std::move(r_shock));
}
