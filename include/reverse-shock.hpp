//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _REVERSESHOCK_
#define _REVERSESHOCK_

#include "shock.h"
#include "utilities.h"

/********************************************************************************************************************
 * CLASS: RState
 * DESCRIPTION: This is a helper class to provide named access to the reverse shock ODE state array (required by ODE
 *              solver) components.
 ********************************************************************************************************************/
struct RState {
    RState() = delete;

    template <typename StateArray>
    RState(StateArray& y)
        : width(y[0]), N3(y[1]), r(y[2]), t_com(y[3]), theta(y[4]), M_sw(y[5]), M_ej(y[6]), E_ej(y[7]) {}
    Real& width;
    Real& N3;
    Real& r;
    Real& t_com;
    Real& theta;
    Real& M_sw;
    Real& M_ej;
    Real& E_ej;
};

/********************************************************************************************************************
 * CLASS: constRState
 * DESCRIPTION: This is a helper class to provide named access to the reverse shock ODE state array (required by ODE
 *              solver) components.
 ********************************************************************************************************************/
struct constRState {
    constRState() = delete;

    template <typename StateArray>
    constRState(StateArray const& y)
        : width(y[0]), N3(y[1]), r(y[2]), t_com(y[3]), theta(y[4]), M_sw(y[5]), M_ej(y[6]), E_ej(y[7]) {}
    Real const& width;
    Real const& N3;
    Real const& r;
    Real const& t_com;
    Real const& theta;
    Real const& M_sw;
    Real const& M_ej;
    Real const& E_ej;
};

/********************************************************************************************************************
 * CLASS: FRShockEqn
 * DESCRIPTION: Represents the reverse shock (or forward-reverse shock) equation for a given Jet and medium.
 *              It defines a state vector (an array of 8 Reals) and overloads operator() to compute the
 *              derivatives of the state with respect to radius r. It also declares a helper function to compute
 *              the derivative of N3 (number per solid angle) with respect to r.
 ********************************************************************************************************************/
template <typename Ejecta>
class FRShockEqn {
   public:
    // State vector for reverse shock variables [D_jet, N3, r, t_com, theta]
    // - D_jet: comoving shell width,
    // - N3: electron number per unit solid angle in region 3
    // - r: radius
    // - t_com: comoving time
    // - theta: jet opening angle
    // - M_sw: swept mass per solid angle
    // - M_ej: ejecta mass per solid angle
    // - E_ej: ejecta energy per solid angle
    using StateArray = std::array<Real, 8>;

    FRShockEqn(Medium const& medium, Ejecta const& jet, Real phi, Real theta, Real eps_e);

    Medium const& medium;  // Reference to the medium properties
    Ejecta const& ejecta;  // Reference to the jet properties
    Real const phi{0};     // Angular coordinate phi
    Real const theta0{0};  // Angular coordinate theta
    Real const eps_e{0};   // Electron energy fraction
    Real u_x{0};
    Real r_x{0};

    // Reverse shock ODE equation
    void operator()(StateArray const& y, StateArray& dydt, Real t);

    // Set the shock state when the reverse shock crosses the jet.
    template <typename State>
    void setCrossState(State const& state, Real B, Real t);

    // calculate the Gamma3 during the shock crossing phase.
    template <typename State>
    Real crossingGamma3(State const& state, Real t);

    // calculate the Gamma_43 post shock crossing.
    template <typename State>
    Real crossedGamma_rel(State const& state);

    // calculate the magnetic field post shock crossing.
    template <typename State>
    Real crossedB(State const& state);

    // calculate the Gamma3 post shock crossing.
    Real crossedGamma3(Real gamma_rel, Real r);

   private:
    Real N0{0};               // normalized total electron (for post crossing scaling calculation).
    Real adiabatic_const{1};  // normalized adiabatic constant where C = rho^idx/p.
    Real Emag_const{1};       // normalized magnetic energy constant where C = B^2/p.
    Real ad_idx0{4. / 3};     // adiabatic index at the shock crossing.
    bool crossed{false};
};

/********************************************************************************************************************
 * FUNCTION: gParameter
 * DESCRIPTION: post shock crossing power law index for the Gamma3*beta3 \propto r^-g.
 ********************************************************************************************************************/
inline Real gParameter(Real gamma_rel, Real k = 0) {
    const Real g_low = 3 / 2.0;   // k is the medium power law index
    constexpr Real g_high = 3.5;  // Blandford-McKee limit// TODO: need to be modified for non ISM medium
    Real p = std::sqrt(std::sqrt(gamma_rel - 1));
    return g_low + (g_high - g_low) * p / (1 + p);
}

inline Real Gamma_eff(Real adx, Real Gamma) { return (adx * Gamma * Gamma - adx + 1) / Gamma; }

inline Real calc_gamma3(Real gamma4, Real N2, Real N3, Real N4, Real E40, Real sigma) {
    Real Einit = E40 + N2;
    Real E4 = N4 * gamma4 * (1 + sigma);
    // TODO
    auto func = [=](Real gamma3) -> Real {
        Real gamma34 = relativeLorentz(gamma4, gamma3);
        Real adx3 = adiabaticIndex(gamma34);
        Real adx2 = adiabaticIndex(gamma3);
        Real g_eff2 = Gamma_eff(adx2, gamma3);
        Real g_eff3 = Gamma_eff(adx3, gamma3);

        Real E2 = N2 * (gamma3 + g_eff2 * (gamma3 - 1));
        Real E3 = N3 * (gamma3 + g_eff3 * (gamma34 - 1)) * (1 + sigma);

        return E2 + E3 + E4 - Einit;
    };

    return rootBisection(func, 1, gamma4, 1e-4);
}

/********************************************************************************************************************
 * FUNCTION: calcInitN3
 * DESCRIPTION: calculate the initial N3 at t0 (to ensure the energy conservation)
 ********************************************************************************************************************/
template <typename Eqn>
inline auto calcInitN3(Eqn& eqn, RState const& state, Real gamma4, Real t0) {
    Real n1 = eqn.medium.rho(eqn.phi, eqn.theta0, state.r) / con::mp;
    Real N4 = state.M_ej / con::mp;
    Real n4 = N4 / (state.r * state.r * state.width);
    Real sigma_p1 = state.E_ej / (gamma4 * state.M_ej * con::c2);

    Real drdt_ = drdt(gammaTobeta(gamma4));
    Real dN3dt_ = dN3dt(state.r, n1, n4, gamma4, drdt_, gamma4, 0) * std::sqrt(sigma_p1);
    Real N3 = dN3dt_ * t0;
    return std::make_tuple(N3 >= N4, std::min(N3, N4));
}

inline Real initShellWidth(Real gamma, Real T, Real r) { return gamma * T * con::c + r / gamma / std::sqrt(3); }

/********************************************************************************************************************
 * FUNCTION: setReverseInit
 * DESCRIPTION: Set the initial conditions for the reverse shock ODE.
 ********************************************************************************************************************/
template <typename Eqn>
bool setReverseInit(Eqn& eqn, typename Eqn::StateArray& y, Real t0) {
    RState state(y);
    Real gamma40 = eqn.ejecta.Gamma(eqn.phi, eqn.theta0, 0);  // Obtain initial Gamma from the jet
    Real beta4 = gammaTobeta(gamma40);
    state.r = beta4 * con::c * t0 / (1 - beta4);
    state.t_com = state.r / std::sqrt(gamma40 * gamma40 - 1) / con::c;
    // initial width + spreading from r=0 to r=r0 (with sound speed of c/sqrt(3))
    state.width = initShellWidth(gamma40, eqn.ejecta.T0, state.r);
    auto [crossed, N3] = calcInitN3(eqn, state, gamma40, t0);
    state.N3 = N3;
    state.theta = eqn.theta0;
    state.M_sw = 0;
    state.E_ej = eqn.ejecta.dE0dOmega(eqn.phi, eqn.theta0);
    state.M_ej = state.E_ej / (gamma40 * (1 + eqn.ejecta.sigma0(eqn.phi, eqn.theta0)) * con::c2);

    return crossed;
}

/********************************************************************************************************************
 * CONSTRUCTOR: FRShockEqn::FRShockEqn
 * DESCRIPTION: constructor for the FRShockEqn class.
 ********************************************************************************************************************/
template <typename Ejecta>
FRShockEqn<Ejecta>::FRShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, Real eps_e)
    : medium(medium), ejecta(ejecta), phi(phi), theta0(theta), eps_e(eps_e) {}

/********************************************************************************************************************
 * METHOD: FRShockEqn::operator()(State const& y, State& dydt, Real t)
 * DESCRIPTION: Reverse shock ODE.
 ********************************************************************************************************************/
template <typename Ejecta>
void FRShockEqn<Ejecta>::operator()(StateArray const& y, StateArray& dydt, Real t) {
    constRState state(y);
    RState diff(dydt);

    Real gamma3 = 1;
    Real beta3 = 1;
    Real gamma_rel = 1;

    if (!crossed) {
        Real n1 = medium.rho(phi, state.theta, state.r) / con::mp;

        Real N40 = state.M_ej / con::mp;
        Real n4 = N40 / (state.r * state.r * state.width);
        Real N2 = state.M_sw / con::mp;
        Real gamma4 = ejecta.Gamma(phi, theta0, t);
        Real sigma4 = state.E_ej / (gamma4 * state.M_ej * con::c2) - 1;

        gamma3 = calc_gamma3(gamma4, N2, state.N3, N40 - state.N3, state.E_ej / (con::mp * con::c2), sigma4);
        beta3 = gammaTobeta(gamma3);

        gamma_rel = relativeLorentz(gamma4, gamma3);
        diff.r = drdt(beta3);
        diff.N3 = dN3dt(state.r, n1, n4, gamma3, diff.r, gamma4, sigma4);
    } else {
        gamma_rel = crossedGamma_rel(state);
        gamma3 = crossedGamma3(gamma_rel, state.r);
        beta3 = gammaTobeta(gamma3);
        diff.r = drdt(beta3);
        diff.N3 = 0;
    }

    diff.t_com = dtdt_CoMoving(gamma3);
    diff.width = dDdt_Jet(gamma_rel, diff.t_com);
    Real u3 = gamma3 * beta3;
    if (ejecta.spreading && state.theta < 0.5 * con::pi && u3 * state.theta < 0.5) {
        diff.theta = dtheta_dt(u3, diff.r, state.r, gamma3);
    } else {
        diff.theta = 0;
    }
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::setCrossState(State const& state, Real B, Real t)
 * DESCRIPTION: Set some constants at the shock crossed point for later scaling calculation and switch the ODE from
 *              shock crossing state to crossed state.
 ********************************************************************************************************************/
template <typename Ejecta>
template <typename State>
void FRShockEqn<Ejecta>::setCrossState(State const& state, Real B, Real t) {
    this->r_x = state.r;
    constexpr Real n3_norm = 1;  // normalized n3/n3_x = 1
    this->N0 = n3_norm * state.width * state.r * state.r;

    Real gamma3 = crossingGamma3(state, t);
    this->u_x = std::sqrt(gamma3 * gamma3 - 1);

    Real gamma4 = ejecta.Gamma(phi, theta0, t);
    Real gamma_rel = relativeLorentz(gamma4, gamma3);
    Real gamma_electron = (gamma_rel - 1) * con::mp / con::me * eps_e + 1;  // electron Lorentz factor

    // use electron adiabatic index see section 3.2 https://arxiv.org/pdf/astro-ph/9910241
    this->ad_idx0 = adiabaticIndex(gamma_electron);
    Real p_norm = (ad_idx0 - 1) * (gamma_rel - 1) * n3_norm;      // p = (ad-1)(gamma-1)n m_p c^2
    this->adiabatic_const = std::pow(n3_norm, ad_idx0) / p_norm;  // p \propto n^ad
    this->Emag_const = p_norm / (B * B);
    this->crossed = true;
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::crossedGamma3(Real gamma_rel, Real r)
 * DESCRIPTION: Post crossing gamma3*beta3 profile that is \propto r^-g. Using gamma3*beta3 for Newtonian regime as
 *              well.
 ********************************************************************************************************************/
template <typename Ejecta>
Real FRShockEqn<Ejecta>::crossedGamma3(Real gamma_rel, Real r) {
    Real g = gParameter(gamma_rel);
    Real u = u_x * std::pow(r / r_x, -g);
    return std::sqrt(u * u + 1);
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::crossingGamma3(State const& state, Real t)
 * DESCRIPTION: Shock crossing gamma3 calculation.
 ********************************************************************************************************************/
template <typename Ejecta>
template <typename State>
Real FRShockEqn<Ejecta>::crossingGamma3(State const& state, Real t) {
    Real N2 = state.M_sw / con::mp;
    Real N40 = state.M_ej / con::mp;
    Real gamma4 = ejecta.Gamma(phi, theta0, t);
    Real sigma4 = state.E_ej / (gamma4 * state.M_ej * con::c2) - 1;
    return calc_gamma3(gamma4, N2, state.N3, N40 - state.N3, state.E_ej / (con::mp * con::c2), sigma4);
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::crossedGamma_rel(State const& state)
 * DESCRIPTION: Shock crossing gamma_43 calculation.
 ********************************************************************************************************************/
template <typename Ejecta>
template <typename State>
Real FRShockEqn<Ejecta>::crossedGamma_rel(State const& state) {
    Real n3_norm = N0 / (state.width * state.r * state.r);        // proton number conservation
    Real p3_norm = std::pow(n3_norm, ad_idx0) / adiabatic_const;  // adiabatic expansion
    return p3_norm / ((ad_idx0 - 1) * n3_norm) + 1;
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::crossedB(State const& state)
 * DESCRIPTION: Post shock crossing magnetic field calculation.
 ********************************************************************************************************************/
template <typename Ejecta>
template <typename State>
Real FRShockEqn<Ejecta>::crossedB(State const& state) {
    Real n3_norm = N0 / (state.width * state.r * state.r);        // proton number conservation
    Real p3_norm = std::pow(n3_norm, ad_idx0) / adiabatic_const;  // adiabatic expansion

    return std::sqrt(p3_norm / Emag_const);
}

/********************************************************************************************************************
 * INLINE FUNCTION: updateCrossedReverseShock
 * DESCRIPTION: Updates the post shock crossing state based on current ODE solution at t.
 ********************************************************************************************************************/
template <typename Eqn>
inline void updateCrossedReverseShock(size_t i, size_t j, size_t k, size_t k0, Eqn& eqn, RState const& state,
                                      Shock& shock) {
    Real r0 = shock.r[i][j][k0];
    shock.t_com[i][j][k] = state.t_com;
    shock.r[i][j][k] = state.r;
    shock.theta[i][j][k] = state.theta;
    shock.Gamma_rel[i][j][k] = eqn.crossedGamma_rel(state);
    shock.Gamma[i][j][k] = eqn.crossedGamma3(shock.Gamma_rel[i][j][k], state.r);
    shock.column_num_den[i][j][k] = shock.column_num_den[i][j][k0] * (r0 * r0) / (state.r * state.r);
    shock.B[i][j][k] = eqn.crossedB(state);
}

/********************************************************************************************************************
 * FUNCTION: reverseShockExists
 * DESCRIPTION: Check if reverse shock can be generated.
 ********************************************************************************************************************/
template <typename Eqn>
bool reverseShockExists(Eqn const& eqn, Real t0) {           // TODO: for general t case
    Real gamma = eqn.ejecta.Gamma(eqn.phi, eqn.theta0, t0);  // Initial Lorentz factor
    Real beta0 = gammaTobeta(gamma);
    Real r0 = beta0 * con::c * t0 / (1 - beta0);
    Real D_jet0 = initShellWidth(gamma, eqn.ejecta.duration, r0);

    Real n4 = calc_n4(eqn.ejecta.dEdOmega(eqn.phi, eqn.theta0, t0), gamma, r0, D_jet0, eqn.jet_sigma);
    Real n1 = eqn.medium.rho(eqn.phi, eqn.theta0, r0) / con::mp;
    return eqn.jet_sigma < 8. / 3 * (gamma * gamma - 1) * n1 / n4;
}

/********************************************************************************************************************
 * FUNCTION: set_f_state
 * DESCRIPTION: Set the initial condition for forward shock ODE solver based on the reverse shock state at crossed
 *              point.
 ********************************************************************************************************************/
template <typename Eqn, typename FStateArray, typename RStateArray>
void set_f_state(Eqn& eqn, FStateArray& y_fwd, RStateArray const& y_rvs, Real gamma2) {
    FState state_fwd(y_fwd);
    constRState state_rvs(y_rvs);
    state_fwd.r = state_rvs.r;
    state_fwd.t_com = state_rvs.t_com;
    state_fwd.theta = state_rvs.theta;
    state_fwd.Gamma = gamma2;
    state_fwd.M_sw = state_rvs.M_sw;
    state_fwd.M_ej = state_rvs.M_ej;
    state_fwd.E_ej = state_rvs.E_ej;
    state_fwd.u = (gamma2 - 1) * state_fwd.M_sw * con::c2;
}

template <typename Eqn, typename State>
bool updateForwardReverseShock(size_t i, size_t j, int k, Eqn& eqn_rvs, State const& state, Real t, Shock& shock_fwd,
                               Shock& shock_rvs) {
    Real N4 = state.M_ej / con::mp;
    Real n4 = N4 / (state.r * state.r * state.width);
    Real gamma4 = eqn_rvs.ejecta.Gamma(eqn_rvs.phi, eqn_rvs.theta0, t);
    Real sigma4 = state.E_ej / (gamma4 * state.M_ej * con::c2) - 1;

    Real N2 = state.M_sw / con::mp;
    Real n1 = eqn_rvs.medium.rho(eqn_rvs.phi, state.theta, state.r) / con::mp;

    Real gamma3 = calc_gamma3(gamma4, N2, state.N3, N4 - state.N3, state.E_ej / (con::mp * con::c2), sigma4);

    updateShockState(shock_fwd, i, j, k, state, gamma3, 1., N2, n1, 0.);

    updateShockState(shock_rvs, i, j, k, state, gamma3, gamma4, state.N3, n4, sigma4);
    return state.N3 >= N4;
}

/********************************************************************************************************************
 * FUNCTION: solveFRShell
 * DESCRIPTION: Solve the reverse/forward shock ODE at grid (phi[i], theta[j]) as a function of t (on-axis observation
 *              time).
 ********************************************************************************************************************/
template <typename FShockEqn, typename RShockEqn>
void solveFRShell(size_t i, size_t j, Array const& t, Shock& shock_fwd, Shock& shock_rvs, FShockEqn& eqn_fwd,
                  RShockEqn& eqn_rvs, Real rtol = 1e-9) {
    using namespace boost::numeric::odeint;
    // auto stepper_fwd = bulirsch_stoer_dense_out<typename FShockEqn::StateArray>{0, rtol};
    auto stepper_fwd = make_dense_output(0, rtol, runge_kutta_dopri5<typename FShockEqn::StateArray>());

    bool crossed = false;
    Real t0 = t[0];
    Real dt = (t[1] - t[0]) / 100;
    typename RShockEqn::StateArray y_rvs;
    typename FShockEqn::StateArray y_fwd;

    FState state_fwd(y_fwd);
    RState state_rvs(y_rvs);

    setForwardInit(eqn_fwd, y_fwd, t0);

    crossed = setReverseInit(eqn_rvs, y_rvs, t0);

    if (state_fwd.Gamma < con::Gamma_cut) {
        setStoppingShock(i, j, shock_fwd, t, state_rvs.r, eqn_rvs.theta0);
        setStoppingShock(i, j, shock_rvs, t, state_rvs.r, eqn_rvs.theta0);
        return;
    }

    if (crossed) {
        solveForwardShell(i, j, t, shock_fwd, eqn_fwd, rtol);
    } else {
        // auto stepper_rvs = bulirsch_stoer_dense_out<typename RShockEqn::StateArray>{0, rtol};
        auto stepper_rvs = make_dense_output(0, rtol, runge_kutta_dopri5<typename RShockEqn::StateArray>());
        stepper_rvs.initialize(y_rvs, t0, dt);
        Real t_back = t[t.size() - 1];

        // shock crossing
        size_t k0 = 0;
        for (size_t k = 0; !crossed && stepper_rvs.current_time() <= t_back;) {
            stepper_rvs.do_step(eqn_rvs);
            for (; k < t.size() && stepper_rvs.current_time() > t[k]; k++) {
                stepper_rvs.calc_state(t[k], y_rvs);
                crossed = updateForwardReverseShock(i, j, k, eqn_rvs, state_rvs, t[k], shock_fwd, shock_rvs);
                if (crossed) {
                    k0 = k;  // k0 is the index of the first element in the array that the reverse shock crosses
                    eqn_rvs.setCrossState(state_rvs, shock_rvs.B[i][j][k], t[k]);
                    stepper_rvs.initialize(y_rvs, t[k0], stepper_rvs.current_time_step());
                    break;
                }
            }
        }

        // crossed forward shock evolution

        set_f_state(eqn_fwd, y_fwd, y_rvs, shock_rvs.Gamma[i][j][k0]);
        stepper_fwd.initialize(y_fwd, t[k0], stepper_rvs.current_time_step());
        for (size_t k = k0 + 1; stepper_fwd.current_time() <= t_back;) {
            stepper_fwd.do_step(eqn_fwd);

            for (; k < t.size() && stepper_fwd.current_time() > t[k]; k++) {
                stepper_fwd.calc_state(t[k], y_fwd);
                updateForwardShock(i, j, k, eqn_fwd, y_fwd, shock_fwd);
            }
        }

        // crossed reverse shock evolution
        for (size_t k = k0 + 1; stepper_rvs.current_time() <= t_back;) {
            stepper_rvs.do_step(eqn_rvs);
            for (; k < t.size() && stepper_rvs.current_time() > t[k]; k++) {
                stepper_rvs.calc_state(t[k], y_rvs);
                updateCrossedReverseShock(i, j, k, k0, eqn_rvs, y_rvs, shock_rvs);
            }
        }
    }
}

#endif