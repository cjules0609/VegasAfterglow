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
 * CLASS: FRShockEqn
 * DESCRIPTION: Represents the reverse shock (or forward-reverse shock) equation for a given Jet and Injector.
 *              It defines a state vector (an array of 5 Reals) and overloads operator() to compute the
 *              derivatives of the state with respect to radius r. It also declares a helper function to compute
 *              the derivative of N3 (number per solid angle) with respect to r.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
class FRShockEqn {
   public:
    using State =
        std::array<Real, 6>;  // State vector for reverse shock variables [ignored, N3, r, t_com, D_jet, theta]

    FRShockEqn(Medium const& medium, Jet const& jet, Injector const& inject, Real phi, Real theta);

    Medium const& medium;     // Reference to the medium properties
    Jet const& jet;           // Reference to the jet properties
    Injector const& inject;   // Reference to the injector properties
    Real const phi{0};        // Angular coordinate phi
    Real const theta0{0};     // Angular coordinate theta
    Real const jet_sigma{0};  // Jet magnetization parameter
    Real gamma4{1};           // Initial Gamma parameter from the jet

    // Overloaded operator() to compute the derivatives of the state vector with respect to time t.
    void operator()(State const& y, State& dydt, Real t);

   private:
    // Helper function: computes the derivative of N3 (number per solid angle) with respect to t.
    inline Real dN3dtPerOmega(Real r, Real n1, Real n4, Real gamma3, Real drdt);
};

inline Real sigmoid(Real x, Real x0, Real w) { return 1.0 / (1.0 + std::exp(-w * (x - x0))); }

inline Real mechanicalModelCorrection(Real gamma34, Real sigma, Real k) {
    Real S = sigmoid(gamma34, 10, 5);
    return (1 - S) * ((5 - k) / (3 - k)) + S * ((22 - 6 * k) / (12 - 3 * k));
}
/********************************************************************************************************************
 * INLINE FUNCTION: calc_gamma3
 * DESCRIPTION: Computes gamma3 for the reverse shock based on radius, upstream and downstream densities,
 *              the jet Lorentz factor (gamma4), and magnetization (sigma).
 ********************************************************************************************************************/
inline Real calc_gamma3(Real n1, Real n4, Real gamma4, Real sigma, Real k) {
    Real C = (1 + sigma) * n4 / n1;

    auto func = [=](Real g3) -> Real {
        Real gamma34 = relativeLorentz(gamma4, g3);
        Real G = mechanicalModelCorrection(gamma34, sigma, k);
        return G * C * (gamma34 * gamma34 - 1) - (g3 * g3 - 1);
    };

    return rootBisection(func, 1, gamma4, 1e-3);
    /* Real gamma3 =
         gamma4 * std::sqrt(((C - 2) - 2 * std::sqrt(C * (gamma4 * gamma4 - 1) + 1)) / (C - 4 * gamma4 * gamma4));
     return gamma3;*/
}

/********************************************************************************************************************
 * FUNCTION: setReverseInit
 * DESCRIPTION: Initializes the state vector for reverse shock evolution at the given radius t0.
 ********************************************************************************************************************/
template <typename ShockEqn>
void setReverseInit(ShockEqn& eqn, typename ShockEqn::State& state, Real t0) {
    Real gamma4 = eqn.jet.Gamma0(eqn.phi, eqn.theta0, t0);  // Obtain initial Gamma from the jet
    Real beta0 = gammaTobeta(gamma4);
    Real r0 = beta0 * con::c * t0 / (1 - beta0);
    Real t_com0 = r0 / std::sqrt(gamma4 * gamma4 - 1) / con::c;
    Real D_jet0 = con::c * eqn.jet.duration;
    Real dN3dOmega = 0;  // Initialize number per unit solid angle to zero
    state = {0., dN3dOmega, r0, t_com0, D_jet0, eqn.theta0};
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::dN3dtPerOmega
 * DESCRIPTION: Computes the derivative of N3 (number per unit solid angle) with respect to time.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
Real FRShockEqn<Jet, Injector>::dN3dtPerOmega(Real r, Real n1, Real n4, Real gamma3, Real drdt) {
    Real gamma34 = relativeLorentz(gamma4, gamma3);
    Real ratio_u = u_UpStr2u_DownStr(gamma34, this->jet_sigma);
    Real n3 = n4 * ratio_u;
    Real dxdr = 1. / (gamma4 * std::sqrt((1 + this->jet_sigma) * n4 / n1) * std::abs(1 - gamma4 * n4 / gamma3 / n3));
    return n3 * r * r * gamma3 * dxdr * drdt;
}

/********************************************************************************************************************
 * CONSTRUCTOR: FRShockEqn::FRShockEqn
 * DESCRIPTION: Initializes an FRShockEqn object with references to the medium, jet, and injector, and sets the
 *              angular coordinates, jet magnetization, and initial Gamma.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
FRShockEqn<Jet, Injector>::FRShockEqn(Medium const& medium, Jet const& jet, Injector const& inject, Real phi,
                                      Real theta)
    : medium(medium),
      jet(jet),
      inject(inject),
      phi(phi),
      theta0(theta),
      jet_sigma(jet.sigma0(phi, theta, 0)),
      gamma4(jet.Gamma0(phi, theta, 0)) {}

/********************************************************************************************************************
 * METHOD: FRShockEqn::operator()(State const& y, State& dydt, Real t)
 * DESCRIPTION: Computes the derivatives for the reverse shock evolution.
 *              The state vector for FRShockEqn is similar to that of ForwardShockEqn.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
void FRShockEqn<Jet, Injector>::operator()(State const& y, State& dydt, Real t) {
    // y[0] is left blank.
    // Real N_per_Omega = y[1];
    Real r = y[2];
    // Real t_com = y[3];
    Real D_jet_lab = y[4];
    Real theta = y[5];

    Real n4 = calc_n4(jet.dEdOmega(phi, theta0, t), gamma4, r, D_jet_lab, jet_sigma);
    Real n1 = medium.rho(r) / con::mp;

    Real gamma3 = calc_gamma3(n1, n4, gamma4, jet_sigma, medium.k);
    Real beta3 = gammaTobeta(gamma3);
    Real uv = gamma3 * beta3;
    Real beta4 = gammaTobeta(gamma4);
    dydt[0] = 0;
    dydt[2] = drdt(beta3);
    dydt[1] = dN3dtPerOmega(r, n1, n4, gamma3, dydt[2]);
    dydt[3] = dtdt_CoMoving(gamma3, beta3);
    dydt[4] = dDdt_Jet(gamma4, beta4);
    if (jet.spreading && theta < 0.5 * con::pi && uv * theta < 0.5) {
        dydt[5] = dtheta_dt(uv, dydt[2], r, gamma3);
    } else {
        dydt[5] = 0;
    }
}

inline Real calcFromRelativeLF(Real gamma34, Real gamma4) {
    Real b = 2 * gamma34 * gamma4;
    Real c = gamma34 * gamma34 + gamma4 * gamma4 - 1;
    return (b - std::sqrt(b * b - 4 * c)) / 2;
}

/********************************************************************************************************************
 * INLINE FUNCTION: postCrossingProfile
 * DESCRIPTION: Updates the shock state post shock crossing
 ********************************************************************************************************************/
template <typename State>
inline void postCrossingProfile(size_t i, size_t j, size_t k, size_t k0, State const& state, Shock& shock, Real ) {
    Real r = state[2];
    Real t_com = state[3];
    Real theta = state[4];

    Real const g = 2.;
    Real r0 = shock.r[i][j][k0];
    // Real Gamma0 = calcFromRelativeLF(shock.Gamma_rel[i][j][k0], gamma4);
    // Real Gamma3 = (Gamma0 - 1) * std::pow(r / r0, -g) + 1;
    shock.t_com[i][j][k] = t_com;
    shock.r[i][j][k] = r;
    shock.theta[i][j][k] = theta;
    shock.Gamma_rel[i][j][k] = (shock.Gamma_rel[i][j][k0] - 1) * std::pow(r / r0, -(6 + 2 * g) / 7) + 1;
    shock.column_num_den[i][j][k] = shock.column_num_den[i][j][k0] * (r0 * r0) / (r * r);
    shock.B[i][j][k] = shock.B[i][j][k0] * std::pow(r / r0, (6. * (3. + g) - 14.) / 7.);
}

// Determines whether a reverse shock exists based on current shock parameters.
template <typename ShockEqn>
bool reverseShockExists(ShockEqn const& eqn, Real t0) {
    Real gamma = eqn.jet.Gamma0(eqn.phi, eqn.theta0, t0);  // Initial Lorentz factor
    Real beta0 = gammaTobeta(gamma);
    Real r0 = beta0 * con::c * t0 / (1 - beta0);
    Real D_jet0 = con::c * eqn.jet.duration;

    Real n4 = calc_n4(eqn.jet.dEdOmega(eqn.phi, eqn.theta0, t0), gamma, r0, D_jet0, eqn.jet_sigma);
    Real n1 = eqn.medium.rho(r0) / con::mp;
    return eqn.jet_sigma < 8. / 3 * gamma * gamma * n1 / n4;
}

template <typename State>
inline auto unpack_r_state(State const& r_state) {
    return std::make_tuple(r_state[0], r_state[1], r_state[2], r_state[3], r_state[4], r_state[5]);
}

template <typename ShockEqn, typename FState, typename RState>
void set_f_state(ShockEqn& eqn, FState& f_state, RState const& r_state, Real gamma2) {
    Real r0 = r_state[2];
    Real theta = r_state[5];
    Real u0 = (gamma2 - 1) * eqn.medium.mass(r0) / (4 * con::pi) * con::c2;
    Real t_com0 = r_state[3];
    f_state = {gamma2, u0, r0, t_com0, theta};
    eqn.gamma4 = gamma2;
}

template <typename ShockEqn, typename State>
bool updateForwardReverseShock(size_t i, size_t j, int k, ShockEqn& eqn, State const& state, Shock& f_shock,
                               Shock& r_shock, Real dEdOmega, Real dN4dOmega) {
    auto [ignore, dN3dOmega, r, t_com, D_jet, theta] = unpack_r_state(state);

    Real n1 = eqn.medium.rho(r) / con::mp;
    Real n4 = calc_n4(dEdOmega, eqn.gamma4, r, D_jet, eqn.jet_sigma);
    Real gamma3 = calc_gamma3(n1, n4, eqn.gamma4, eqn.jet_sigma, eqn.medium.k);
    Real dN2dOmega = eqn.medium.mass(r) / (4 * con::pi * con::mp);

    updateShockState(f_shock, i, j, k, r, theta, gamma3, t_com, dN2dOmega, n1, eqn.jet_sigma);

    Real gamma34 = relativeLorentz(eqn.gamma4, gamma3);
    updateShockState(r_shock, i, j, k, r, theta, gamma34, t_com, dN3dOmega, n4, eqn.jet_sigma);
    return dN3dOmega >= dN4dOmega;
}

/********************************************************************************************************************
 * FUNCTION: solveFRShell
 * DESCRIPTION: Solves the evolution of the forward–reverse shock shell over the radius array r, updating both
 *              the forward shock (f_shock) and reverse shock (r_shock) objects accordingly.
 ********************************************************************************************************************/
template <typename FShockEqn, typename RShockEqn>
void solveFRShell(size_t i, size_t j, Array const& t, Shock& f_shock, Shock& r_shock, FShockEqn& eqn_f,
                  RShockEqn& eqn_r, Real rtol = 1e-9) {
    using namespace boost::numeric::odeint;

    if (!reverseShockExists(eqn_r, t[0])) {
        solveForwardShell(i, j, t, f_shock, eqn_f, rtol);
    } else {
        Real t0 = t[0];
        Real dt = (t[1] - t[0]) / 100;
        typename RShockEqn::State r_state;
        typename FShockEqn::State f_state;
        setReverseInit(eqn_r, r_state, t0);

        if (eqn_r.gamma4 < con::Gamma_cut) {
            setStoppingShock(i, j, f_shock, t, r_state[2], r_state[5]);
            setStoppingShock(i, j, r_shock, t, r_state[2], r_state[5]);
            return;
        }

        Real dEdOmega = eqn_r.jet.dEdOmega(eqn_r.phi, eqn_r.theta0, t0);
        Real dN4dOmega = dEdOmega / (eqn_r.gamma4 * con::mp * con::c2 * (1 + eqn_r.jet_sigma));
        bool crossed = false;

        // auto r_stepper = bulirsch_stoer_dense_out<typename RShockEqn::State>{0, rtol};
        auto r_stepper = make_dense_output(0, rtol, runge_kutta_dopri5<typename RShockEqn::State>());

        r_stepper.initialize(r_state, t0, dt);

        Real t_back = t[t.size() - 1];

        int k = 0, k0 = 0;
        for (; !crossed && r_stepper.current_time() <= t_back;) {
            r_stepper.do_step(eqn_r);
            for (; k < t.size() && r_stepper.current_time() > t[k]; k++) {
                r_stepper.calc_state(t[k], r_state);
                crossed = updateForwardReverseShock(i, j, k, eqn_r, r_state, f_shock, r_shock, dEdOmega, dN4dOmega);
                if (crossed) {
                    k0 = k;  // k0 is the index of the first element in the array that the reverse shock crosses
                    break;
                }
            }
        }

        // auto f_stepper = bulirsch_stoer_dense_out<typename FShockEqn::State>{0, rtol};
        auto f_stepper = make_dense_output(0, rtol, runge_kutta_dopri5<typename FShockEqn::State>());
        set_f_state(eqn_f, f_state, r_state, f_shock.Gamma_rel[i][j][k0]);
        f_stepper.initialize(f_state, t[k], r_stepper.current_time_step());

        for (; f_stepper.current_time() <= t_back;) {
            f_stepper.do_step(eqn_f);

            for (; k < t.size() && f_stepper.current_time() > t[k]; k++) {
                f_stepper.calc_state(t[k], f_state);
                updateForwardShock(i, j, k, eqn_f, f_state, f_shock);
                postCrossingProfile(i, j, k, k0, f_state, r_shock);
            }
        }
    }
}

#endif