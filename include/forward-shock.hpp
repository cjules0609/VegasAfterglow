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
    using State = std::array<Real, 5>;  // State vector: typically [Gamma, u, t_eng, t_com, D_jet]

    ForwardShockEqn(Medium const& medium, Jet const& jet, Injector const& inject, Real phi, Real theta, Real eps_e);

    Medium const& medium;      // Reference to the medium properties
    Jet const& jet;            // Reference to the jet properties
    Injector const& inject;    // Reference to the injector properties
    Real const phi{0};         // Angular coordinate phi
    Real const theta{0};       // Angular coordinate theta
    Real const eps_e{0};       // Electron energy fraction
    Real const jet_sigma{0};   // Jet magnetization parameter
    Real gamma4{1};            // Initial Lorentz factor (or a related parameter)
    Real spreading_factor{1};  // Factor to account for jet spreading

    // Overloaded operator() to compute the derivatives of the state vector with respect to radius r.
    void operator()(State const& y, State& dydr, Real r);

   private:
    // Helper function: computes the derivative of Gamma with respect to r.
    inline Real dGammadr(Real r, Real Gamma, Real u, Real t_eng, Real ad_idx, Real rho, Real dtdr);
    // Helper function: computes the derivative of u with respect to r.
    inline Real dUdr(Real r, Real Gamma, Real u, Real t_eng, Real ad_idx, Real rho, Real dGdr);

    Real const inj_Gamma0{0};  // Initial Gamma from the injector
    Real const inj_sigma{0};   // Injector magnetization parameter
    Real const dM0{0};         // Initial mass per unit solid angle
};

/********************************************************************************************************************
 * METHOD: ForwardShockEqn::operator()(State const& y, State& dydr, Real r)
 * DESCRIPTION: Computes the derivatives of the state variables with respect to radius r.
 *              The state vector components are:
 *                  y[0] - Gamma (Lorentz factor)
 *                  y[1] - u (internal energy-related variable)
 *                  y[2] - t_eng (engine time)
 *                  y[3] - t_com (co-moving time) [unused here]
 *                  y[4] - D_jet (jet shell width) [unused here]
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
void ForwardShockEqn<Jet, Injector>::operator()(State const& y, State& dydr, Real r) {
    Real Gamma = y[0];
    Real u = y[1];
    Real t_eng = y[2];  // engine time
    // Real t_com = y[3];  // co-moving time (unused)
    // Real D_jet = y[4];  // co-moving jet shell width (unused)

    Real ad_idx = adiabaticIndex(Gamma);  // Compute adiabatic index based on Gamma
    Real rho = medium.rho(r);             // Get medium density at radius r
    Real beta = gammaTobeta(Gamma);       // Convert Gamma to beta (velocity/c)
    Real beta4 = gammaTobeta(gamma4);     // Convert gamma4 to beta

    dydr[2] = dtdr_Engine(beta);  // Compute derivative of engine time with respect to r
    dydr[0] = dGammadr(r, Gamma, u, t_eng, ad_idx, rho, dydr[2]);  // d(Gamma)/dr
    dydr[1] = dUdr(r, Gamma, u, t_eng, ad_idx, rho, dydr[0]);      // d(u)/dr
    dydr[3] = dtdr_CoMoving(Gamma, beta);                          // d(t_com)/dr
    dydr[4] = dDdr_Jet(gamma4, beta4);                             // d(D_jet)/dr
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
      theta(theta),
      eps_e(eps_e),
      jet_sigma(jet.sigma0(phi, theta, 0)),
      gamma4(jet.Gamma0(phi, theta, 0)),
      spreading_factor(1),
      inj_Gamma0(inject.Gamma0(phi, theta, 0)),
      inj_sigma(inject.sigma0(phi, theta, 0)),
      dM0(jet.dEdOmega(phi, theta, 0) / (gamma4 * (1 + jet_sigma) * con::c2)) {
    // dM0dOmega(jet.dE0dOmega(theta) / (jet.Gamma0(theta) * con::c2)) is commented out.
}

/********************************************************************************************************************
 * METHOD: ForwardShockEqn::dGammadr
 * DESCRIPTION: Computes the derivative of Gamma with respect to radius r.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
Real ForwardShockEqn<Jet, Injector>::dGammadr(Real r, Real Gamma, Real u, Real t_eng, Real ad_idx, Real rho,
                                              Real dtdr) {
    Real ad_idx_m1 = ad_idx - 1;
    Real Gamma2_m1 = Gamma * Gamma - 1;
    Real term1 = ad_idx * Gamma2_m1 + 1;

    Real dm = medium.mass(r) / (4 * con::pi);  // Mass per unit solid angle from medium
    Real dm_inj = inject.dEdOmega(phi, theta, t_eng) / (inj_Gamma0 * (1 + inj_sigma) * con::c2);  // Injected mass
    Real L_inj = inject.dLdOmega(phi, theta, t_eng);  // Injected luminosity per unit solid angle

    Real a1 = -Gamma2_m1 * (ad_idx * Gamma - ad_idx + 1) * r * r * rho * con::c2;
    Real a2 = ad_idx_m1 * term1 * 3 * u / r;
    Real a3 = Gamma * dtdr * L_inj * (1 - Gamma / (inj_Gamma0 * (1 + inj_sigma)));

    Real b1 = Gamma * (dM0 + dm + dm_inj) * con::c2;
    Real b2 = (ad_idx * term1 + 2 * ad_idx_m1) / Gamma * u;

    return (a1 + a2 + a3) / (b1 + b2);
}

/********************************************************************************************************************
 * METHOD: ForwardShockEqn::dUdr
 * DESCRIPTION: Computes the derivative of u with respect to radius r.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
Real ForwardShockEqn<Jet, Injector>::dUdr(Real r, Real Gamma, Real u, Real t_eng, Real ad_idx, Real rho, Real dGdr) {
    Real E = r * r * rho * con::c2;
    return (1 - eps_e) * (Gamma - 1) * E - (ad_idx - 1) * (3 / r - dGdr / Gamma) * u * spreading_factor;
}

// Updates the forward shock state at grid index (i, j, k) using the current state vector and the medium properties.
template <typename ShockEqn>
void updateForwardShock(size_t i, size_t j, int k, Real r_k, ShockEqn& eqn, const typename ShockEqn::State& state,
                        Shock& f_shock) {
    Real n1 = eqn.medium.rho(r_k) / con::mp;                // Compute upstream number density
    Real Gamma = state[0];                                  // Lorentz factor from state
    Real t_eng = state[2];                                  // Engine time from state
    Real t_com = state[3];                                  // Comoving time from state
    Real dM1dOmega = eqn.medium.mass(r_k) / (4 * con::pi);  // Mass per unit solid angle

    updateShockState(f_shock, i, j, k, r_k, Gamma, t_com, t_eng, dM1dOmega, n1, eqn.jet_sigma);
}

// Initializes the forward shock state vector at radius r0.
template <typename ShockEqn>
void setForwardInit(ShockEqn& eqn, typename ShockEqn::State& state, Real r0) {
    Real gamma2 = eqn.jet.Gamma0(eqn.phi, eqn.theta, 0);  // Initial Lorentz factor
    Real u0 = (gamma2 - 1) * eqn.medium.mass(r0) / (4 * con::pi) * con::c2;
    Real beta0 = gammaTobeta(gamma2);
    Real t_eng0 = r0 * (1 - beta0) / beta0 / con::c;
    Real t_com0 = r0 / std::sqrt(gamma2 * gamma2 - 1) / con::c;
    Real D_jet0 = con::c * eqn.jet.duration;
    state = {gamma2, u0, t_eng0, t_com0, D_jet0};
}

// Solves the forward shock evolution for a given shell (across radius values in array r) and updates the Shock object.
template <typename ShockEqn>
void solveForwardShell(size_t i, size_t j, const Array& r, Shock& f_shock, ShockEqn& eqn, Real t_max) {
    using namespace boost::numeric::odeint;

    Real atol = 0, rtol = 1e-6, r0 = r[0];
    Real dr = (r[1] - r[0]) / 100;

    typename ShockEqn::State state;
    setForwardInit(eqn, state, r0);  // Initialize state at starting radius

    if (state[0] <= con::Gamma_cut) {  // If initial Lorentz factor is too low, exit early
        return;
    }

    // Create a dense output stepper for integrating the shock equations
    // auto stepper = bulirsch_stoer_dense_out<typename ShockEqn::State>{atol, rtol};

    auto stepper = make_dense_output(atol, rtol, runge_kutta_dopri5<typename ShockEqn::State>());
    stepper.initialize(state, r0, dr);

    Real r_back = r[r.size() - 1];  // Last radius in the array

    // Iterate over the radius array, updating the state and the Shock object as needed
    for (int k = 0; stepper.current_time() <= r_back && state[2] <= t_max;) {
        stepper.do_step(eqn);
        while (k < r.size() && stepper.current_time() > r[k]) {
            stepper.calc_state(r[k], state);
            // std::cout << "k: " << k << state[0] << std::endl;
            updateForwardShock(i, j, k, r[k], eqn, state, f_shock);
            ++k;
        }
    }
}
#endif