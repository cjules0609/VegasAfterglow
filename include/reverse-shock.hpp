//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <array>
/********************************************************************************************************************
 * CLASS: FRShockEqn
 * DESCRIPTION: Represents the reverse shock (or forward-reverse shock) equation for a given Jet and medium.
 *              It defines a state vector (an array of 8 Reals) and overloads operator() to compute the
 *              derivatives of the state with respect to radius r. It also declares a helper function to compute
 *              the derivative of N3 (number per solid angle) with respect to t.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
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
    using State = RState<StateArray>;             // Mutable state wrapper
    using constState = RState<const StateArray>;  // Const state wrapper

    FRShockEqn(Medium const& medium, Ejecta const& jet, Real phi, Real theta, Real eps_e);

    Medium const& medium;  // Reference to the medium properties
    Ejecta const& ejecta;  // Reference to the jet properties
    Real const phi{0};     // Angular coordinate phi
    Real const theta0{0};  // Angular coordinate theta
    Real const eps_e{0};   // Electron energy fraction
    Real gamma4{1};        // initial Lorentz factor of the jet
    Real u_x{0};           // reverse shock crossed four velocity
    Real r_x{0};           // reverse shock crossed radius

    // Reverse shock ODE equation
    void operator()(StateArray const& y, StateArray& dydt, Real t);

    // Set the shock state when the reverse shock crosses the jet.
    void setCrossState(State const& state, Real B, Real t);

    // calculate the Gamma3 during the shock crossing phase.
    Real crossingGamma3(State const& state, Real t) const;

    // calculate the Gamma_43 post shock crossing.
    template <typename GenState>
    Real crossedGamma_rel(GenState const& state) const;

    // calculate the magnetic field post shock crossing.
    Real crossedB(State const& state) const;

    // calculate the Gamma3 post shock crossing.
    Real crossedGamma3(Real gamma_rel, Real r) const;

   private:
    Real N0{0};               // normalized total electron (for post crossing scaling calculation).
    Real adiabatic_const{1};  // normalized adiabatic constant where C = rho^idx/p.
    Real Emag_const{1};       // normalized magnetic energy constant where C = B^2/p.
    Real ad_idx0{4. / 3};     // adiabatic index at the shock crossing.
    bool crossed{false};
};

#include "../src/dynamics/reverse-shock.tpp"