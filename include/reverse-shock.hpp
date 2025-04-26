//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <array>

template <typename Ejecta, typename Medium>
struct ReverseState {
    static constexpr bool mass_inject = HasDmdt<Ejecta>;
    static constexpr bool energy_inject = HasDedt<Ejecta>;
    static constexpr bool mass_profile = HasMass<Medium>;
    static constexpr size_t array_size = 7 + (mass_profile ? 0 : 1);

    MAKE_THIS_ODEINT_STATE(data, array_size)

    union {
        struct {
            Real width{0};
            Real M3{0};
            Real r{0};
            Real t_com{0};
            Real theta{0};
            Real E_ej{0};
            Real M_ej{0};

            [[no_unique_address]] std::conditional_t<mass_profile, class Empty, Real> M_sw{};
        };
        array_type data;
    };
};

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
    using State = ReverseState<Ejecta, Medium>;

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
    void operator()(State const& state, State& diff, Real t);

    bool setInitState(State& state, Real t0) const noexcept;

    // Set the shock state when the reverse shock crosses the jet.
    void setCrossState(State const& state, Real B, Real t);

    // calculate the Gamma3 during the shock crossing phase.
    Real crossingGamma3(State const& state) const;

    // calculate the Gamma_43 post shock crossing.
    Real crossedGamma_rel(State const& state) const;

    // calculate the magnetic field post shock crossing.
    Real crossedB(State const& state) const;

    // calculate the Gamma3 post shock crossing.
    Real crossedGamma3(Real gamma_rel, Real r) const;

    Real shellMagnetization(State const& state) const;

   private:
    std::pair<Real, Real> getInjection(Real t) const;

    bool isInjecting(Real t) const;

    Real N0{0};               // normalized total electron (for post crossing scaling calculation).
    Real adiabatic_const{1};  // normalized adiabatic constant where C = rho^idx/p.
    Real Emag_const{1};       // normalized magnetic energy constant where C = B^2/p.
    Real ad_idx0{4. / 3};     // adiabatic index at the shock crossing.
    Real dE0dt{0};            // ejecta energy injection rate
    Real dM0dt{0};            // ejecta mass injection rate
    Real u4{0};
    bool crossed{false};
};

#include "../src/dynamics/reverse-shock.tpp"