//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <array>

#include "macros.h"
/********************************************************************************************************************
 * STRUCT: ReverseState
 * DESCRIPTION: Represents the state variables for the reverse shock simulation. It defines a state vector
 *              containing properties like shell width, mass, radius, time, and energy. The struct dynamically
 *              adapts its size based on template parameters to include mass/energy injection capabilities.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
struct ReverseState {
    static constexpr bool mass_inject = HasDmdt<Ejecta>;    // Whether ejecta has mass injection
    static constexpr bool energy_inject = HasDedt<Ejecta>;  // Whether ejecta has energy injection
    static constexpr size_t array_size = 7;

    MAKE_THIS_ODEINT_STATE(ReverseState, data, array_size)

    union {
        struct {
            Real width_shell{0};  // Width of the shell
            Real m3{0};           // Shocked ejecta mass per solid angle
            Real r{0};            // Radius
            Real t_comv{0};       // Comoving time
            Real theta{0};        // Angular coordinate theta
            Real eps_shell{0};    // energy of shell per solid angle
            Real m_shell{0};      // shell mass per solid angle
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

    // Constructor: Initialize with medium, ejecta, angular coordinates, and electron energy fraction
    FRShockEqn(Medium const& medium, Ejecta const& jet, Real phi, Real theta, Real eps_e);

    Medium const& medium;  // Reference to the medium properties
    Ejecta const& ejecta;  // Reference to the jet properties
    Real const phi{0};     // Angular coordinate phi
    Real const theta0{0};  // Angular coordinate theta
    Real const eps_e{0};   // Electron energy fraction
    Real Gamma4{1};        // Initial Lorentz factor of the jet
    Real u_x{0};           // Reverse shock crossed four velocity
    Real r_x{0};           // Reverse shock crossed radius

    // Reverse shock ODE equation - callable interface for ODE solver
    void operator()(State const& state, State& diff, Real t);

    // Set initial state for the ODE solver
    bool set_init_state(State& state, Real t0) const noexcept;

    // Set the shock state when the reverse shock crosses the jet
    void set_cross_state(State const& state, Real B, Real t);

    // Calculate the Gamma3 during the shock crossing phase
    Real compute_crossing_Gamma3(State const& state) const;

    // Calculate the Gamma_rel (relative Lorentz factor) post shock crossing
    Real compute_crossed_Gamma_rel(State const& state) const;

    // Calculate the magnetic field post shock crossing
    Real compute_crossed_B(State const& state) const;

    // Calculate the Gamma3 post shock crossing
    Real compute_crossed_Gamma3(Real Gamma_rel, Real r) const;

    // Calculate the magnetization parameter of the shell
    Real compute_shell_sigma(State const& state) const;

    // Check if the shell is still being injected at time t
    bool is_injecting(Real t) const;

   private:
    // Get the energy and mass injection rates at time t
    std::pair<Real, Real> get_injection_rate(Real t) const;

    Real N_electron{0};        // Normalized total electron (for post crossing scaling calculation)
    Real adiabatic_const{1};   // Normalized adiabatic constant where C = rho^idx/p
    Real e_mag_const{1};       // Normalized magnetic energy constant where C = B^2/p
    Real gamma_hat_x{4. / 3};  // Adiabatic index at the shock crossing
    Real deps0_dt{0};          // Ejecta energy injection rate
    Real dm0_dt{0};            // Ejecta mass injection rate
    Real u4{0};                // Four-velocity of the unshocked ejecta
    bool crossed{false};       // Flag indicating if shock has crossed the shell
};

#include "../src/dynamics/reverse-shock.tpp"