//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <array>

#include "shock.h"

template <typename Ejecta, typename Medium>
struct ForwardState {
    static constexpr bool mass_inject = HasDmdt<Ejecta>;
    static constexpr bool energy_inject = HasDedt<Ejecta>;
    static constexpr bool mass_profile = HasMass<Medium>;
    static constexpr size_t array_size = 5 + (mass_inject ? 1 : 0) + (energy_inject ? 1 : 0) + (mass_profile ? 0 : 1);

    MAKE_THIS_ODEINT_STATE(data, array_size)

    union {
        struct {
            Real Gamma{1};
            Real u{0};
            Real r{0};
            Real t_com{0};
            Real theta{0};

            [[no_unique_address]] std::conditional_t<energy_inject, Real, class Empty> E_ej{};
            [[no_unique_address]] std::conditional_t<mass_inject, Real, class Empty> M_ej{};
            [[no_unique_address]] std::conditional_t<mass_profile, class Empty, Real> M_sw{};
        };
        array_type data;
    };
};

/********************************************************************************************************************
 * CLASS: ForwardShockEqn
 * DESCRIPTION: Represents the forward shock equation for a given jet and medium. It defines a state vector
 *              (an array of 8 Reals) and overloads operator() to compute the derivatives of the state with
 *              respect to t. It also declares helper functions for the derivatives.
 *              This class implements the physical equations governing the forward shock evolution.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
class ForwardShockEqn {
   public:
    using State = ForwardState<Ejecta, Medium>;

    // Constructor: Initialize the forward shock equation with physical parameters
    ForwardShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, Real eps_e, Real theta_s);

    // Forward shock ODE equation - callable interface for ODE solver
    // Computes the time derivatives of the state vector
    void operator()(State const& state, State& diff, Real t) const noexcept;

    void setInitState(State& state, Real t0) const noexcept;

    // References to model components
    Medium const& medium;  // Reference to the ambient medium properties
    Ejecta const& ejecta;  // Reference to the ejecta properties

    // Model parameters
    Real const phi{0};     // Angular coordinate phi in the jet frame
    Real const theta0{0};  // Initial angular coordinate theta
    Real const eps_e{0};   // Fraction of energy given to electrons

   private:
    // Computes the time derivative of the Lorentz factor
    // with respect to on-axis observer time
    inline Real dGammadt(Real M_sw, Real dMdt_sw, State const& state, State const& diff, Real ad_idx) const noexcept;

    // Computes the time derivative of internal energy
    // with respect to on-axis observer time
    inline Real dUdt(Real M_sw, Real dMdt_sw, State const& state, State const& diff, Real ad_idx) const noexcept;

    // Private member variables
    Real const dOmega0{0};  // Initial solid angle element
    Real const theta_s{0};  // Jet structure parameter controlling angular dependence
    Real M_ej{0};           // Ejecta mass per solid angle
};

#include "../src/dynamics/forward-shock.tpp"