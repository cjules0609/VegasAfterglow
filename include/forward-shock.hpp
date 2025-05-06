//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <array>

#include "shock.h"

/********************************************************************************************************************
 * @struct ForwardState
 * @brief State vector structure for forward shock calculations.
 * @details Template structure that defines the state vector for forward shock calculations, adapting its size
 *          based on whether the Ejecta class supports mass and energy injection methods.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
struct ForwardState {
    static constexpr bool mass_inject = HasDmdt<Ejecta>::value;    // whether Ejecta class has dmdt method
    static constexpr bool energy_inject = HasDedt<Ejecta>::value;  // whether Ejecta class has dedt method
    // use least fixed array size for integrator efficiency
    static constexpr size_t array_size = 5 + (mass_inject ? 1 : 0) + (energy_inject ? 1 : 0);

    MAKE_THIS_ODEINT_STATE(ForwardState, data, array_size)

    union {
        struct {
            Real Gamma;   // Lorentz factor
            Real u;       // internal energy density
            Real r;       // radius
            Real t_comv;  // comoving time
            Real theta;   // angle

            // shell energy density per solid angle
            [[no_unique_address]] std::conditional_t<energy_inject, Real, class Empty> eps_shell;

            // shell mass per solid angle
            [[no_unique_address]] std::conditional_t<mass_inject, Real, class Empty> m_shell;
        };
        array_type data;
    };
};

/********************************************************************************************************************
 * @class ForwardShockEqn
 * @brief Represents the forward shock equation for a given jet and medium.
 * @details It defines a state vector (with variable size based on template parameters) and overloads operator()
 *          to compute the derivatives of the state with respect to t. It also declares helper functions for the
 *          derivatives. This class implements the physical equations governing the forward shock evolution.
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

    // Set initial state for the ODE solver
    void set_init_state(State& state, Real t0) const noexcept;

    // References to model components
    Medium const& medium;  // Reference to the ambient medium properties
    Ejecta const& ejecta;  // Reference to the ejecta properties

    // Model parameters
    Real const phi{0};     // Angular coordinate phi in the jet frame
    Real const theta0{0};  // Initial angular coordinate theta
    Real const eps_e{0};   // Fraction of energy given to electrons

   private:
    // Computes the time derivative of the Lorentz factor with respect to engine time
    inline Real dGamma_dt(Real m_swept, Real dm_dt_swept, State const& state, State const& diff,
                          Real ad_idx) const noexcept;

    // Computes the time derivative of internal energy with respect to engine time
    inline Real dU_dt(Real m_swept, Real dm_dt_swept, State const& state, State const& diff,
                      Real ad_idx) const noexcept;

    // Private member variables
    Real const dOmega0{0};  // Initial solid angle element
    Real const theta_s{0};  // Jet structure parameter controlling angular dependence
    Real m_shell{0};        // Ejecta mass per solid angle
};

#include "../src/dynamics/forward-shock.tpp"