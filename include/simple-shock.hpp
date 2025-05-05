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
 * CLASS: SimpleState
 * DESCRIPTION: Represents the state vector for the simple shock equation.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
struct SimpleState {
    static constexpr bool mass_inject = HasDmdt<Ejecta>;    // whether Ejecta class has dmdt method
    static constexpr bool energy_inject = HasDedt<Ejecta>;  // whether Ejecta class has dedt method
    // use least fixed array size for integrator efficiency
    static constexpr size_t array_size = 4 + (mass_inject ? 1 : 0) + (energy_inject ? 1 : 0);

    MAKE_THIS_ODEINT_STATE(SimpleState, data, array_size)

    union {
        struct {
            Real Gamma{1};   // Lorentz factor
            Real r{0};       // radius
            Real t_comv{0};  // comoving time
            Real theta{0};   // angle

            // shell energy density per solid angle
            [[no_unique_address]] std::conditional_t<energy_inject, Real, class Empty> eps_shell;

            // shell mass per solid angle
            [[no_unique_address]] std::conditional_t<mass_inject, Real, class Empty> m_shell{};
        };
        array_type data;
    };
};

/********************************************************************************************************************
 * CLASS: SimpleShockEqn
 * DESCRIPTION: Represents the forward shock equation for a given Jet. It defines a state vector
 *              and overloads operator() to compute the derivatives of the state with
 *              respect to time t. It also declares helper functions for the derivatives. Simple version from
 *              Huang et al. 2000
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
class SimpleShockEqn {
   public:
    using State = SimpleState<Ejecta, Medium>;

    SimpleShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, Real eps_e, Real theta_s);

    // Overloaded operator() to compute the derivatives of the state vector with respect to engine time t.
    void operator()(State const& state, State& diff, Real t) const noexcept;

    // Initialize the state vector at time t0
    void set_init_state(State& state, Real t0) const noexcept;

    Medium const& medium;  // Reference to the medium properties
    Ejecta const& ejecta;  // Reference to the ejecta properties
    Real const phi{0};     // Angular coordinate phi
    Real const theta0{0};  // Angular coordinate theta
    Real const eps_e{0};   // Electron energy fraction

   private:
    // Helper function: computes the derivative of Gamma with respect to engine time t.
    Real dGamma_dt(Real dm_dt_swept, State const& state, State const& diff) const noexcept;
    Real const dOmega0{0};  // Initial solid angle
    Real const theta_s{0};  // Critical angle for jet spreading
    Real m_shell{0};        // Ejecta mass per solid angle
};

#include "../src/dynamics/simple-shock.tpp"