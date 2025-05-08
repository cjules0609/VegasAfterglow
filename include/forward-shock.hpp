//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <array>

#include "shock.h"

/**
 * <!-- ************************************************************************************** -->
 * @struct ForwardState
 * @brief State vector structure for forward shock calculations.
 * @details Template structure that defines the state vector for forward shock calculations, adapting its size
 *          based on whether the Ejecta class supports mass and energy injection methods.
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium>
struct ForwardState {
    static constexpr bool mass_inject = HasDmdt<Ejecta>;    ///< whether Ejecta class has dmdt method
    static constexpr bool energy_inject = HasDedt<Ejecta>;  ///< whether Ejecta class has dedt method
    /// use least fixed array size for integrator efficiency
    static constexpr size_t array_size = 5 + (mass_inject ? 1 : 0) + (energy_inject ? 1 : 0);

    MAKE_THIS_ODEINT_STATE(ForwardState, data, array_size)

    union {
        struct {
            Real Gamma;   ///< Lorentz factor
            Real u;       ///< internal energy density
            Real r;       ///< radius
            Real t_comv;  ///< comoving time
            Real theta;   ///< angle

            /// shell energy density per solid angle
            [[no_unique_address]] std::conditional_t<energy_inject, Real, class Empty> eps_shell;

            /// shell mass per solid angle
            [[no_unique_address]] std::conditional_t<mass_inject, Real, class Empty> m_shell;
        };
        array_type data;
    };
};

/**
 * <!-- ************************************************************************************** -->
 * @class ForwardShockEqn
 * @brief Represents the forward shock equation for a given jet and medium.
 * @details It defines a state vector (with variable size based on template parameters) and overloads operator()
 *          to compute the derivatives of the state with respect to t. It also declares helper functions for the
 *          derivatives. This class implements the physical equations governing the forward shock evolution.
 *
 * @tparam Ejecta The ejecta class template parameter
 * @tparam Medium The medium class template parameter
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium>
class ForwardShockEqn {
   public:
    using State = ForwardState<Ejecta, Medium>;

    /**
     * <!-- ************************************************************************************** -->
     * @brief ForwardShockEqn constructor
     * @details Initializes the forward shock equation with the given medium, ejecta, and parameters.
     * @param medium The medium through which the shock propagates
     * @param ejecta The ejecta driving the shock
     * @param phi Azimuthal angle
     * @param theta Polar angle
     * @param eps_e Electron energy fraction
     * @param theta_s Critical angle for jet spreading
     * <!-- ************************************************************************************** -->
     */
    ForwardShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, Real eps_e, Real theta_s);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the derivatives of the state variables with respect to engine time t.
     * @details Implements the system of ODEs that describe the evolution of the forward shock.
     * @param state Current state of the system
     * @param diff Output derivatives to be populated
     * @param t Current time
     * <!-- ************************************************************************************** -->
     */
    void operator()(State const& state, State& diff, Real t) const noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Set the initial conditions for the forward shock ODE solver.
     * @details Computes the initial state of the shock based on the ejecta properties and ambient medium.
     * @param state State vector to initialize
     * @param t0 Initial time
     * <!-- ************************************************************************************** -->
     */
    void set_init_state(State& state, Real t0) const noexcept;

    Medium const& medium;  ///< Reference to the ambient medium properties
    Ejecta const& ejecta;  ///< Reference to the ejecta properties

    Real const phi{0};     ///< Angular coordinate phi in the jet frame
    Real const theta0{0};  ///< Initial angular coordinate theta
    Real const eps_e{0};   ///< Fraction of energy given to electrons

   private:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the derivative of Gamma with respect to engine time t.
     * @details Calculates the rate of change of the Lorentz factor based on various physical factors.
     * @param m_swept Total swept-up mass
     * @param dm_dt_swept Rate of swept-up mass
     * @param state Current state of the system
     * @param diff Current derivatives
     * @param ad_idx Adiabatic index
     * @return The time derivative of Gamma
     * <!-- ************************************************************************************** -->
     */
    inline Real dGamma_dt(Real m_swept, Real dm_dt_swept, State const& state, State const& diff,
                          Real ad_idx) const noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the derivative of internal energy with respect to time t.
     * @details Calculates the rate of change of the internal energy density considering adiabatic expansion
     *          and energy injection from newly swept-up material.
     * @param m_swept Total swept-up mass
     * @param dm_dt_swept Rate of swept-up mass
     * @param state Current state of the system
     * @param diff Current derivatives
     * @param ad_idx Adiabatic index
     * @return The time derivative of internal energy
     * <!-- ************************************************************************************** -->
     */
    inline Real dU_dt(Real m_swept, Real dm_dt_swept, State const& state, State const& diff,
                      Real ad_idx) const noexcept;

    Real const dOmega0{0};  ///< Initial solid angle element
    Real const theta_s{0};  ///< Jet structure parameter controlling angular dependence
    Real m_shell{0};        ///< Ejecta mass per solid angle
};

/**
 * <!-- ************************************************************************************** -->
 * @brief Updates the forward shock state at a grid point using the current ODE solution.
 * @details Computes physical quantities for the shock state and updates the shock object.
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param k Grid index for time
 * @param eqn The equation system containing physical parameters
 * @param state Current state of the system
 * @param shock Shock object to update
 * <!-- ************************************************************************************** -->
 */
template <typename Eqn, typename State>
void save_fwd_shock_state(size_t i, size_t j, size_t k, Eqn const& eqn, State const& state, Shock& shock);

/**
 * <!-- ************************************************************************************** -->
 * @brief Solves the forward shock ODE at a grid point as a function of time.
 * @details Uses an adaptive step size ODE solver to evolve the shock state through time.
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param t View of time points at which to evaluate the solution
 * @param shock Shock object to store the results
 * @param eqn The equation system defining the ODE
 * @param rtol Relative tolerance for the ODE solver
 * <!-- ************************************************************************************** -->
 */
template <typename FwdEqn, typename View>
void grid_solve_fwd_shock(size_t i, size_t j, View const& t, Shock& shock, FwdEqn const& eqn, double rtol);

/**
 * <!-- ************************************************************************************** -->
 * @brief Generates a forward shock model for given coordinates, medium, and jet parameters.
 * @details Creates a Shock object and solves the shock evolution for each grid point.
 * @param coord Coordinate system definition
 * @param medium The medium through which the shock propagates
 * @param jet The jet (ejecta) driving the shock
 * @param eps_e Electron energy fraction
 * @param eps_B Magnetic energy fraction
 * @param rtol Relative tolerance for the ODE solver (default: 1e-5)
 * @return A Shock object containing the evolution data
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium>
Shock generate_fwd_shock(Coord const& coord, Medium const& medium, Ejecta const& jet, Real eps_e, Real eps_B,
                         Real rtol = 1e-5);

//========================================================================================================
//                                  template function implementation
//========================================================================================================
#include "../src/dynamics/forward-shock.tpp"