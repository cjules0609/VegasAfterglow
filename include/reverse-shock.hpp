//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <array>

#include "macros.h"
/**
 * <!-- ************************************************************************************** -->
 * @struct ReverseState
 * @brief Represents the state variables for the reverse shock simulation.
 * @details It defines a state vector containing properties like shell width, mass, radius, time, and energy.
 *          The struct dynamically adapts its size based on template parameters to include mass/energy
 *          injection capabilities.
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium>
struct ReverseState {
    static constexpr bool mass_inject = HasDmdt<Ejecta>;    ///< Whether ejecta has mass injection
    static constexpr bool energy_inject = HasDedt<Ejecta>;  ///< Whether ejecta has energy injection
    static constexpr size_t array_size = 7;

    MAKE_THIS_ODEINT_STATE(ReverseState, data, array_size)

    union {
        struct {
            Real width_shell;  ///< Width of the shell
            Real m3;           ///< Shocked ejecta mass per solid angle
            Real r;            ///< Radius
            Real t_comv;       ///< Comoving time
            Real theta;        ///< Angular coordinate theta
            Real eps_shell;    ///< energy of shell per solid angle
            Real m_shell;      ///< shell mass per solid angle
        };
        array_type data;
    };
};

/**
 * <!-- ************************************************************************************** -->
 * @class FRShockEqn
 * @brief Represents the reverse shock (or forward-reverse shock) equation for a given Jet and medium.
 * @details It defines a state vector (an array of 8 Reals) and overloads operator() to compute the
 *          derivatives of the state with respect to radius r.
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium>
class FRShockEqn {
   public:
    using State = ReverseState<Ejecta, Medium>;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Constructor for the FRShockEqn class.
     * @details Initializes the forward-reverse shock equation with the given medium, ejecta, and parameters.
     * @param medium The medium through which the shock propagates
     * @param ejecta The ejecta driving the shock
     * @param phi Azimuthal angle
     * @param theta Polar angle
     * <!-- ************************************************************************************** -->
     */
    FRShockEqn(Medium const& medium, Ejecta const& jet, Real phi, Real theta, RadParams const& rad_params);

    Medium const& medium;  ///< Reference to the medium properties
    Ejecta const& ejecta;  ///< Reference to the jet properties
    RadParams const rad;   ///< Radiation parameters
    Real const phi{0};     ///< Angular coordinate phi
    Real const theta0{0};  ///< Angular coordinate theta
    Real Gamma4{1};        ///< Initial Lorentz factor of the jet
    Real u_x{0};           ///< Reverse shock crossed four velocity
    Real r_x{0};           ///< Reverse shock crossed radius

    /**
     * <!-- ************************************************************************************** -->
     * @brief Implements the reverse shock ODE system.
     * @details Computes the derivatives of state variables with respect to time.
     * @param state Current state of the system
     * @param diff Output derivatives to be populated
     * @param t Current time
     * <!-- ************************************************************************************** -->
     */
    void operator()(State const& state, State& diff, Real t);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Set the initial conditions for the reverse shock ODE.
     * @details Sets up initial state values and determines if the shock has already crossed.
     * @param state State vector to initialize
     * @param t0 Initial time
     * @return True if the shock has already crossed at the initial time, false otherwise
     * <!-- ************************************************************************************** -->
     */
    bool set_init_state(State& state, Real t0) const noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Sets constants at the shock crossing point for later calculations.
     * @details Stores key parameters and switches the ODE from crossing to crossed state.
     * @param state State at shock crossing
     * @param B Magnetic field at shock crossing
     * <!-- ************************************************************************************** -->
     */
    void set_cross_state(State const& state, Real B);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the Lorentz factor (Gamma3) during shock crossing.
     * @details Uses energy conservation to determine the appropriate Gamma3 value.
     * @param state Current state of the system
     * @return The computed Lorentz factor for region 3 (shocked ejecta)
     * <!-- ************************************************************************************** -->
     */
    Real compute_crossing_Gamma3(State const& state) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the relative Lorentz factor after shock crossing.
     * @details Uses adiabatic expansion to determine the relative Lorentz factor.
     * @param state Current state of the system
     * @return The relative Lorentz factor post-crossing
     * <!-- ************************************************************************************** -->
     */
    Real compute_crossed_Gamma_rel(State const& state) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the magnetic field after shock crossing.
     * @details Calculates magnetic field based on pressure and energy conservation.
     * @param state Current state of the system
     * @return Magnetic field strength in the shocked region
     * <!-- ************************************************************************************** -->
     */
    Real compute_crossed_B(State const& state) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the Lorentz factor for region 3 after shock crossing.
     * @details Uses a power-law profile that applies both in relativistic and Newtonian regimes.
     * @param gamma_rel Relative Lorentz factor
     * @param r Current radius
     * @return The Lorentz factor for region 3 post-crossing
     * <!-- ************************************************************************************** -->
     */
    Real compute_crossed_Gamma3(Real Gamma_rel, Real r) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the magnetization parameter of the shell.
     * @details Sigma is defined as (ε/Γmc²) - 1, where ε is the energy per solid angle.
     * @param state Current state of the system
     * @return The magnetization parameter of the shell
     * <!-- ************************************************************************************** -->
     */
    Real compute_shell_sigma(State const& state) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Checks if the ejecta is still injecting mass or energy at time t.
     * @details Returns true if either energy or mass is being injected.
     * @param t Current time
     * @return Boolean indicating if injection is happening at time t
     * <!-- ************************************************************************************** -->
     */
    bool is_injecting(Real t) const;

   private:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the energy and mass injection rates at time t.
     * @details Returns a pair containing energy injection rate and mass injection rate.
     * @param t Current time
     * @return A pair {energy_rate, mass_rate} with injection rates at time t
     * <!-- ************************************************************************************** -->
     */
    std::pair<Real, Real> get_injection_rate(Real t) const;

    Real N_electron{0};        ///< Normalized total electron (for post crossing scaling calculation)
    Real adiabatic_const{1};   ///< Normalized adiabatic constant where C = rho^idx/p
    Real e_mag_const{1};       ///< Normalized magnetic energy constant where C = B^2/p
    Real gamma_hat_x{4. / 3};  ///< Adiabatic index at the shock crossing
    Real deps0_dt{0};          ///< Ejecta energy injection rate
    Real dm0_dt{0};            ///< Ejecta mass injection rate
    Real u4{0};                ///< Four-velocity of the unshocked ejecta
    bool crossed{false};       ///< Flag indicating if shock has crossed the shell
};

/**
 * <!-- ************************************************************************************** -->
 * @brief Function templates for shock generation
 * @details Declares interfaces to generate forward shocks (2D and 3D) and forward/reverse shock pairs.
 * <!-- ************************************************************************************** -->
 */
using ShockPair = std::pair<Shock, Shock>;

/**
 * <!-- ************************************************************************************** -->
 * @brief Generates a pair of forward and reverse shocks.
 * @details Creates two Shock objects and solves the shock evolution for each grid point.
 * @param coord Coordinate system definition
 * @param medium The medium through which the shock propagates
 * @param jet The jet (ejecta) driving the shock
 * @param eps_e_f Electron energy fraction for forward shock
 * @param eps_B_f Magnetic energy fraction for forward shock
 * @param eps_e_r Electron energy fraction for reverse shock
 * @param eps_B_r Magnetic energy fraction for reverse shock
 * @param rtol Relative tolerance for ODE solver
 * @return A pair of Shock objects {forward_shock, reverse_shock}
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium>
ShockPair generate_shock_pair(Coord const& coord, Medium const& medium, Ejecta const& jet, RadParams const& rad_fwd,
                              RadParams const& rad_rvs, Real rtol = 1e-5);

//========================================================================================================
//                                  template function implementation
//========================================================================================================

#include "../src/dynamics/reverse-shock.tpp"