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
    static constexpr bool mass_inject = HasDmdt<Ejecta>;
    static constexpr bool energy_inject = HasDedt<Ejecta>;
    static constexpr bool mass_profile = HasMass<Medium>;
    static constexpr size_t array_size = 4 + (mass_inject ? 1 : 0) + (energy_inject ? 1 : 0) + (mass_profile ? 0 : 1);

    using array_type = std::array<Real, array_size>;
    using value_type = typename array_type::value_type;
    using iterator = typename array_type::iterator;
    using const_iterator = typename array_type::const_iterator;

    union {
        struct {
            Real Gamma{1};
            Real r{0};
            Real t_com{0};
            Real theta{0};

            [[no_unique_address]] std::conditional_t<energy_inject, Real, class Empty> E_ej{};
            [[no_unique_address]] std::conditional_t<mass_inject, Real, class Empty> M_ej{};
            [[no_unique_address]] std::conditional_t<mass_profile, class Empty, Real> M_sw{};
        };
        array_type data;
    };
    constexpr size_t size() const noexcept { return array_size; }
    constexpr iterator begin() noexcept { return data.begin(); }
    constexpr iterator end() noexcept { return data.end(); }
    constexpr const_iterator begin() const noexcept { return data.begin(); }
    constexpr const_iterator end() const noexcept { return data.end(); }
    constexpr Real& operator[](size_t i) noexcept { return data[i]; }
    constexpr const Real& operator[](size_t i) const noexcept { return data[i]; }
};

/********************************************************************************************************************
 * CLASS: SimpleShockEqn
 * DESCRIPTION: Represents the forward shock equation for a given Jet. It defines a state vector
 *              (an array of 8 Reals) and overloads operator() to compute the derivatives of the state with
 *              respect to radius t. It also declares helper functions for the derivatives. Simple version from
 *              Huang et al. 2000
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
class SimpleShockEqn {
   public:
    using State = SimpleState<Ejecta, Medium>;

    SimpleShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, Real eps_e, Real theta_s);

    // Overloaded operator() to compute the derivatives of the state vector with respect to engine time t.
    void operator()(State const& state, State& diff, Real t) const noexcept;

    void setInitState(State& state, Real t0) const noexcept;

    Medium const& medium;  // Reference to the medium properties
    Ejecta const& ejecta;  // Reference to the ejecta properties
    Real const phi{0};     // Angular coordinate phi
    Real const theta0{0};  // Angular coordinate theta
    Real const eps_e{0};   // Electron energy fraction

   private:
    // Helper function: computes the derivative of Gamma with respect to t.
    Real dGammadt(Real dMdt_sw, State const& state, State const& diff) const noexcept;
    Real const dOmega0{0};  // Initial solid angle
    Real const theta_s{0};  // Critical angle for jet spreading
    Real M_ej{0};           // Ejecta mass per solid angle
};

#include "../src/dynamics/simple-shock.tpp"