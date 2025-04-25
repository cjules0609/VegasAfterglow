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
 * CLASS: ForwardShockEqn
 * DESCRIPTION: Represents the forward shock equation for a given jet and medium. It defines a state vector
 *              (an array of 8 Reals) and overloads operator() to compute the derivatives of the state with
 *              respect to t. It also declares helper functions for the derivatives.
 *              This class implements the physical equations governing the forward shock evolution.
 ********************************************************************************************************************/
template <typename Ejecta, typename Medium>
class ForwardShockEqn {
   public:
    // State vector: [Gamma, u, r, t_com, theta, M_sw, M_ej, E_ej]
    using StateArray = std::array<Real, 8>;       // 8-dimensional state vector
    using State = FState<StateArray>;             // Mutable state wrapper
    using constState = FState<const StateArray>;  // Const state wrapper

    // Constructor: Initialize the forward shock equation with physical parameters
    ForwardShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, Real eps_e, Real theta_s);

    // References to model components
    Medium const& medium;  // Reference to the ambient medium properties
    Ejecta const& ejecta;  // Reference to the ejecta properties

    // Model parameters
    Real const phi{0};     // Angular coordinate phi in the jet frame
    Real const theta0{0};  // Initial angular coordinate theta
    Real const eps_e{0};   // Fraction of energy given to electrons

    // Forward shock ODE equation - callable interface for ODE solver
    // Computes the time derivatives of the state vector
    void operator()(StateArray const& y, StateArray& dydt, Real t) const noexcept;

   private:
    // Helper methods for computing derivatives

    // Computes the time derivative of the Lorentz factor
    // with respect to on-axis observer time
    inline Real dGammadt(Real t, constState const& state, State const& diff, Real ad_idx) const noexcept;

    // Computes the time derivative of internal energy
    // with respect to on-axis observer time
    inline Real dUdt(constState const& state, State const& diff, Real ad_idx) const noexcept;

    // Private member variables
    Real const dOmega0{0};  // Initial solid angle element
    Real const theta_s{0};  // Jet structure parameter controlling angular dependence
};

#include "../src/dynamics/forward-shock.tpp"