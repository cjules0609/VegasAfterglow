//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _SIMPLESHOCK_
#define _SIMPLESHOCK_
#include "shock.h"
/********************************************************************************************************************
 * CLASS: SimpleShockEqn
 * DESCRIPTION: Represents the forward shock equation for a given Jet. It defines a state vector
 *              (an array of 8 Reals) and overloads operator() to compute the derivatives of the state with
 *              respect to radius t. It also declares helper functions for the derivatives. Simple version from
 *              Huang et al. 2000
 ********************************************************************************************************************/
class SimpleShockEqn {
   public:
    using StateArray =
        std::array<Real, 8>;  // State vector: typically [Gamma, u, r, t_com, theta_jet, m, M_ej, E_ejecta]

    SimpleShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta_lo, Real theta, Real eps_e,
                   Real theta_s);

    Medium const& medium;  // Reference to the medium properties
    Ejecta const& ejecta;  // Reference to the ejecta properties
    Real const phi{0};     // Angular coordinate phi
    Real const theta0{0};  // Angular coordinate theta
    Real const eps_e{0};   // Electron energy fraction

    // Overloaded operator() to compute the derivatives of the state vector with respect to radius r.
    void operator()(StateArray const& y, StateArray& dydt, Real t);

   private:
    // Helper function: computes the derivative of Gamma with respect to t.
    inline Real dGammadt(Real t, FState<const StateArray> const& state, FState<StateArray> const& diff);
    Real const dOmega0{0};  // Initial solid angle
    Real const theta_s{0};  // Critical angle for jet spreading
    Real const theta_lo{0};
};


#endif