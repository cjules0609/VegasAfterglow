//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _SHOCKDYNAMICS_
#define _SHOCKDYNAMICS_

#include <boost/numeric/odeint.hpp>
#include <tuple>

#include "jet.h"
#include "medium.h"
#include "mesh.h"
#include "physics.h"

/********************************************************************************************************************
 * CLASS: Shock
 * DESCRIPTION: Represents a shock structure that stores grid-based data for the shock evolution, including
 *              comoving time, engine time, relative Lorentz factor, magnetic field, and proton column density.
 *              It also stores constant energy fractions (eps_e and eps_B) and provides a method to return
 *              the grid dimensions.
 ********************************************************************************************************************/
class Shock {
   public:
    Shock(size_t phi_size, size_t theta_size, size_t r_size, Real eps_e, Real eps_B);
    Shock();

    MeshGrid3d t_com;           // comoving time
    MeshGrid3d t_eng;           // engine time
    MeshGrid3d Gamma_rel;       // relative lorentz factor between down stream and up stream
    MeshGrid3d B;               // comoving magnetic field
    MeshGrid3d column_num_den;  // down stream proton column number density
    Real eps_e{0};              // electron energy fraction
    Real eps_B{0};              // magnetic energy fraction

    auto shape() const { return std::make_tuple(phi_size, theta_size, r_size); }  // Returns grid dimensions

   private:
    size_t phi_size{0};    // Number of grid points in phi direction
    size_t theta_size{0};  // Number of grid points in theta direction
    size_t r_size{0};      // Number of grid points in radial direction
};

/********************************************************************************************************************
 * FUNCTION PROTOTYPES: Shock Generation Interfaces
 * DESCRIPTION: These function templates declare interfaces to generate forward shocks (2D and 3D) and
 *              forward/reverse shock pairs.
 ********************************************************************************************************************/
using ShockPair = std::pair<Shock, Shock>;

template <typename Jet, typename Injector>
Shock genForwardShock(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, Real eps_e,
                      Real eps_B);

template <typename Jet, typename Injector>
Shock genForwardShock3D(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, Real eps_e,
                        Real eps_B);

template <typename Jet, typename Injector>
ShockPair genFRShocks(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, Real eps_e,
                      Real eps_B);

template <typename Jet, typename Injector>
ShockPair genFRShocks3D(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, Real eps_e,
                        Real eps_B);

/********************************************************************************************************************
 * INLINE FUNCTIONS: Shock Utilities
 * DESCRIPTION: This section defines a set of inline functions used in shock calculations. These functions compute
 *              various physical quantities such as the comoving magnetic field (via the Weibel instability),
 *              thermal energy density, time derivatives, jet width derivative, downstream number density, fluid
 *              velocities, and update the shock state.
 ********************************************************************************************************************/

Real u_DownStr(Real gamma_rel, Real sigma);
Real u_UpStr2u_DownStr(Real gamma_rel, Real sigma);
void updateShockState(Shock& shock, size_t i, size_t j, size_t k, Real r, Real Gamma_rel, Real t_com, Real t_eng,
                      Real dMdOmega_up, Real n_up_str, Real sigma);

inline Real coMovingWeibelB(Real eps_B, Real e_thermal) { return sqrt(8 * con::pi * eps_B * e_thermal); }
inline Real dtdr_Engine(Real beta) { return std::abs(1 - beta) / (beta * con::c); }
inline Real dtdr_CoMoving(Real Gamma, Real beta) { return 1 / (Gamma * beta * con::c); };
inline Real calc_pB4(Real n4, Real sigma) { return sigma * n4 * con::mp * con::c2 / 2; }
inline Real u_UpStr(Real u_down, Real gamma_rel) {
    return std::sqrt((1 + u_down * u_down) * (gamma_rel * gamma_rel - 1)) + u_down * gamma_rel;
}
inline Real n_DownStr(Real n_up_str, Real gamma_rel, Real sigma) {
    return n_up_str * u_UpStr2u_DownStr(gamma_rel, sigma);
}
inline Real e_ThermalDownStr(Real gamma_rel, Real n_down_str) {
    return n_down_str * (gamma_rel - 1) * con::mp * con::c2;
}
inline Real dDdr_Jet(Real Gamma, Real beta) {
    Real constexpr cs = 0.5773502691896258 * con::c;  // sound speed approximation factor
    return cs * dtdr_CoMoving(Gamma, beta) / Gamma;
}
inline Real calc_n4(Real dEdOmega, Real Gamma0, Real r, Real D_jet_lab, Real sigma) {
    return dEdOmega / (Gamma0 * con::mp * con::c2 * r * r * Gamma0 * D_jet_lab) / (1 + sigma);
}

//

template <typename Jet, typename Injector>
class SimpleShockEqn;
template <typename Jet, typename Injector>
class ForwardShockEqn;
template <typename Jet, typename Injector>
class FRShockEqn;

/********************************************************************************************************************
 * FUNCTION: genForwardShock
 * DESCRIPTION: Generates a forward shock (2D) using the provided coordinates, medium, jet, injector, and energy
 *fractions. It creates a Shock object for a single phi value and iterates over theta values, solving the shock
 *              evolution for each theta slice.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
Shock genForwardShock(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, Real eps_e,
                      Real eps_B) {
    auto [phi_size, theta_size, r_size] = coord.shape();  // Unpack coordinate dimensions
    Shock f_shock(1, theta_size, r_size, eps_e, eps_B);   // Create Shock with 1 phi slice
    for (size_t j = 0; j < theta_size; ++j) {
        // Create a ForwardShockEqn for each theta slice (phi is set to 0)
        auto eqn = ForwardShockEqn(medium, jet, inject, 0, coord.theta[j], eps_e);
        // auto eqn = SimpleShockEqn(medium, jet, inject, 0, coord.theta[j], eps_e);
        //     Solve the shock shell for this theta slice
        solveForwardShell(0, j, coord.r, f_shock, eqn, coord.t_max);
    }

    return f_shock;
}

/********************************************************************************************************************
 * FUNCTION: genForwardShock3D
 * DESCRIPTION: Generates a forward shock (3D) using the provided coordinates, medium, jet, injector, and energy
 *fractions. It creates a Shock object covering all phi and theta slices and iterates over both dimensions.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
Shock genForwardShock3D(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, Real eps_e,
                        Real eps_B) {
    auto [phi_size, theta_size, r_size] = coord.shape();

    Shock f_shock(phi_size, theta_size, r_size, eps_e, eps_B);  // Create Shock with full 3D dimensions
    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            // Create a ForwardShockEqn for each (phi, theta) pair
            auto eqn = ForwardShockEqn(medium, jet, inject, coord.phi[i], coord.theta[j], eps_e);
            // Solve the shock shell for this (phi, theta) slice
            solveForwardShell(i, j, coord.r, f_shock, eqn, coord.t_max);
        }
    }
    return f_shock;
}

/********************************************************************************************************************
 * FUNCTION: genFRShocks
 * DESCRIPTION: Generates a pair of forward and reverse shocks (2D) using the provided coordinates, medium, jet,
 *              injector, and energy fractions. It creates two Shock objects (one forward, one reverse) and solves
 *              the shock shells for each theta slice.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
ShockPair genFRShocks(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, Real eps_e,
                      Real eps_B) {
    auto [phi_size, theta_size, r_size] = coord.shape();

    Shock f_shock(1, theta_size, r_size, eps_e, eps_B);  // Forward shock for 1 phi slice
    Shock r_shock(1, theta_size, r_size, eps_e, eps_B);  // Reverse shock for 1 phi slice

    for (size_t j = 0; j < theta_size; ++j) {
        // Create equations for forward and reverse shocks for each theta slice (phi is 0)
        auto eqn_f = ForwardShockEqn(medium, jet, inject, 0, coord.theta[j], eps_e);
        auto eqn_r = FRShockEqn(medium, jet, inject, 0, coord.theta[j]);
        // Solve the forward-reverse shock shell
        solveFRShell(0, j, coord.r, f_shock, r_shock, eqn_f, eqn_r, coord.t_max);
    }

    return std::make_pair(std::move(f_shock), std::move(r_shock));
}

/********************************************************************************************************************
 * FUNCTION: genFRShocks3D
 * DESCRIPTION: Generates a pair of forward and reverse shocks (3D) using the provided coordinates, medium, jet,
 *              injector, and energy fractions. It creates two Shock objects covering all phi and theta slices and
 *              solves the shock shells for each (phi, theta) pair.
 ********************************************************************************************************************/
template <typename Jet, typename Injector>
ShockPair genFRShocks3D(Coord const& coord, Medium const& medium, Jet const& jet, Injector const& inject, Real eps_e,
                        Real eps_B) {
    auto [phi_size, theta_size, r_size] = coord.shape();

    Shock f_shock(phi_size, theta_size, r_size, eps_e, eps_B);  // Forward shock for full 3D dimensions
    Shock r_shock(phi_size, theta_size, r_size, eps_e, eps_B);  // Reverse shock for full 3D dimensions
    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            // Create equations for forward and reverse shocks for each (phi, theta) pair
            auto eqn_f = ForwardShockEqn(medium, jet, inject, coord.phi[i], coord.theta[j], eps_e);
            auto eqn_r = FRShockEqn(medium, jet, inject, coord.phi[i], coord.theta[j]);
            // Solve the forward-reverse shock shell for this (phi, theta) pair
            solveFRShell(i, j, coord.r, f_shock, r_shock, eqn_f, eqn_r, coord.t_max);
        }
    }
    return std::make_pair(std::move(f_shock), std::move(r_shock));
}

#include "forward-shock.hpp"
#include "reverse-shock.hpp"
#include "simple-shock.hpp"
#endif