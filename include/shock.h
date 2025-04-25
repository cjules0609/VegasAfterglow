//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <tuple>

#include "boost/numeric/odeint.hpp"
#include "jet.h"
#include "medium.h"
#include "mesh.h"
#include "physics.h"

/********************************************************************************************************************
 * CLASS: Shock
 * DESCRIPTION: Represents a shock structure that stores grid-based data for the shock evolution, including
 *              comoving time, radius, theta(for jet spreading), bulk Lorentz factor, relative Lorentz factor,
 *              magnetic field, and proton column density.
 *              It also stores constant energy fractions (eps_e and eps_B) and provides a method to return
 *              the grid dimensions.
 ********************************************************************************************************************/
class Shock {
   public:
    Shock(size_t phi_size, size_t theta_size, size_t t_size, Real eps_e, Real eps_B);
    Shock() = default;

    MeshGrid3d t_com;           // comoving time
    MeshGrid3d r;               // radius
    MeshGrid3d theta;           // theta for jet spreading
    MeshGrid3d Gamma;           // bulk lorentz factor
    MeshGrid3d Gamma_rel;       // relative lorentz factor between down stream and up stream
    MeshGrid3d B;               // comoving magnetic field
    MeshGrid3d column_num_den;  // down stream proton column number density
    MeshGrid injection_idx;     // beyond which grid index there is no electron injection
    MaskGrid required;          // Grid that actually required for final flux calculation
    Real eps_e{0};              // electron energy fraction
    Real eps_B{0};              // magnetic energy fraction
    auto shape() const { return std::make_tuple(phi_size, theta_size, t_size); }  // Returns grid dimensions

    void resize(size_t phi_size, size_t theta_size, size_t t_size);

   private:
    size_t phi_size{0};    // Number of grid points in phi direction
    size_t theta_size{0};  // Number of grid points in theta direction
    size_t t_size{0};      // Number of grid points in time direction
};

/********************************************************************************************************************
 * CLASS: FState
 * DESCRIPTION: This is a helper class to provide named access to the forward ODE state array (required by ODE solver)
 *              components. It maps array indices to meaningful physical quantities and provides type-safe access.
 ********************************************************************************************************************/
template <typename StateArray>
struct FState {
    using T = decltype(std::declval<StateArray>()[0]);  // Extract the type of array elements

    // Prevent default construction to ensure proper initialization
    FState() = delete;

    // Constructor that maps array elements to named references
    // The state vector contains: [Gamma, u, r, t_com, theta, M_sw, M_ej, E_ej]
    inline explicit constexpr FState(StateArray& y) noexcept
        : Gamma(y[0]),   // Bulk Lorentz factor of the shock
          u(y[1]),       // Internal energy per unit solid angle
          r(y[2]),       // Radial position of the shock
          t_com(y[3]),   // Comoving time in the shock frame
          theta(y[4]),   // Jet opening angle
          M_sw(y[5]),    // Swept-up mass per unit solid angle
          M_ej(y[6]),    // Ejecta mass per unit solid angle
          E_ej(y[7]) {}  // Ejecta energy per unit solid angle

    // References to state variables - these provide direct access to the underlying array
    T& Gamma;  // Bulk Lorentz factor of the shock
    T& u;      // Internal energy per unit solid angle
    T& r;      // Radial position of the shock
    T& t_com;  // Comoving time in the shock frame
    T& theta;  // Jet opening angle
    T& M_sw;   // Swept-up mass per unit solid angle
    T& M_ej;   // Ejecta mass per unit solid angle
    T& E_ej;   // Ejecta energy per unit solid angle
};

/********************************************************************************************************************
 * CLASS: RState
 * DESCRIPTION: This is a helper class to provide named access to the reverse shock ODE state array (required by ODE
 *              solver) components.
 ********************************************************************************************************************/
template <typename StateArray>
struct RState {
    using T = decltype(std::declval<StateArray>()[0]);
    RState() = delete;

    inline explicit constexpr RState(StateArray& y) noexcept
        : width(y[0]), N3(y[1]), r(y[2]), t_com(y[3]), theta(y[4]), M_sw(y[5]), M_ej(y[6]), E_ej(y[7]) {}
    T& width;
    T& N3;
    T& r;
    T& t_com;
    T& theta;
    T& M_sw;
    T& M_ej;
    T& E_ej;
};

/********************************************************************************************************************
 * INLINE FUNCTIONS: Shock Utilities
 * DESCRIPTION: This section defines a set of inline functions used in shock calculations. These functions compute
 *              various physical quantities such as the comoving magnetic field (via the Weibel instability),
 *              thermal energy density, time derivatives, jet width derivative, downstream number density, fluid
 *              velocities, and update the shock state.
 ********************************************************************************************************************/

Real u_DownStr(Real gamma_rel, Real sigma);

Real u_UpStr2u_DownStr(Real gamma_rel, Real sigma);

inline Real soundSpeed(Real Gamma_rel) {
    Real ad_idx = adiabaticIndex(Gamma_rel);
    return std::sqrt(ad_idx * (ad_idx - 1) * (Gamma_rel - 1) / (1 + (Gamma_rel - 1) * ad_idx));
}

inline Real Gamma_eff(Real adx, Real Gamma) { return (adx * Gamma * Gamma - adx + 1) / Gamma; }

inline Real coMovingWeibelB(Real eps_B, Real e_thermal) { return std::sqrt(8 * con::pi * eps_B * e_thermal); }

inline Real drdt(Real beta) { return (beta * con::c) / (1 - beta); }

inline Real dtheta_dt(Real theta_s, Real theta, Real drdt, Real r, Real Gamma) {
    Real ratio = theta / theta_s;
    // return drdt / (4 * Gamma * theta_s * r * (ratio + 1 / ratio));
    //  return drdt / (4 * Gamma * theta_s * r) / ratio;
    return con::c / (4 * Gamma * theta_s * r * (ratio + 1 / ratio));
    // return con::c / (4 * Gamma * theta_s * r) * ratio;
}

inline Real dtdt_CoMoving(Real Gamma) { return 1 / (Gamma - std::sqrt(Gamma * Gamma - 1)); };

inline Real upStrMagPressure(Real n_up, Real sigma) { return sigma * n_up * con::mp * con::c2 / 2; }

inline Real u_UpStr(Real u_down, Real gamma_rel) {
    return std::sqrt((1 + u_down * u_down) * (gamma_rel * gamma_rel - 1)) + u_down * gamma_rel;
}

inline Real n_DownStr(Real n_up_str, Real gamma_rel, Real sigma) {
    return n_up_str * u_UpStr2u_DownStr(gamma_rel, sigma);
}

inline Real relativeLorentz(Real gamma1, Real gamma2) {
    return gamma1 * gamma2 - std::sqrt((gamma1 * gamma1 - 1) * (gamma2 * gamma2 - 1));
}

inline Real relativeLorentz(Real gamma1, Real gamma2, Real beta1, Real beta2) {
    return gamma1 * gamma2 * (1 - beta1 * beta2);
}

inline Real GammaFromReltive(Real gamma4, Real gamma_rel) {
    Real b = -2 * gamma4 * gamma_rel;
    Real c = gamma4 * gamma4 + gamma_rel * gamma_rel - 1;
    return (-b - std::sqrt(b * b - 4 * c)) / 2;
}

inline Real DownStrThermEnergy(Real gamma_rel, Real n_down_str) {
    return n_down_str * (gamma_rel - 1) * con::mp * con::c2;
}

// D_jet co-moving shell width
inline Real ShellWidthSpredingRate(Real Gamma_rel, Real dtdt_com) {
    Real cs = soundSpeed(Gamma_rel);
    return cs * dtdt_com;
}

inline Real dN3dt(Real r, Real n4, Real gamma3, Real gamma4, Real sigma) {
    if (gamma3 == gamma4) {
        return 0.;
    }
    Real beta3 = gammaTobeta(gamma3);
    Real beta4 = gammaTobeta(gamma4);
    Real gamma34 = relativeLorentz(gamma4, gamma3, beta4, beta3);
    Real ratio_u = u_UpStr2u_DownStr(gamma34, sigma);
    Real n3 = n4 * ratio_u;
    Real dxdt = (beta4 - beta3) * con::c / ((1 - beta3) * (1 - gamma4 / (gamma3 * ratio_u)));
    return n3 * r * r * gamma3 * dxdt;
}

inline Real calc_n4(Real dEdOmega, Real Gamma0, Real r, Real D_jet, Real sigma) {
    return dEdOmega / ((Gamma0 * con::mp * con::c2 * r * r * D_jet) * (1 + sigma));
}

template <typename State>
inline void setStoppingShock(size_t i, size_t j, Shock& shock, State const& state0) {
    xt::view(shock.t_com, i, j, xt::all()) = state0.t_com;
    xt::view(shock.r, i, j, xt::all()) = state0.r;
    xt::view(shock.theta, i, j, xt::all()) = state0.theta;
    xt::view(shock.Gamma, i, j, xt::all()) = 1;
    xt::view(shock.Gamma_rel, i, j, xt::all()) = 1;
    xt::view(shock.B, i, j, xt::all()) = 0;
    xt::view(shock.column_num_den, i, j, xt::all()) = 0;
}

template <typename State>
void updateShockState(Shock& shock, size_t i, size_t j, size_t k, State const& state, Real Gamma_down_str,
                      Real Gamma_up_str, Real N_down_str, Real n_up_str, Real sigma_up_str) {
    Real Gamma_rel = relativeLorentz(Gamma_up_str, Gamma_down_str);
    Real ratio_u = u_UpStr2u_DownStr(Gamma_rel, sigma_up_str);
    Real pB_up_str = upStrMagPressure(n_up_str, sigma_up_str);
    Real pB_down_str = pB_up_str * ratio_u * ratio_u;
    Real n_down_str = n_up_str * ratio_u;
    Real e_th = DownStrThermEnergy(Gamma_rel, n_down_str);
    shock.t_com(i, j, k) = state.t_com;
    shock.r(i, j, k) = state.r;
    shock.theta(i, j, k) = state.theta;
    shock.Gamma(i, j, k) = Gamma_down_str;
    shock.Gamma_rel(i, j, k) = Gamma_rel;
    shock.B(i, j, k) = coMovingWeibelB(shock.eps_B, e_th) + std::sqrt(pB_down_str * 8 * con::pi);
    shock.column_num_den(i, j, k) = N_down_str / (state.r * state.r);
}
