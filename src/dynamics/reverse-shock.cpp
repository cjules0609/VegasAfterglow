//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "reverse-shock.hpp"

/********************************************************************************************************************
 * CONSTRUCTOR: FRShockEqn::FRShockEqn
 * DESCRIPTION: constructor for the FRShockEqn class.
 ********************************************************************************************************************/
FRShockEqn::FRShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, Real eps_e)
    : medium(medium), ejecta(ejecta), phi(phi), theta0(theta), eps_e(eps_e) {}

/********************************************************************************************************************
 * METHOD: FRShockEqn::operator()(State const& y, State& dydt, Real t)
 * DESCRIPTION: Reverse shock ODE.
 ********************************************************************************************************************/
void FRShockEqn::operator()(StateArray const& y, StateArray& dydt, Real t) {
    RState state(y);
    RState diff(dydt);

    Real gamma3 = 1;
    Real beta3 = 1;
    Real gamma_rel = 1;
    if (!crossed) {
        Real n4 = state.M_ej / (state.r * state.r * state.width * con::mp);
        Real N2 = state.M_sw / con::mp;
        Real gamma4 = ejecta.Gamma0(phi, theta0);
        Real sigma4 = state.E_ej / (gamma4 * state.M_ej * con::c2) - 1;
        gamma3 = calc_gamma3(gamma4, N2, state.N3, sigma4);
        beta3 = gammaTobeta(gamma3);
        gamma_rel = relativeLorentz(gamma4, gamma3);
        diff.N3 = dN3dt(state.r, n4, gamma3, gamma4, sigma4);
    } else {
        gamma_rel = crossedGamma_rel(state);
        gamma3 = crossedGamma3(gamma_rel, state.r);
        beta3 = gammaTobeta(gamma3);
        diff.N3 = 0;
    }
    diff.r = drdt(beta3);
    diff.t_com = dtdt_CoMoving(gamma3);
    diff.width = dDdt_Jet(gamma_rel, diff.t_com);
    Real rho = medium.rho(phi, state.theta, state.r);
    diff.M_sw = state.r * state.r * rho * diff.r;
    diff.M_ej = ejecta.dMdtdOmega(phi, theta0, t);
    diff.E_ej = ejecta.dEdtdOmega(phi, theta0, t);
    diff.theta = 0;  // no spreading
}

/********************************************************************************************************************
 * FUNCTION: gParameter
 * DESCRIPTION: post shock crossing power law index for the Gamma3*beta3 \propto r^-g.
 ********************************************************************************************************************/
inline Real gParameter(Real gamma_rel, Real k = 0) {
    constexpr Real g_low = 1.5;   // k is the medium power law index
    constexpr Real g_high = 3.5;  // Blandford-McKee limit// TODO: need to be modified for non ISM medium
    Real p = std::sqrt(std::sqrt(gamma_rel - 1));
    return g_low + (g_high - g_low) * p / (1 + p);
}

/********************************************************************************************************************
 * METHOD: FRShockEqn::crossedGamma3(Real gamma_rel, Real r)
 * DESCRIPTION: Post crossing gamma3*beta3 profile that is \propto r^-g. Using gamma3*beta3 for Newtonian regime as
 *              well.
 ********************************************************************************************************************/
Real FRShockEqn::crossedGamma3(Real gamma_rel, Real r) const {
    Real g = gParameter(gamma_rel);
    Real u = u_x * std::pow(r / r_x, -g);
    return std::sqrt(u * u + 1);
}

/********************************************************************************************************************
 * FUNCTION: calc_gamma3
 * DESCRIPTION: Calculate the reverse shock crossing gamma3
 ********************************************************************************************************************/
Real calc_gamma3(Real gamma4, Real N2, Real N3, Real sigma) {
    auto func = [=](Real gamma3) -> Real {
        Real gamma34 = relativeLorentz(gamma4, gamma3);
        Real adx3 = adiabaticIndex(gamma34);
        Real adx2 = adiabaticIndex(gamma3);
        Real g_eff2 = Gamma_eff(adx2, gamma3);
        Real g_eff3 = Gamma_eff(adx3, gamma3);

        Real E2 = N2 * (gamma3 - 1 + g_eff2 * (gamma3 - 1));
        Real E3 = N3 * (gamma3 - gamma4 + g_eff3 * (gamma34 - 1)) * (1 + sigma);

        return E2 + E3;
    };
    constexpr Real r_tol = 1e-3;
    return rootBisection(func, 1, gamma4, r_tol);
}