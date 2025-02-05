//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "synchrotron.h"

#include <boost/math/tools/roots.hpp>
#include <cmath>

#include "afterglow.h"
#include "inverse-compton.h"
#include "macros.h"
#include "physics.h"
#include "utilities.h"
/********************************************************************************************************************
 * FUNCTION: InverseComptonY::InverseComptonY(Real nu_m, Real nu_c, Real B, Real Y_T)
 * DESCRIPTION: Constructor that initializes an InverseComptonY object using the provided minimum and cooling
 *              frequencies (nu_m, nu_c), magnetic field B, and Y_T parameter. It computes characteristic
 *              gamma values and corresponding synchrotron frequencies, and then sets the cooling regime.
 ********************************************************************************************************************/
InverseComptonY::InverseComptonY(Real nu_m, Real nu_c, Real B, Real Y_T) {
    gamma_hat_m = con::me * con::c2 / (con::h * nu_m);  // Compute minimum characteristic Lorentz factor
    gamma_hat_c = con::me * con::c2 / (con::h * nu_c);  // Compute cooling characteristic Lorentz factor
    this->Y_T = Y_T;                                    // Set the Thomson Y parameter
    nu_hat_m = syn_nu(gamma_hat_m, B);                  // Compute corresponding synchrotron frequency for gamma_hat_m
    nu_hat_c = syn_nu(gamma_hat_c, B);                  // Compute corresponding synchrotron frequency for gamma_hat_c

    if (nu_hat_m <= nu_hat_c) {
        regime = 1;  // fast IC cooling regime
    } else {
        regime = 2;  // slow IC cooling regime
    }
}

/********************************************************************************************************************
 * FUNCTION: InverseComptonY::InverseComptonY(Real Y_T)
 * DESCRIPTION: Constructor that initializes an InverseComptonY object with only Y_T provided.
 ********************************************************************************************************************/
InverseComptonY::InverseComptonY(Real Y_T) {
    this->Y_T = Y_T;  // Set the Thomson Y parameter
    regime = 3;       // Set regime to 3 (special case)
}

/********************************************************************************************************************
 * FUNCTION: InverseComptonY::InverseComptonY()
 * DESCRIPTION: Default constructor initializing all member variables to zero.
 ********************************************************************************************************************/
InverseComptonY::InverseComptonY() {
    nu_hat_m = 0;
    nu_hat_c = 0;
    gamma_hat_m = 0;
    gamma_hat_c = 0;
    Y_T = 0;
    regime = 0;
}

/********************************************************************************************************************
 * FUNCTION: InverseComptonY::as_gamma(Real gamma, Real p) const
 * DESCRIPTION: Returns the effective Y parameter as a function of electron Lorentz factor (gamma) and
 *              power-law index (p), based on the cooling regime.
 ********************************************************************************************************************/
Real InverseComptonY::as_gamma(Real gamma, Real p) const {
    switch (regime) {
        case 3:
            return Y_T;  // In regime 3, simply return Y_T
            break;
        case 1:
            if (gamma <= gamma_hat_m) {
                return Y_T;  // For gamma below gamma_hat_m, no modification
            } else if (gamma <= gamma_hat_c) {
                return Y_T / std::sqrt(gamma / gamma_hat_m);  // Intermediate regime scaling
            } else
                return Y_T * pow43(gamma_hat_c / gamma) / std::sqrt(gamma_hat_c / gamma_hat_m);  // High gamma scaling

            break;
        case 2:
            if (gamma <= gamma_hat_c) {
                return Y_T;  // For gamma below gamma_hat_c, no modification
            } else if (gamma <= gamma_hat_m) {
                return Y_T * fastPow(gamma / gamma_hat_c, (p - 3) / 2);  // Scaling in intermediate regime
            } else
                return Y_T * pow43(gamma_hat_m / gamma) *
                       fastPow(gamma_hat_m / gamma_hat_c, (p - 3) / 2);  // High gamma scaling

            break;
        default:
            return 0;
            break;
    }
}

/********************************************************************************************************************
 * FUNCTION: InverseComptonY::as_nu(Real nu, Real p) const
 * DESCRIPTION: Returns the effective Y parameter as a function of frequency (nu) and power-law index (p),
 *              based on the cooling regime.
 ********************************************************************************************************************/
Real InverseComptonY::as_nu(Real nu, Real p) const {
    switch (regime) {
        case 3:
            return Y_T;  // In regime 3, simply return Y_T
            break;
        case 1:
            if (nu <= nu_hat_m) {
                return Y_T;  // For frequencies below nu_hat_m, no modification
            } else if (nu <= nu_hat_c) {
                return Y_T * std::sqrt(std::sqrt(nu_hat_m / nu));  // Intermediate frequency scaling
            } else
                return Y_T * pow23(nu_hat_c / nu) * std::sqrt(std::sqrt(nu_hat_m / nu));  // High frequency scaling

            break;
        case 2:
            if (nu <= nu_hat_c) {
                return Y_T;  // For frequencies below nu_hat_c, no modification
            } else if (nu <= nu_hat_m) {
                return Y_T * fastPow(nu / nu_hat_c, (p - 3) / 4);  // Intermediate frequency scaling
            } else
                return Y_T * pow23(nu_hat_m / nu) *
                       fastPow(nu_hat_m / nu_hat_c, (p - 3) / 4);  // High frequency scaling

            break;
        default:
            return 0;
            break;
    }
}

/********************************************************************************************************************
 * FUNCTION: InverseComptonY::Y_Thompson(std::vector<InverseComptonY> const& Ys)
 * DESCRIPTION: Computes the total Y parameter by summing Y_T of all InverseComptonY objects in the vector.
 ********************************************************************************************************************/
Real InverseComptonY::Y_Thompson(InverseComptonY const& Ys) {
    /*Real Y_tilt = 0;
    for (auto& Y : Ys) {
        Y_tilt += Y.Y_T;  // Sum each object's Y_T
    }
    return Y_tilt;*/
    return Ys.Y_T;
}

/********************************************************************************************************************
 * FUNCTION: InverseComptonY::Y_tilt_gamma(std::vector<InverseComptonY> const& Ys, Real gamma, Real p)
 * DESCRIPTION: Computes the total effective Y parameter as a function of electron Lorentz factor (gamma) and
 *              power-law index (p) by summing contributions from all InverseComptonY objects.
 ********************************************************************************************************************/
Real InverseComptonY::Y_tilt_gamma(InverseComptonY const& Ys, Real gamma, Real p) {
    /*Real Y_tilt = 0;
    for (auto& Y : Ys) {
        Y_tilt += Y.as_gamma(gamma, p);  // Sum effective Y parameters based on gamma
    }
    return Y_tilt;*/
    return Ys.as_gamma(gamma, p);
}

/********************************************************************************************************************
 * FUNCTION: InverseComptonY::Y_tilt_nu(std::vector<InverseComptonY> const& Ys, Real nu, Real p)
 * DESCRIPTION: Computes the total effective Y parameter as a function of frequency (nu) and power-law index (p)
 *              by summing contributions from all InverseComptonY objects.
 ********************************************************************************************************************/
Real InverseComptonY::Y_tilt_nu(InverseComptonY const& Ys, Real nu, Real p) {
    /* Real Y_tilt = 0;
     for (auto& Y : Ys) {
         Y_tilt += Y.as_nu(nu, p);  // Sum effective Y parameters based on frequency
     }
     return Y_tilt;*/
    return Ys.as_nu(nu, p);
}

/********************************************************************************************************************
 * FUNCTION: createSynPhotonGrid(size_t phi_size, size_t theta_size, size_t t_size)
 * DESCRIPTION: Creates and returns a SynPhotonGrid with the specified dimensions.
 ********************************************************************************************************************/
SynPhotonGrid createSynPhotonGrid(size_t phi_size, size_t theta_size, size_t t_size) {
    SynPhotonGrid grid(boost::extents[phi_size][theta_size][t_size]);
    return grid;
}

/********************************************************************************************************************
 * FUNCTION: createSynElectronGrid(size_t phi_size, size_t theta_size, size_t t_size)
 * DESCRIPTION: Creates and returns a SynElectronGrid with the specified dimensions.
 ********************************************************************************************************************/
SynElectronGrid createSynElectronGrid(size_t phi_size, size_t theta_size, size_t t_size) {
    SynElectronGrid grid(boost::extents[phi_size][theta_size][t_size]);
    return grid;
}

/********************************************************************************************************************
 * FUNCTION: SynElectrons::columnNumDen(Real gamma) const
 * DESCRIPTION: Computes the column number density for electrons at a given Lorentz factor (gamma) by scaling
 *              the base column density with the electron spectrum. It accounts for inverse Compton effects.
 ********************************************************************************************************************/
Real SynElectrons::columnNumDen(Real gamma) const {
    if (gamma < gamma_c) {
        return column_num_den * gammaSpectrum(gamma);  // Below cooling Lorentz factor: direct scaling
    } else {
        return column_num_den * gammaSpectrum(gamma) * (1 + Y_c) / (1 + InverseComptonY::Y_tilt_gamma(Ys, gamma, p));
        // Above cooling Lorentz factor: include inverse Compton Y tilt correction
    }
}

/********************************************************************************************************************
 * FUNCTION: order(Real a, Real b, Real c)
 * DESCRIPTION: Helper inline function to check if three numbers are in non-decreasing order.
 ********************************************************************************************************************/
inline bool order(Real a, Real b, Real c) { return a <= b && b <= c; };

/********************************************************************************************************************
 * FUNCTION: getRegime(Real a, Real c, Real m)
 * DESCRIPTION: Determines the regime based on the ordering of three parameters (typically gamma_a, gamma_c, gamma_m).
 *              Returns a regime identifier (1-6) or 0 if none match.
 ********************************************************************************************************************/
size_t getRegime(Real a, Real c, Real m) {
    if (order(a, m, c)) {
        return 1;
    } else if (order(m, a, c)) {
        return 2;
    } else if (order(a, c, m)) {
        return 3;
    } else if (order(c, a, m)) {
        return 4;
    } else if (order(m, c, a)) {
        return 5;
    } else if (order(c, m, a)) {
        return 6;
    } else
        return 0;
}

/********************************************************************************************************************
 * FUNCTION: SynElectrons::gammaSpectrum(Real gamma) const
 * DESCRIPTION: Returns the electron energy spectrum as a function of Lorentz factor (gamma) based on the
 *              current regime. Different scaling laws apply in different regimes.
 ********************************************************************************************************************/
Real SynElectrons::gammaSpectrum(Real gamma) const {
    switch (regime) {
        case 1:  // same as case 2
        case 2:
            if (gamma <= gamma_m) {
                return 0;  // Below minimum Lorentz factor, spectrum is zero
            } else if (gamma <= gamma_c) {
                return (p - 1) * fastPow(gamma / gamma_m, -p) /
                       gamma_m;  // Power-law spectrum between gamma_m and gamma_c
            } else
                return (p - 1) * fastPow(gamma / gamma_m, -p) * gamma_c / (gamma * gamma_m) * fastExp(-gamma / gamma_M);
            // Above cooling Lorentz factor: exponential cutoff applied

            break;
        case 3:
            if (gamma <= gamma_c) {
                return 0;  // Below cooling Lorentz factor, spectrum is zero
            } else if (gamma <= gamma_m) {
                return gamma_c / (gamma * gamma);  // Intermediate regime scaling
            } else
                return gamma_c / (gamma * gamma_m) * fastPow(gamma / gamma_m, -p) * fastExp(-gamma / gamma_M);
            // Above minimum Lorentz factor: power-law with exponential cutoff

            break;
        case 4:  // Gao, Lei, Wu and Zhang 2013 Eq 18
            if (gamma <= gamma_a) {
                return 3 * gamma * gamma / (gamma_a * gamma_a * gamma_a);  // Rising part of the spectrum
            } else if (gamma <= gamma_m) {
                return gamma_c / (gamma * gamma);  // Transition region
            } else
                return gamma_c / (gamma * gamma_m) * fastPow(gamma / gamma_m, -p) * fastExp(-gamma / gamma_M);
            // High energy tail with exponential cutoff

            break;
        case 5:  // Gao, Lei, Wu and Zhang 2013 Eq 19
            if (gamma <= gamma_a) {
                return 3 * gamma * gamma / (gamma_a * gamma_a * gamma_a);  // Rising part of the spectrum
            } else
                return (p - 1) * gamma_c / (gamma * gamma_m) * fastPow(gamma / gamma_m, -p) * fastExp(-gamma / gamma_M);
            // Power-law decay with exponential cutoff

            break;
        case 6:  // Gao, Lei, Wu and Zhang 2013 Eq 20
            if (gamma <= gamma_a) {
                return 3 * gamma * gamma / (gamma_a * gamma_a * gamma_a);  // Rising part of the spectrum
            } else
                return fastPow(gamma_m, p - 1) * gamma_c * fastPow(gamma, -(p + 1)) * fastExp(-gamma / gamma_M);
            // Steeper decay in this regime

            break;
        default:
            return 0;
    }
}

/********************************************************************************************************************
 * FUNCTION: SynPhotons::I_nu(Real nu) const
 * DESCRIPTION: Computes the synchrotron photon intensity at frequency nu based on the electron peak intensity
 *              and the computed spectrum. Adjusts for inverse Compton effects if nu exceeds nu_c.
 ********************************************************************************************************************/
Real SynPhotons::I_nu(Real nu) const {
    if (nu < nu_c) {
        return e->I_nu_peak * spectrum(nu);  // Below cooling frequency, simple scaling
    } else {
        return e->I_nu_peak * spectrum(nu) * (1 + e->Y_c) / (1 + InverseComptonY::Y_tilt_nu(e->Ys, nu, e->p));
        // Above cooling frequency, include inverse Compton correction
    }
}

/********************************************************************************************************************
 * FUNCTION: SynPhotons::updateConstant()
 * DESCRIPTION: Updates internal constants used for computing the synchrotron photon spectrum based on the
 *              current spectral parameters from the associated electrons.
 ********************************************************************************************************************/
void SynPhotons::updateConstant() {
    // Update constants based on spectral parameters
    Real p = e->p;
    if (e->regime == 1) {
        // a_m_1_3 = std::cbrt(nu_a / nu_m);  // (nu_a / nu_m)^(1/3)
        // c_m_1_2 = std::sqrt(nu_c / nu_m);  // (nu_c / nu_m)^(1/2)
        C1_ = std::cbrt(nu_a / nu_m);
        C2_ = std::sqrt(nu_c / nu_m);
    } else if (e->regime == 2) {
        // m_a_pa4_2 = fastPow(nu_m / nu_a, (p + 4) / 2);    // (nu_m / nu_a)^((p+4)/2)
        // a_m_mpa1_2 = fastPow(nu_a / nu_m, (-p + 1) / 2);  // (nu_a / nu_m)^((-p+1)/2)
        // c_m_1_2 = std::sqrt(nu_c / nu_m);
        C1_ = fastPow(nu_m / nu_a, (p + 4) / 2);
        C2_ = fastPow(nu_a / nu_m, (-p + 1) / 2);
        C3_ = std::sqrt(nu_c / nu_m);
    } else if (e->regime == 3) {
        // a_c_1_3 = std::cbrt(nu_a / nu_c);  // (nu_a / nu_c)^(1/3)
        // c_m_1_2 = std::sqrt(nu_c / nu_m);  // (nu_c / nu_m)^(1/2)
        C1_ = std::cbrt(nu_a / nu_c);
        C2_ = std::sqrt(nu_c / nu_m);
    } else if (e->regime == 4) {
        // a_m_1_2 = std::sqrt(nu_a / nu_m);  // (nu_a / nu_m)^(1/2)
        // R4 = std::sqrt(nu_c / nu_a) / 3;   // (nu_c / nu_a)^(1/2) / 3; // R4: scaling factor for regime 4
        C1_ = std::sqrt(nu_a / nu_m);
        C2_ = std::sqrt(nu_c / nu_a) / 3;
    } else if (e->regime == 5 || e->regime == 6) {
        // R4 = std::sqrt(nu_c / nu_a) / 3;              // (nu_c / nu_a)^(1/2) / 3; // R4: scaling factor for regime 4
        // R6 = R4 * fastPow(nu_m / nu_a, (p - 1) / 2);  // R6: scaling factor for regime 6
        C1_ = std::sqrt(nu_c / nu_a) / 3;
        C2_ = C1_ * fastPow(nu_m / nu_a, (p - 1) / 2);
    }
    // R5 = (p - 1) * R6;  // R5 scales R6 (commented out)
}

/********************************************************************************************************************
 * FUNCTION: SynPhotons::spectrum(Real nu) const
 * DESCRIPTION: Computes the synchrotron photon spectrum at frequency nu based on the regime of the associated
 *              electrons. Different formulae apply in each regime.
 ********************************************************************************************************************/
Real SynPhotons::spectrum(Real nu) const {
    Real p = e->p;
    switch (e->regime) {
        case 1:
            if (nu <= nu_a) {
                return C1_ * (nu / nu_a) * (nu / nu_a);
            }
            if (nu <= nu_m) {
                return std::cbrt(nu / nu_m);
            }
            if (nu <= nu_c) {
                return fastPow(nu / nu_m, -(p - 1) / 2);
            }
            return C2_ * fastPow(nu / nu_m, -p / 2) * fastExp(-nu / nu_M);

            break;
        case 2:
            if (nu <= nu_m) {
                return C1_ * (nu / nu_m) * (nu / nu_m);
            }
            if (nu <= nu_a) {
                return C2_ * pow52(nu / nu_a);  // Using pow52 for (nu / nu_a)^(5/2)
            }
            if (nu <= nu_c) {
                return fastPow(nu / nu_m, -(p - 1) / 2);
            }
            return C3_ * fastPow(nu / nu_m, -p / 2) * fastExp(-nu / nu_M);

            break;
        case 3:
            if (nu <= nu_a) {
                return C1_ * (nu / nu_a) * (nu / nu_a);
            }
            if (nu <= nu_c) {
                return std::cbrt(nu / nu_c);
            }
            if (nu <= nu_m) {
                return std::sqrt(nu_c / nu);
            }
            return C2_ * fastPow(nu / nu_m, -p / 2) * fastExp(-nu / nu_M);

            break;
        case 4:
            if (nu <= nu_a) {
                return (nu / nu_a) * (nu / nu_a);
            }
            if (nu <= nu_m) {
                return C2_ * std::sqrt(nu_a / nu);
            }
            return C2_ * C1_ * fastPow(nu / nu_m, -p / 2) * fastExp(-nu / nu_M);

            break;
        case 5:
            if (nu <= nu_a) {
                return (nu / nu_a) * (nu / nu_a);
            }
            return (p - 1) * C2_ * fastPow(nu / nu_a, -p / 2) * fastExp(-nu / nu_M);

            break;
        case 6:
            if (nu <= nu_a) {
                return (nu / nu_a) * (nu / nu_a);
            }
            return C2_ * fastPow(nu / nu_a, -p / 2) * fastExp(-nu / nu_M);

            break;

        default:
            return 0;
            break;
    }
}

/********************************************************************************************************************
 * FUNCTION: syn_p_nu_peak(Real B, Real p)
 * DESCRIPTION: Computes the peak power per electron in the comoving frame based on magnetic field B and power-law
 *              index p.
 ********************************************************************************************************************/
constexpr double sqrt3_half = 1.73205080757 / 2;
Real syn_p_nu_peak(Real B, Real p) { return (p - 1) * B * (sqrt3_half * con::e3 / (con::me * con::c2)); }

/********************************************************************************************************************
 * FUNCTION: syn_nu(Real gamma, Real B)
 * DESCRIPTION: Computes the synchrotron frequency for an electron with Lorentz factor gamma in a magnetic field B.
 ********************************************************************************************************************/
Real syn_nu(Real gamma, Real B) {
    Real nu = 3 * con::e / (4 * con::pi * con::me * con::c) * B * gamma * gamma;
    return nu;
}

/********************************************************************************************************************
 * FUNCTION: syn_gamma(Real nu, Real B)
 * DESCRIPTION: Computes the electron Lorentz factor corresponding to a synchrotron frequency nu in a magnetic field B.
 ********************************************************************************************************************/
Real syn_gamma(Real nu, Real B) {
    Real gamma = std::sqrt((4 * con::pi * con::me * con::c / (3 * con::e)) * (nu / B));
    return gamma;
}

/********************************************************************************************************************
 * FUNCTION: syn_gamma_M(Real B, std::vector<InverseComptonY> const& Ys, Real p)
 * DESCRIPTION: Computes the maximum electron Lorentz factor (gamma_M) for synchrotron emission, iterating until
 *              the inverse Compton Y parameter converges.
 ********************************************************************************************************************/
Real syn_gamma_M(Real B, InverseComptonY const& Ys, Real p) {
    if (B == 0) {
        return std::numeric_limits<Real>::infinity();
    }
    Real Y0 = InverseComptonY::Y_Thompson(Ys);
    Real gamma_M = std::sqrt(6 * con::pi * con::e / (con::sigmaT * B * (1 + Y0)));
    Real Y1 = InverseComptonY::Y_tilt_gamma(Ys, gamma_M, p);

    for (; std::fabs((Y1 - Y0) / Y0) > 1e-5;) {
        gamma_M = std::sqrt(6 * con::pi * con::e / (con::sigmaT * B * (1 + Y1)));
        Y0 = Y1;
        Y1 = InverseComptonY::Y_tilt_gamma(Ys, gamma_M, p);
    }

    return gamma_M;
}

/********************************************************************************************************************
 * FUNCTION: syn_gamma_m(Real Gamma_rel, Real gamma_M, Real eps_e, Real p, Real xi)
 * DESCRIPTION: Computes the minimum electron Lorentz factor (gamma_m) for synchrotron emission based on the
 *              available energy, power-law index p, and fraction of electrons xi. Uses rootBisection for p=2.
 ********************************************************************************************************************/
Real syn_gamma_m(Real Gamma_rel, Real gamma_M, Real eps_e, Real p, Real xi) {
    Real gamma_bar_minus_1 = eps_e * (Gamma_rel - 1) * (con::mp / con::me) / xi;
    Real gamma_m_minus_1 = 1;
    if (p > 2) {
        gamma_m_minus_1 = (p - 2) / (p - 1) * gamma_bar_minus_1;
    } else if (p < 2) {
        // Handle non-relativistic limit when p < 2
        gamma_m_minus_1 = std::pow((2 - p) / (p - 1) * gamma_bar_minus_1 * std::pow(gamma_M, p - 2), 1 / (p - 1));
    } else {
        gamma_m_minus_1 = rootBisection(
            [=](Real x) -> Real {
                return (x * std::log(gamma_M) - (x + 1) * std::log(x) - gamma_bar_minus_1 - std::log(gamma_M));
            },
            0, gamma_M);
    }
    return gamma_m_minus_1 + 1;
}

/********************************************************************************************************************
 * FUNCTION: syn_gamma_c(Real t_com, Real B, std::vector<InverseComptonY> const& Ys, Real p)
 * DESCRIPTION: Computes the cooling electron Lorentz factor (gamma_c) based on the comoving time t_com, magnetic
 *              field B, and inverse Compton corrections. Iterates until convergence.
 ********************************************************************************************************************/
Real syn_gamma_c(Real t_com, Real B, InverseComptonY const& Ys, Real p) {
    // t_com = (6*pi*gamma*me*c^2) /(gamma^2*beta^2*sigma_T*c*B^2*(1 + Y_tilt))
    // Real gamma_c = 6 * con::pi * con::me * con::c / (con::sigmaT * B * B * (1 + Y_tilt) * t_com);

    Real Y0 = InverseComptonY::Y_Thompson(Ys);
    Real gamma_bar = (6 * con::pi * con::me * con::c / con::sigmaT) / (B * B * (1 + Y0) * t_com);
    Real gamma_c = (gamma_bar + std::sqrt(gamma_bar * gamma_bar + 4)) / 2;
    Real Y1 = InverseComptonY::Y_tilt_gamma(Ys, gamma_c, p);

    for (; std::fabs((Y1 - Y0) / Y0) > 1e-3;) {
        gamma_bar = (6 * con::pi * con::me * con::c / con::sigmaT) / (B * B * (1 + Y1) * t_com);
        gamma_c = (gamma_bar + std::sqrt(gamma_bar * gamma_bar + 4)) / 2;
        Y0 = Y1;
        Y1 = InverseComptonY::Y_tilt_gamma(Ys, gamma_c, p);
    }

    return gamma_c;
}

/********************************************************************************************************************
 * FUNCTION: syn_gamma_a(Real Gamma_rel, Real B, Real I_syn_peak, Real gamma_m, Real gamma_c)
 * DESCRIPTION: Computes the self-absorption electron Lorentz factor (gamma_a) by equating the synchrotron intensity
 *              to a black-body like spectrum. Uses rootBisection if necessary.
 ********************************************************************************************************************/
Real syn_gamma_a(Real Gamma_rel, Real B, Real I_syn_peak, Real gamma_m, Real gamma_c) {
    Real gamma_peak = std::min(gamma_m, gamma_c);
    Real nu_peak = syn_nu(gamma_peak, B);
    Real ad_idx = adiabaticIndex(Gamma_rel);

    Real kT = (gamma_peak - 1) * (con::me * con::c2) * (ad_idx - 1);
    // 2kT(nu_a/c)^2 = I_peak*(nu_a/nu_peak)^(1/3)
    Real nu_a = fastPow(I_syn_peak * con::c2 / (std::cbrt(nu_peak) * 2 * kT), 0.6);

    // nu_peak is not real peak, peak at nu_a; kT = (gamma_a-1) * me *c^2*(ad_idx-1), I_syn = I_peak;
    // strong absorption
    if (nu_a > nu_peak) {
        nu_a = fastPow(
            I_syn_peak / (2 * con::me * (ad_idx - 1) * sqrt((4 * con::pi * con::me * con::c / (3 * con::e)) / B)), 0.4);
        /* Real gamma_a = syn_gamma(nu_a, B);
        if (gamma_a > 10) {
            return gamma_a;
        } else {
            double nu_max = syn_nu(10, B);
            double nu_min = syn_nu(1, B);
            Real C0 = std::sqrt((4 * con::pi * con::me * con::c / (3 * con::e)) / B);
            Real C1 = I_syn_peak / (2 * con::me * (ad_idx - 1));
            nu_a = rootBisection([=](Real x) -> Real { return C0 * x * x * x * x * x - x * x * x * x - C1; },
                                 std::sqrt(nu_min), std::sqrt(nu_max), 1e-3);
            nu_a *= nu_a;
        }*/
    }
    return syn_gamma(nu_a, B) + 1;
}

/********************************************************************************************************************
 * FUNCTION: syn_gamma_N_peak(Real gamma_a, Real gamma_m, Real gamma_c)
 * DESCRIPTION: Determines the electron Lorentz factor at which the column density peaks.
 ********************************************************************************************************************/
Real syn_gamma_N_peak(Real gamma_a, Real gamma_m, Real gamma_c) {
    Real gamma_peak = std::min(gamma_m, gamma_c);
    if (gamma_a > gamma_c) {
        return gamma_a;
    } else {
        return gamma_peak;
    }
}

/********************************************************************************************************************
 * FUNCTION: syn_gamma_N_peak(SynElectrons const& e)
 * DESCRIPTION: Overloaded function to determine the electron Lorentz factor at which the column density peaks,
 *              using a SynElectrons object.
 ********************************************************************************************************************/
Real syn_gamma_N_peak(SynElectrons const& e) { return syn_gamma_N_peak(e.gamma_a, e.gamma_m, e.gamma_c); }

/********************************************************************************************************************
 * FUNCTION: updateElectrons4Y(SynElectronGrid& e, Shock const& shock)
 * DESCRIPTION: Updates electron properties in the SynElectronGrid based on new inverse Compton Y parameter values
 *              and shock parameters.
 ********************************************************************************************************************/
void updateElectrons4Y(SynElectronGrid& e, Shock const& shock) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < t_size; ++k) {
                Real Gamma_rel = shock.Gamma_rel[i][j][k];
                Real t_com = shock.t_com[i][j][k];
                Real B = shock.B[i][j][k];
                Real p = e[i][j][k].p;
                auto& Ys = e[i][j][k].Ys;

                auto& electron = e[i][j][k];

                electron.gamma_M = syn_gamma_M(B, Ys, p);         // Update maximum electron Lorentz factor
                electron.gamma_c = syn_gamma_c(t_com, B, Ys, p);  // Update cooling electron Lorentz factor
                electron.gamma_a = syn_gamma_a(Gamma_rel, B, electron.I_nu_peak, electron.gamma_m, electron.gamma_c);
                electron.regime = getRegime(electron.gamma_a, electron.gamma_c, electron.gamma_m);
                electron.Y_c = InverseComptonY::Y_tilt_gamma(Ys, electron.gamma_c, p);
            }
        }
    }
}

/********************************************************************************************************************
 * FUNCTION: genSynElectrons(Shock const& shock, Real p, Real xi)
 * DESCRIPTION: Generates a SynElectronGrid based on the shock parameters, power-law index p, and electron
 *              partition factor xi.
 ********************************************************************************************************************/
SynElectronGrid genSynElectrons(Shock const& shock, Real p, Real xi) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    SynElectronGrid electrons = createSynElectronGrid(phi_size, theta_size, t_size);

    constexpr Real gamma_syn_limit = 3;

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < t_size; ++k) {
                Real Gamma_rel = shock.Gamma_rel[i][j][k];
                Real t_com = shock.t_com[i][j][k];
                Real B = shock.B[i][j][k];
                Real Sigma = shock.column_num_den[i][j][k];

                auto& e = electrons[i][j][k];

                e.gamma_M = syn_gamma_M(B, electrons[i][j][k].Ys, p);
                e.gamma_m = syn_gamma_m(Gamma_rel, e.gamma_M, shock.eps_e, p, xi);
                // Fraction of synchrotron electrons; the rest are cyclotron
                Real f = 1.;
                if (1 < e.gamma_m && e.gamma_m < gamma_syn_limit) {
                    f = std::min(fastPow((gamma_syn_limit - 1) / (e.gamma_m - 1), 1 - p), 1_r);
                    e.gamma_m = gamma_syn_limit;
                }
                e.column_num_den = Sigma * f;
                e.I_nu_peak = syn_p_nu_peak(B, p) * e.column_num_den / (4 * con::pi);
                e.gamma_c = syn_gamma_c(t_com, B, electrons[i][j][k].Ys, p);
                e.gamma_a = syn_gamma_a(Gamma_rel, B, e.I_nu_peak, e.gamma_m, e.gamma_c);
                e.regime = getRegime(e.gamma_a, e.gamma_c, e.gamma_m);
                e.p = p;
            }
        }
    }
    return electrons;
}

/********************************************************************************************************************
 * FUNCTION: genSynPhotons(Shock const& shock, SynElectronGrid const& e)
 * DESCRIPTION: Generates a SynPhotonGrid based on the shock parameters and the precomputed SynElectronGrid.
 *              For each grid cell, synchrotron photon frequencies are computed and internal constants are updated.
 ********************************************************************************************************************/
SynPhotonGrid genSynPhotons(Shock const& shock, SynElectronGrid const& e) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    SynPhotonGrid ph = createSynPhotonGrid(phi_size, theta_size, t_size);

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < t_size; ++k) {
                ph[i][j][k].e = &(e[i][j][k]);
                Real B = shock.B[i][j][k];

                ph[i][j][k].nu_M = syn_nu(e[i][j][k].gamma_M, B);
                ph[i][j][k].nu_m = syn_nu(e[i][j][k].gamma_m, B);
                ph[i][j][k].nu_c = syn_nu(e[i][j][k].gamma_c, B);
                ph[i][j][k].nu_a = syn_nu(e[i][j][k].gamma_a, B);

                ph[i][j][k].updateConstant();
            }
        }
    }
    return ph;
}
