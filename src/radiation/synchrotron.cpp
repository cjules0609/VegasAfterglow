//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "synchrotron.h"

#include <cmath>

#include "afterglow.h"
#include "inverse-compton.h"
#include "macros.h"
#include "physics.h"
#include "utilities.h"

Real SynElectrons::compute_column_num_den(Real gamma) const {
    if (gamma < gamma_c) {  // Below cooling Lorentz factor: direct scaling
        return column_num_den * compute_gamma_spectrum(gamma) * fast_exp2(-gamma / gamma_M);
    } else {
        return column_num_den * compute_gamma_spectrum(gamma) * fast_exp2(-gamma / gamma_M) * (1 + Y_c) /
               (1 + InverseComptonY::compute_Y_tilt_at_gamma(Ys, gamma, p));
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Helper function that checks if three values are in non-decreasing order.
 * @param a First value
 * @param b Middle value
 * @param c Last value
 * @return True if a ≤ b ≤ c, false otherwise
 * <!-- ************************************************************************************** -->
 */
inline bool order(Real a, Real b, Real c) { return a <= b && b <= c; };

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Determines the spectral regime (1-6) based on the ordering of characteristic Lorentz factors.
 * @details Classifies the regime based on the ordering of absorption (a), cooling (c),
 *          and minimum (m) Lorentz factors.
 * @param a Absorption Lorentz factor
 * @param c Cooling Lorentz factor
 * @param m Minimum Lorentz factor
 * @return Regime number (1-6) or 0 if no valid regime is found
 * <!-- ************************************************************************************** -->
 */
size_t determine_regime(Real a, Real c, Real m) {
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

Real SynElectrons::compute_gamma_spectrum(Real gamma) const {
    switch (regime) {
        case 1:  // same as case 2
        case 2:
            if (gamma <= gamma_m) {
                return 0;  // Below minimum Lorentz factor, spectrum is zero
            } else if (gamma <= gamma_c) {
                return (p - 1) * fast_pow(gamma / gamma_m, -p) /
                       gamma_m;  // Power-law spectrum between gamma_m and gamma_c
            } else
                return (p - 1) * fast_pow(gamma / gamma_m, -p) * gamma_c / (gamma * gamma_m);
            // Above cooling Lorentz factor: exponential cutoff applied

            break;
        case 3:
            if (gamma <= gamma_c) {
                return 0;  // Below cooling Lorentz factor, spectrum is zero
            } else if (gamma <= gamma_m) {
                return gamma_c / (gamma * gamma);  // Intermediate regime scaling
            } else
                return gamma_c / (gamma * gamma_m) * fast_pow(gamma / gamma_m, -p);
            // Above minimum Lorentz factor: power-law with exponential cutoff

            break;
        case 4:  // Gao, Lei, Wu and Zhang 2013 Eq 18
            if (gamma <= gamma_a) {
                return 3 * gamma * gamma / (gamma_a * gamma_a * gamma_a);  // Rising part of the spectrum
            } else if (gamma <= gamma_m) {
                return gamma_c / (gamma * gamma);  // Transition region
            } else
                return gamma_c / (gamma * gamma_m) * fast_pow(gamma / gamma_m, -p);
            // High energy tail with exponential cutoff

            break;
        case 5:  // Gao, Lei, Wu and Zhang 2013 Eq 19
            if (gamma <= gamma_a) {
                return 3 * gamma * gamma / (gamma_a * gamma_a * gamma_a);  // Rising part of the spectrum
            } else
                return (p - 1) * gamma_c / (gamma * gamma_m) * fast_pow(gamma / gamma_m, -p);
            // Power-law decay with exponential cutoff

            break;
        case 6:  // Gao, Lei, Wu and Zhang 2013 Eq 20
            if (gamma <= gamma_a) {
                return 3 * gamma * gamma / (gamma_a * gamma_a * gamma_a);  // Rising part of the spectrum
            } else
                return fast_pow(gamma_m, p - 1) * gamma_c * fast_pow(gamma, -(p + 1));
            // Steeper decay in this regime

            break;
        default:
            return 0;
    }
}

Real SynPhotons::compute_I_nu(Real nu) const {
    if (nu < nu_c) {  // Below cooling frequency, simple scaling
        return elec->I_nu_peak * compute_spectrum(nu) * fast_exp2(-nu / nu_M);
    } else {
        return elec->I_nu_peak * compute_spectrum(nu) * fast_exp2(-nu / nu_M) * (1 + elec->Y_c) /
               (1 + InverseComptonY::compute_Y_tilt_at_nu(elec->Ys, nu, elec->p));
    }
}

Real SynPhotons::compute_log2_I_nu(Real log2_nu) const {
    if (log2_nu < log2_nu_c) {  // Below cooling frequency, simple scaling
        return log2_I_nu_peak + compute_log2_spectrum(log2_nu) - fast_exp2(log2_nu) / nu_M;
    } else {
        return log2_I_nu_peak + compute_log2_spectrum(log2_nu) - fast_exp2(log2_nu) / nu_M + fast_log2(1 + elec->Y_c) -
               fast_log2(1 + InverseComptonY::compute_Y_tilt_at_nu(elec->Ys, std::exp2(log2_nu), elec->p));
    }
}

void SynPhotons::update_constant() {
    // Update constants based on spectral parameters
    Real p = elec->p;
    if (elec->regime == 1) {
        // a_m_1_3 = std::cbrt(nu_a / nu_m);  // (nu_a / nu_m)^(1/3)
        // c_m_mpa1_2 = fastPow(nu_c / nu_m, (-p + 1) / 2);  // (nu_c / nu_m)^((-p+1)/2)
        C1_ = std::cbrt(nu_a / nu_m);
        C2_ = fast_pow(nu_c / nu_m, (-p + 1) / 2);

        log2_C1_ = fast_log2(nu_a / nu_m) / 3 - 2 * fast_log2(nu_a);
        log2_C2_ = -fast_log2(nu_m) / 3;
        log2_C3_ = (p - 1) / 2 * fast_log2(nu_m);
        log2_C4_ = (p - 1) / 2 * fast_log2(nu_m / nu_c) + p / 2 * fast_log2(nu_c);
    } else if (elec->regime == 2) {
        // m_a_pa4_2 = fastPow(nu_m / nu_a, (p + 4) / 2);    // (nu_m / nu_a)^((p+4)/2)
        // a_m_mpa1_2 = fastPow(nu_a / nu_m, (-p + 1) / 2);  // (nu_a / nu_m)^((-p+1)/2)
        // c_m_mpa1_2 = fastPow(nu_c / nu_m, (-p + 1) / 2);  // (nu_c / nu_m)^((-p+1)/2)
        C1_ = fast_pow(nu_m / nu_a, (p + 4) / 2);
        C2_ = fast_pow(nu_a / nu_m, (-p + 1) / 2);
        C3_ = fast_pow(nu_c / nu_m, (-p + 1) / 2);

        log2_C1_ = (p + 4) / 2 * fast_log2(nu_m / nu_a) - 2 * fast_log2(nu_m);
        log2_C2_ = (p - 1) / 2 * fast_log2(nu_m / nu_a) - 2.5 * fast_log2(nu_a);
        log2_C3_ = (p - 1) / 2 * fast_log2(nu_m);
        log2_C4_ = (p - 1) / 2 * fast_log2(nu_m / nu_c) + p / 2 * fast_log2(nu_c);
    } else if (elec->regime == 3) {
        // a_c_1_3 = std::cbrt(nu_a / nu_c);  // (nu_a / nu_c)^(1/3)
        // c_m_1_2 = std::sqrt(nu_c / nu_m);  // (nu_c / nu_m)^(1/2)
        C1_ = std::cbrt(nu_a / nu_c);
        C2_ = std::sqrt(nu_c / nu_m);

        log2_C1_ = fast_log2(nu_a / nu_c) / 3 - 2 * fast_log2(nu_a);
        log2_C2_ = -fast_log2(nu_c) / 3;
        log2_C3_ = fast_log2(nu_c) / 2;
        log2_C4_ = fast_log2(nu_c / nu_m) / 2 + p / 2 * fast_log2(nu_m);
    } else if (elec->regime == 4) {
        // a_m_1_2 = std::sqrt(nu_a / nu_m);  // (nu_a / nu_m)^(1/2)
        // R4 = std::sqrt(nu_c / nu_a) / 3;   // (nu_c / nu_a)^(1/2) / 3; // R4: scaling factor for regime 4
        C1_ = std::sqrt(nu_a / nu_m);
        C2_ = std::sqrt(nu_c / nu_a) / 3;

        log2_C1_ = -2 * fast_log2(nu_a);
        log2_C2_ = fast_log2(C2_) + fast_log2(nu_a) / 2;
        log2_C3_ = fast_log2(C2_) + fast_log2(nu_a / nu_m) / 2 + p / 2 * fast_log2(nu_m);
    } else if (elec->regime == 5 || elec->regime == 6) {
        // R4 = std::sqrt(nu_c / nu_a) / 3;              // (nu_c / nu_a)^(1/2) / 3; // R4: scaling factor for
        // regime 4 R6 = R4 * fastPow(nu_m / nu_a, (p - 1) / 2);  // R6: scaling factor for regime 6
        C1_ = std::sqrt(nu_c / nu_a) / 3;
        C2_ = C1_ * fast_pow(nu_m / nu_a, (p - 1) / 2);

        log2_C1_ = -2 * fast_log2(nu_a);
        log2_C2_ = fast_log2((p - 1) * C2_) + p / 2 * fast_log2(nu_a);
        log2_C3_ = fast_log2(C2_) + p / 2 * fast_log2(nu_a);
    }
    // R5 = (p - 1) * R6;  // R5 scales R6 (commented out)
}

Real SynPhotons::compute_spectrum(Real nu) const {
    Real p = elec->p;
    switch (elec->regime) {
        case 1:
            if (nu <= nu_a) {
                return C1_ * (nu / nu_a) * (nu / nu_a);
            }
            if (nu <= nu_m) {
                return std::cbrt(nu / nu_m);
            }
            if (nu <= nu_c) {
                return fast_pow(nu / nu_m, -(p - 1) / 2);
            }

            return C2_ * fast_pow(nu / nu_c, -p / 2);

            break;
        case 2:
            if (nu <= nu_m) {
                return C1_ * (nu / nu_m) * (nu / nu_m);
            }
            if (nu <= nu_a) {
                return C2_ * pow52(nu / nu_a);  // Using pow52 for (nu / nu_a)^(5/2)
            }
            if (nu <= nu_c) {
                return fast_pow(nu / nu_m, -(p - 1) / 2);
            }

            return C3_ * fast_pow(nu / nu_c, -p / 2);

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
            return C2_ * fast_pow(nu / nu_m, -p / 2);

            break;
        case 4:
            if (nu <= nu_a) {
                return (nu / nu_a) * (nu / nu_a);
            }
            if (nu <= nu_m) {
                return C2_ * std::sqrt(nu_a / nu);
            }
            return C2_ * C1_ * fast_pow(nu / nu_m, -p / 2);

            break;
        case 5:
            if (nu <= nu_a) {
                return (nu / nu_a) * (nu / nu_a);
            }
            return (p - 1) * C2_ * fast_pow(nu / nu_a, -p / 2);

            break;
        case 6:
            if (nu <= nu_a) {
                return (nu / nu_a) * (nu / nu_a);
            }
            return C2_ * fast_pow(nu / nu_a, -p / 2);

            break;

        default:
            return 0;
            break;
    }
}

Real SynPhotons::compute_log2_spectrum(Real log2_nu) const {
    Real p = elec->p;
    switch (elec->regime) {
        case 1:
            if (log2_nu <= log2_nu_a) {
                return log2_C1_ + 2 * log2_nu;
            }
            if (log2_nu <= log2_nu_m) {
                return log2_C2_ + log2_nu / 3;
            }
            if (log2_nu <= log2_nu_c) {
                return log2_C3_ - (p - 1) / 2 * log2_nu;
            }
            return log2_C4_ - p / 2 * log2_nu;

            break;
        case 2:
            if (log2_nu <= log2_nu_m) {
                return log2_C1_ + 2 * log2_nu;
            }
            if (log2_nu <= log2_nu_a) {
                return log2_C2_ + 2.5 * log2_nu;
            }
            if (log2_nu <= log2_nu_c) {
                return log2_C3_ - (p - 1) / 2 * log2_nu;
            }

            return log2_C4_ - p / 2 * log2_nu;

            break;
        case 3:
            if (log2_nu <= log2_nu_a) {
                return log2_C1_ + 2 * log2_nu;
            }
            if (log2_nu <= log2_nu_c) {
                return log2_C2_ + log2_nu / 3;
            }
            if (log2_nu <= log2_nu_m) {
                return log2_C3_ - log2_nu / 2;
            }

            return log2_C4_ - p / 2 * log2_nu;

            break;
        case 4:
            if (log2_nu <= log2_nu_a) {
                return log2_C1_ + 2 * log2_nu;
            }
            if (log2_nu <= log2_nu_m) {
                return log2_C2_ - log2_nu / 2;
            }

            return log2_C3_ - p / 2 * log2_nu;

            break;
        case 5:
            if (log2_nu <= log2_nu_a) {
                return log2_C1_ + 2 * log2_nu;
            }

            return log2_C2_ - p / 2 * log2_nu;

            break;
        case 6:
            if (log2_nu <= log2_nu_a) {
                return log2_C1_ + 2 * log2_nu;
            }

            return log2_C3_ - p / 2 * log2_nu;

            break;

        default:
            return -con::inf;
            break;
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the peak synchrotron power per electron in the comoving frame.
 * @details Based on magnetic field strength B and power-law index p of the electron distribution.
 * @param B Magnetic field strength
 * @param p Power-law index of electron distribution
 * @return Peak synchrotron power per electron
 * <!-- ************************************************************************************** -->
 */
Real compute_syn_peak_power(Real B, Real p) {
    constexpr double sqrt3_half = 1.73205080757 / 2;
    return (p - 1) * B * (sqrt3_half * con::e3 / (con::me * con::c2));
    // return B * (1.73205080757 * con::e3 / (con::me * con::c2));
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the characteristic synchrotron frequency for electrons with a given Lorentz factor.
 * @details Uses the standard synchrotron formula.
 * @param gamma Electron Lorentz factor
 * @param B Magnetic field strength
 * @return Characteristic synchrotron frequency
 * <!-- ************************************************************************************** -->
 */
Real compute_syn_freq(Real gamma, Real B) {
    Real nu = 3 * con::e / (4 * con::pi * con::me * con::c) * B * gamma * gamma;
    return nu;
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the electron Lorentz factor corresponding to a synchrotron frequency.
 * @details Inverse of the compute_syn_freq function.
 * @param nu Synchrotron frequency
 * @param B Magnetic field strength
 * @return Corresponding electron Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real compute_syn_gamma(Real nu, Real B) {
    Real gamma = std::sqrt((4 * con::pi * con::me * con::c / (3 * con::e)) * (nu / B));
    return gamma;
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the maximum electron Lorentz factor for synchrotron emission.
 * @details Uses an iterative approach to account for inverse Compton cooling effects.
 * @param B Magnetic field strength
 * @param Ys InverseComptonY object
 * @param p Spectral index of electron distribution
 * @return Maximum electron Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real compute_syn_gamma_M(Real B, InverseComptonY const& Ys, Real p) {
    if (B == 0) {
        return std::numeric_limits<Real>::infinity();
    }
    Real Y0 = InverseComptonY::compute_Y_Thompson(Ys);
    Real gamma_M = std::sqrt(6 * con::pi * con::e / con::sigmaT / (B * (1 + Y0)));
    Real Y1 = InverseComptonY::compute_Y_tilt_at_gamma(Ys, gamma_M, p);

    for (; std::fabs((Y1 - Y0) / Y0) > 1e-3;) {
        gamma_M = std::sqrt(6 * con::pi * con::e / con::sigmaT / (B * (1 + Y1)));
        Y0 = Y1;
        Y1 = InverseComptonY::compute_Y_tilt_at_gamma(Ys, gamma_M, p);
    }

    return gamma_M;
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the minimum electron Lorentz factor for synchrotron emission.
 * @details Accounts for different power-law indices with special handling for the p=2 case.
 *          Uses the fraction of shock energy given to electrons (eps_e) and electron fraction (xi).
 * @param Gamma_rel Bulk Lorentz factor of the shock
 * @param gamma_M Maximum electron Lorentz factor
 * @param eps_e Fraction of shock energy given to electrons
 * @param p Power-law index of electron distribution
 * @param xi Fraction of electrons accelerated
 * @return Minimum electron Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real compute_syn_gamma_m(Real Gamma_rel, Real gamma_M, Real eps_e, Real p, Real xi) {
    Real gamma_bar_minus_1 = eps_e * (Gamma_rel - 1) * (con::mp / con::me) / xi;
    Real gamma_m_minus_1 = 1;
    if (p > 2) {
        gamma_m_minus_1 = (p - 2) / (p - 1) * gamma_bar_minus_1;
    } else if (p < 2) {
        gamma_m_minus_1 = std::pow((2 - p) / (p - 1) * gamma_bar_minus_1 * std::pow(gamma_M, p - 2), 1 / (p - 1));
    } else {
        gamma_m_minus_1 = root_bisect(
            [=](Real x) -> Real {
                return (x * std::log(gamma_M) - (x + 1) * std::log(x) - gamma_bar_minus_1 - std::log(gamma_M));
            },
            0, gamma_M);
    }
    return gamma_m_minus_1 + 1;
}

Real compute_gamma_c(Real t_comv, Real B, InverseComptonY const& Ys, Real p) {
    // t_com = (6*pi*gamma*me*c^2) /(gamma^2*beta^2*sigma_T*c*B^2*(1 + Y_tilt))
    // Real gamma_c = 6 * con::pi * con::me * con::c / (con::sigmaT * B * B * (1 + Y_tilt) * t_comv);

    Real Y0 = InverseComptonY::compute_Y_Thompson(Ys);
    Real gamma_bar = (6 * con::pi * con::me * con::c / con::sigmaT) / (B * B * (1 + Y0) * t_comv);
    Real gamma_c = (gamma_bar + std::sqrt(gamma_bar * gamma_bar + 4)) / 2;  // correction on newtonian regime

    Real Y1 = InverseComptonY::compute_Y_tilt_at_gamma(Ys, gamma_c, p);
    for (; std::fabs((Y1 - Y0) / Y0) > 1e-3;) {  // iterate for IC cooling
        gamma_bar = (6 * con::pi * con::me * con::c / con::sigmaT) / (B * B * (1 + Y1) * t_comv);
        gamma_c = (gamma_bar + std::sqrt(gamma_bar * gamma_bar + 4)) / 2;
        Y0 = Y1;
        Y1 = InverseComptonY::compute_Y_tilt_at_gamma(Ys, gamma_c, p);
    }

    return gamma_c;
}

/**
 * <!-- ************************************************************************************** -->
 * @internal
 * @brief Calculates the self-absorption Lorentz factor by equating synchrotron emission to blackbody.
 * @details Uses the peak intensity and shock parameters to determine where absorption becomes important.
 *          Handles both weak and strong absorption regimes.
 * @param Gamma_rel Bulk Lorentz factor of the shock
 * @param B Magnetic field strength
 * @param I_syn_peak Peak synchrotron intensity
 * @param gamma_m Minimum electron Lorentz factor
 * @param gamma_c Cooling electron Lorentz factor
 * @param gamma_M Maximum electron Lorentz factor
 * @param p Power-law index of electron distribution
 * @return Self-absorption Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real compute_syn_gamma_a(Real Gamma_rel, Real B, Real I_syn_peak, Real gamma_m, Real gamma_c, Real gamma_M, Real p) {
    Real gamma_peak = std::min(gamma_m, gamma_c);
    Real nu_peak = compute_syn_freq(gamma_peak, B);

    Real kT = (gamma_peak - 1) * (con::me * con::c2) / 3;
    // 2kT(nu_a/c)^2 = I_peak*(nu_a/nu_peak)^(1/3) // first assume nu_a is in the 1/3 segment
    Real nu_a = fast_pow(I_syn_peak * con::c2 / (std::cbrt(nu_peak) * 2 * kT), 0.6);

    if (nu_a > nu_peak) {  // nu_a is not in the 1/3 segment
        constexpr Real coef = 3 * con::e / (4 * con::pi * con::me * con::c);
        if (gamma_c > gamma_m) {  // then assume nu_a is in the -(p-1)/2 segment
            Real nu_m = compute_syn_freq(gamma_m, B);
            nu_a = fast_pow(I_syn_peak * con::c2 / (2 * kT) * fast_pow(nu_m, p / 2), 2 / (p + 4));
            Real nu_c = compute_syn_freq(gamma_c, B);
            if (nu_a > nu_c) {  // nu_a is not in the -(p-1)/2 segment, strong absorption
                Real C = 1.5 * I_syn_peak / (con::me * coef * coef * B * B);
                Real gamma_a = root_bisect([C](Real x) -> Real { return x * x * x * x * x - x - C; }, gamma_c, gamma_M);
                return gamma_a;
            }
        } else {  // strong absorption
            Real C = 1.5 * I_syn_peak / (con::me * coef * coef * B * B);
            Real gamma_a = root_bisect([C](Real x) -> Real { return x * x * x * x * x - x - C; }, gamma_c, gamma_M);
            return gamma_a;
        }
    }
    return compute_syn_gamma(nu_a, B) + 1;
}

Real compute_gamma_peak(Real gamma_a, Real gamma_m, Real gamma_c) {
    Real gamma_peak = std::min(gamma_m, gamma_c);
    if (gamma_a > gamma_c) {
        return gamma_a;
    } else {
        return gamma_peak;
    }
}

Real cyclotron_correction(Real& gamma_m, Real& gamma_c, Real p) {
    constexpr Real gamma_cyclotron = 2;
    Real f = 1.;

    if (1 < gamma_m && gamma_m <= gamma_cyclotron) {
        f = (gamma_m - 1) / (gamma_cyclotron - 1);
        if (p > 3) {
            f = fast_pow((gamma_m - 1) / (gamma_cyclotron - 1), (p - 1) / 2);
        }
        gamma_m = gamma_cyclotron;
    }

    return f;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Determines the peak Lorentz factor directly from a SynElectrons object.
 * @details Convenient wrapper around the three-parameter version.
 * @param e Synchrotron electron object
 * @return Peak Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real compute_gamma_peak(SynElectrons const& e) { return compute_gamma_peak(e.gamma_a, e.gamma_m, e.gamma_c); }

void update_electrons_4Y(SynElectronGrid& electrons, Shock const& shock) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            size_t k_inj = shock.injection_idx(i, j);
            for (size_t k = 0; k < t_size; ++k) {
                if (shock.required(i, j, k) == 0) {
                    continue;
                }
                Real Gamma_rel = shock.Gamma_rel(i, j, k);
                Real t_com = shock.t_comv(i, j, k);
                Real B = shock.B(i, j, k);
                Real p = electrons(i, j, k).p;
                auto& Ys = electrons(i, j, k).Ys;
                auto& elec = electrons(i, j, k);

                elec.gamma_M = compute_syn_gamma_M(B, Ys, p);  // Update maximum electron Lorentz factor
                if (k <= k_inj) {
                    elec.gamma_c = compute_gamma_c(t_com, B, Ys, p);  // Update cooling electron Lorentz factor
                } else {  // no shocked electron injection, just adiabatic cooling
                    elec.gamma_c = electrons(i, j, k_inj).gamma_c * elec.gamma_m / electrons(i, j, k_inj).gamma_m;
                    elec.gamma_M = elec.gamma_c;
                }
                elec.gamma_a =
                    compute_syn_gamma_a(Gamma_rel, B, elec.I_nu_peak, elec.gamma_m, elec.gamma_c, elec.gamma_M, p);
                elec.regime = determine_regime(elec.gamma_a, elec.gamma_c, elec.gamma_m);
                elec.Y_c = InverseComptonY::compute_Y_tilt_at_gamma(Ys, elec.gamma_c, p);
            }
        }
    }
}

SynElectronGrid generate_syn_electrons(Shock const& shock) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    RadParams rad = shock.rad;

    SynElectronGrid electrons({phi_size, theta_size, t_size});

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            size_t k_inj = shock.injection_idx(i, j);
            for (size_t k = 0; k < t_size; ++k) {
                if (shock.required(i, j, k) == 0) {
                    continue;
                }
                Real Gamma_rel = shock.Gamma_rel(i, j, k);
                Real t_com = shock.t_comv(i, j, k);
                Real B = shock.B(i, j, k);

                auto& elec = electrons(i, j, k);

                elec.gamma_M = compute_syn_gamma_M(B, elec.Ys, rad.p);
                elec.gamma_m = compute_syn_gamma_m(Gamma_rel, elec.gamma_M, rad.eps_e, rad.p, rad.xi_e);

                // Fraction of synchrotron electrons; the rest are cyclotron
                Real f = cyclotron_correction(elec.gamma_m, elec.gamma_c, rad.p);

                elec.column_num_den = shock.proton_column_num_den(i, j, k) * rad.xi_e * f;
                elec.I_nu_peak = compute_syn_peak_power(B, rad.p) * elec.column_num_den / (4 * con::pi);
                if (k <= k_inj) {
                    elec.gamma_c = compute_gamma_c(t_com, B, elec.Ys, rad.p);
                } else {  // no shocked electron injection, just adiabatic cooling
                    elec.gamma_c = electrons(i, j, k_inj).gamma_c * elec.gamma_m / electrons(i, j, k_inj).gamma_m;
                    elec.gamma_M = elec.gamma_c;
                }
                elec.gamma_a =
                    compute_syn_gamma_a(Gamma_rel, B, elec.I_nu_peak, elec.gamma_m, elec.gamma_c, elec.gamma_M, rad.p);
                elec.regime = determine_regime(elec.gamma_a, elec.gamma_c, elec.gamma_m);
                elec.p = rad.p;
            }
        }
    }
    return electrons;
}

void generate_syn_electrons(SynElectronGrid& electrons, Shock const& shock) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    RadParams rad = shock.rad;

    electrons.resize({phi_size, theta_size, t_size});

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            size_t k_inj = shock.injection_idx(i, j);
            for (size_t k = 0; k < t_size; ++k) {
                if (shock.required(i, j, k) == 0) {
                    continue;
                }
                Real Gamma_rel = shock.Gamma_rel(i, j, k);
                Real t_com = shock.t_comv(i, j, k);
                Real B = shock.B(i, j, k);

                auto& elec = electrons(i, j, k);

                elec.gamma_M = compute_syn_gamma_M(B, elec.Ys, rad.p);
                elec.gamma_m = compute_syn_gamma_m(Gamma_rel, elec.gamma_M, rad.eps_e, rad.p, rad.xi_e);

                // Fraction of synchrotron electrons; the rest are cyclotron
                Real f = cyclotron_correction(elec.gamma_m, elec.gamma_c, rad.p);

                elec.column_num_den = shock.proton_column_num_den(i, j, k) * rad.xi_e * f;
                elec.I_nu_peak = compute_syn_peak_power(B, rad.p) * elec.column_num_den / (4 * con::pi);
                if (k <= k_inj) {
                    elec.gamma_c = compute_gamma_c(t_com, B, electrons(i, j, k).Ys, rad.p);
                } else {  // no shocked electron injection, just adiabatic cooling
                    elec.gamma_c = electrons(i, j, k_inj).gamma_c * elec.gamma_m / electrons(i, j, k_inj).gamma_m;
                    elec.gamma_M = elec.gamma_c;
                }
                elec.gamma_a =
                    compute_syn_gamma_a(Gamma_rel, B, elec.I_nu_peak, elec.gamma_m, elec.gamma_c, elec.gamma_M, rad.p);
                elec.regime = determine_regime(elec.gamma_a, elec.gamma_c, elec.gamma_m);

                elec.p = rad.p;
            }
        }
    }
}

SynPhotonGrid generate_syn_photons(Shock const& shock, SynElectronGrid const& electrons) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    SynPhotonGrid photons({phi_size, theta_size, t_size});

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < t_size; ++k) {
                auto& ph = photons(i, j, k);
                auto& elec = electrons(i, j, k);
                ph.elec = &elec;
                if (shock.required(i, j, k) == 0) {
                    continue;
                }
                Real B = shock.B(i, j, k);

                ph.nu_M = compute_syn_freq(elec.gamma_M, B);
                ph.nu_m = compute_syn_freq(elec.gamma_m, B);
                ph.nu_c = compute_syn_freq(elec.gamma_c, B);
                ph.nu_a = compute_syn_freq(elec.gamma_a, B);

                ph.log2_I_nu_peak = fast_log2(elec.I_nu_peak);
                ph.log2_nu_m = fast_log2(ph.nu_m);
                ph.log2_nu_c = fast_log2(ph.nu_c);
                ph.log2_nu_a = fast_log2(ph.nu_a);
                ph.log2_nu_M = fast_log2(ph.nu_M);
                ph.update_constant();
            }
        }
    }
    return photons;
}

void generate_syn_photons(SynPhotonGrid& photons, Shock const& shock, SynElectronGrid const& electrons) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    photons.resize({phi_size, theta_size, t_size});
    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < t_size; ++k) {
                auto& ph = photons(i, j, k);
                auto& elec = electrons(i, j, k);
                ph.elec = &elec;
                if (shock.required(i, j, k) == 0) {
                    continue;
                }
                Real B = shock.B(i, j, k);

                ph.nu_M = compute_syn_freq(elec.gamma_M, B);
                ph.nu_m = compute_syn_freq(elec.gamma_m, B);
                ph.nu_c = compute_syn_freq(elec.gamma_c, B);
                ph.nu_a = compute_syn_freq(elec.gamma_a, B);

                ph.log2_I_nu_peak = fast_log2(elec.I_nu_peak);
                ph.log2_nu_m = fast_log2(ph.nu_m);
                ph.log2_nu_c = fast_log2(ph.nu_c);
                ph.log2_nu_a = fast_log2(ph.nu_a);
                ph.log2_nu_M = fast_log2(ph.nu_M);
                ph.update_constant();
            }
        }
    }
}
