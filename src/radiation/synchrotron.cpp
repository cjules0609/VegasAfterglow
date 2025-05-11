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

InverseComptonY::InverseComptonY(Real nu_m, Real nu_c, Real B, Real Y_T) noexcept {
    gamma_hat_m = con::me * con::c2 / con::h / nu_m;  // Compute minimum characteristic Lorentz factor
    gamma_hat_c = con::me * con::c2 / con::h / nu_c;  // Compute cooling characteristic Lorentz factor
    this->Y_T = Y_T;                                  // Set the Thomson Y parameter
    nu_hat_m = compute_syn_freq(gamma_hat_m, B);      // Compute corresponding synchrotron frequency for gamma_hat_m
    nu_hat_c = compute_syn_freq(gamma_hat_c, B);      // Compute corresponding synchrotron frequency for gamma_hat_c

    if (nu_hat_m <= nu_hat_c) {
        regime = 1;  // fast IC cooling regime
    } else {
        regime = 2;  // slow IC cooling regime
    }
}

InverseComptonY::InverseComptonY(Real Y_T) noexcept {
    this->Y_T = Y_T;  // Set the Thomson Y parameter
    regime = 3;       // Set regime to 3 (special case)
}

InverseComptonY::InverseComptonY() noexcept {
    nu_hat_m = 0;
    nu_hat_c = 0;
    gamma_hat_m = 0;
    gamma_hat_c = 0;
    Y_T = 0;
    regime = 0;
}

Real InverseComptonY::compute_val_at_gamma(Real gamma, Real p) const {
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
                return Y_T * fast_pow(gamma / gamma_hat_c, (p - 3) / 2);  // Scaling in intermediate regime
            } else
                return Y_T * pow43(gamma_hat_m / gamma) *
                       fast_pow(gamma_hat_m / gamma_hat_c, (p - 3) / 2);  // High gamma scaling

            break;
        default:
            return 0;
            break;
    }
}

Real InverseComptonY::compute_val_at_nu(Real nu, Real p) const {
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
                return Y_T * fast_pow(nu / nu_hat_c, (p - 3) / 4);  // Intermediate frequency scaling
            } else
                return Y_T * pow23(nu_hat_m / nu) *
                       fast_pow(nu_hat_m / nu_hat_c, (p - 3) / 4);  // High frequency scaling

            break;
        default:
            return 0;
            break;
    }
}

Real InverseComptonY::compute_Y_Thompson(InverseComptonY const& Ys) {
    /*Real Y_tilt = 0;
    for (auto& Y : Ys) {
        Y_tilt += Y.Y_T;  // Sum each object's Y_T
    }
    return Y_tilt;*/
    return Ys.Y_T;
}

Real InverseComptonY::compute_Y_tilt_at_gamma(InverseComptonY const& Ys, Real gamma, Real p) {
    /*Real Y_tilt = 0;
    for (auto& Y : Ys) {
        Y_tilt += Y.as_gamma(gamma, p);  // Sum effective Y parameters based on gamma
    }
    return Y_tilt;*/
    return Ys.compute_val_at_gamma(gamma, p);
}

Real InverseComptonY::compute_Y_tilt_at_nu(InverseComptonY const& Ys, Real nu, Real p) {
    /* Real Y_tilt = 0;
     for (auto& Y : Ys) {
         Y_tilt += Y.as_nu(nu, p);  // Sum effective Y parameters based on frequency
     }
     return Y_tilt;*/
    return Ys.compute_val_at_nu(nu, p);
}

Real SynElectrons::compute_column_num_den(Real gamma) const {
    if (gamma < gamma_c) {
        return column_num_den * compute_gamma_spectrum(gamma);  // Below cooling Lorentz factor: direct scaling
    } else {
        return column_num_den * compute_gamma_spectrum(gamma) * (1 + Y_c) /
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
                return (p - 1) * fast_pow(gamma / gamma_m, -p) * gamma_c / (gamma * gamma_m) *
                       fast_exp(-gamma / gamma_M);
            // Above cooling Lorentz factor: exponential cutoff applied

            break;
        case 3:
            if (gamma <= gamma_c) {
                return 0;  // Below cooling Lorentz factor, spectrum is zero
            } else if (gamma <= gamma_m) {
                return gamma_c / (gamma * gamma);  // Intermediate regime scaling
            } else
                return gamma_c / (gamma * gamma_m) * fast_pow(gamma / gamma_m, -p) * fast_exp(-gamma / gamma_M);
            // Above minimum Lorentz factor: power-law with exponential cutoff

            break;
        case 4:  // Gao, Lei, Wu and Zhang 2013 Eq 18
            if (gamma <= gamma_a) {
                return 3 * gamma * gamma / (gamma_a * gamma_a * gamma_a);  // Rising part of the spectrum
            } else if (gamma <= gamma_m) {
                return gamma_c / (gamma * gamma);  // Transition region
            } else
                return gamma_c / (gamma * gamma_m) * fast_pow(gamma / gamma_m, -p) * fast_exp(-gamma / gamma_M);
            // High energy tail with exponential cutoff

            break;
        case 5:  // Gao, Lei, Wu and Zhang 2013 Eq 19
            if (gamma <= gamma_a) {
                return 3 * gamma * gamma / (gamma_a * gamma_a * gamma_a);  // Rising part of the spectrum
            } else
                return (p - 1) * gamma_c / (gamma * gamma_m) * fast_pow(gamma / gamma_m, -p) *
                       fast_exp(-gamma / gamma_M);
            // Power-law decay with exponential cutoff

            break;
        case 6:  // Gao, Lei, Wu and Zhang 2013 Eq 20
            if (gamma <= gamma_a) {
                return 3 * gamma * gamma / (gamma_a * gamma_a * gamma_a);  // Rising part of the spectrum
            } else
                return fast_pow(gamma_m, p - 1) * gamma_c * fast_pow(gamma, -(p + 1)) * fast_exp(-gamma / gamma_M);
            // Steeper decay in this regime

            break;
        default:
            return 0;
    }
}

Real SynPhotons::compute_I_nu(Real nu) const {
    if (nu < nu_c) {
        return e->I_nu_peak * compute_spectrum(nu);  // Below cooling frequency, simple scaling
    } else {
        return e->I_nu_peak * compute_spectrum(nu) * (1 + e->Y_c) /
               (1 + InverseComptonY::compute_Y_tilt_at_nu(e->Ys, nu, e->p));
    }
}

Real SynPhotons::compute_log2_I_nu(Real log2_nu) const {
    if (log2_nu < log2_nu_c) {
        return log2_I_nu_peak + compute_log2_spectrum(log2_nu);  // Below cooling frequency, simple scaling
    } else {
        return log2_I_nu_peak + compute_log2_spectrum(log2_nu) + fast_log2(1 + e->Y_c) -
               fast_log2(1 + InverseComptonY::compute_Y_tilt_at_nu(e->Ys, std::exp2(log2_nu), e->p));
    }
}

void SynPhotons::update_constant() {
    // Update constants based on spectral parameters
    Real p = e->p;
    if (e->regime == 1) {
        // a_m_1_3 = std::cbrt(nu_a / nu_m);  // (nu_a / nu_m)^(1/3)
        // c_m_mpa1_2 = fastPow(nu_c / nu_m, (-p + 1) / 2);  // (nu_c / nu_m)^((-p+1)/2)
        C1_ = std::cbrt(nu_a / nu_m);
        C2_ = fast_pow(nu_c / nu_m, (-p + 1) / 2);

        log2_C1_ = fast_log2(nu_a / nu_m) / 3 - 2 * fast_log2(nu_a);
        log2_C2_ = -fast_log2(nu_m) / 3;
        log2_C3_ = (p - 1) / 2 * fast_log2(nu_m);
        log2_C4_ = (p - 1) / 2 * fast_log2(nu_m / nu_c) + p / 2 * fast_log2(nu_c);
    } else if (e->regime == 2) {
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
    } else if (e->regime == 3) {
        // a_c_1_3 = std::cbrt(nu_a / nu_c);  // (nu_a / nu_c)^(1/3)
        // c_m_1_2 = std::sqrt(nu_c / nu_m);  // (nu_c / nu_m)^(1/2)
        C1_ = std::cbrt(nu_a / nu_c);
        C2_ = std::sqrt(nu_c / nu_m);

        log2_C1_ = fast_log2(nu_a / nu_c) / 3 - 2 * fast_log2(nu_a);
        log2_C2_ = -fast_log2(nu_c) / 3;
        log2_C3_ = fast_log2(nu_c) / 2;
        log2_C4_ = fast_log2(nu_c / nu_m) / 2 + p / 2 * fast_log2(nu_m);
    } else if (e->regime == 4) {
        // a_m_1_2 = std::sqrt(nu_a / nu_m);  // (nu_a / nu_m)^(1/2)
        // R4 = std::sqrt(nu_c / nu_a) / 3;   // (nu_c / nu_a)^(1/2) / 3; // R4: scaling factor for regime 4
        C1_ = std::sqrt(nu_a / nu_m);
        C2_ = std::sqrt(nu_c / nu_a) / 3;

        log2_C1_ = -2 * fast_log2(nu_a);
        log2_C2_ = fast_log2(C2_) + fast_log2(nu_a) / 2;
        log2_C3_ = fast_log2(C2_) + fast_log2(nu_a / nu_m) / 2 + p / 2 * fast_log2(nu_m);
    } else if (e->regime == 5 || e->regime == 6) {
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
                return fast_pow(nu / nu_m, -(p - 1) / 2);
            }

            return C2_ * fast_pow(nu / nu_c, -p / 2) * fast_exp2(-nu / nu_M);

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

            return C3_ * fast_pow(nu / nu_c, -p / 2) * fast_exp2(-nu / nu_M);

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
            return C2_ * fast_pow(nu / nu_m, -p / 2) * fast_exp2(-nu / nu_M);

            break;
        case 4:
            if (nu <= nu_a) {
                return (nu / nu_a) * (nu / nu_a);
            }
            if (nu <= nu_m) {
                return C2_ * std::sqrt(nu_a / nu);
            }
            return C2_ * C1_ * fast_pow(nu / nu_m, -p / 2) * fast_exp2(-nu / nu_M);

            break;
        case 5:
            if (nu <= nu_a) {
                return (nu / nu_a) * (nu / nu_a);
            }
            return (p - 1) * C2_ * fast_pow(nu / nu_a, -p / 2) * fast_exp2(-nu / nu_M);

            break;
        case 6:
            if (nu <= nu_a) {
                return (nu / nu_a) * (nu / nu_a);
            }
            return C2_ * fast_pow(nu / nu_a, -p / 2) * fast_exp2(-nu / nu_M);

            break;

        default:
            return 0;
            break;
    }
}

Real SynPhotons::compute_log2_spectrum(Real log2_nu) const {
    Real p = e->p;
    switch (e->regime) {
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
            return log2_C4_ - p / 2 * log2_nu - fast_exp2(log2_nu) / nu_M;

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

            return log2_C4_ - p / 2 * log2_nu - fast_exp2(log2_nu) / nu_M;

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

            return log2_C4_ - p / 2 * log2_nu - fast_exp2(log2_nu) / nu_M;

            break;
        case 4:
            if (log2_nu <= log2_nu_a) {
                return log2_C1_ + 2 * log2_nu;
            }
            if (log2_nu <= log2_nu_m) {
                return log2_C2_ - log2_nu / 2;
            }

            return log2_C3_ - p / 2 * log2_nu - fast_exp2(log2_nu) / nu_M;

            break;
        case 5:
            if (log2_nu <= log2_nu_a) {
                return log2_C1_ + 2 * log2_nu;
            }

            return log2_C2_ - p / 2 * log2_nu - fast_exp2(log2_nu) / nu_M;

            break;
        case 6:
            if (log2_nu <= log2_nu_a) {
                return log2_C1_ + 2 * log2_nu;
            }

            return log2_C3_ - p / 2 * log2_nu - fast_exp2(log2_nu) / nu_M;

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

    for (; std::fabs((Y1 - Y0) / Y0) > 1e-5;) {
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

Real compute_gamma_c(Real t_com, Real B, InverseComptonY const& Ys, Real p) {
    // t_com = (6*pi*gamma*me*c^2) /(gamma^2*beta^2*sigma_T*c*B^2*(1 + Y_tilt))
    // Real gamma_c = 6 * con::pi * con::me * con::c / (con::sigmaT * B * B * (1 + Y_tilt) * t_com);

    Real Y0 = InverseComptonY::compute_Y_Thompson(Ys);
    Real gamma_bar = (6 * con::pi * con::me * con::c / con::sigmaT) / (B * B * (1 + Y0) * t_com);
    Real gamma_c = (gamma_bar + std::sqrt(gamma_bar * gamma_bar + 4)) / 2;
    Real Y1 = InverseComptonY::compute_Y_tilt_at_gamma(Ys, gamma_c, p);

    for (; std::fabs((Y1 - Y0) / Y0) > 1e-3;) {
        gamma_bar = (6 * con::pi * con::me * con::c / con::sigmaT) / (B * B * (1 + Y1) * t_com);
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
 * @return Self-absorption Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real compute_syn_gamma_a(Real Gamma_rel, Real B, Real I_syn_peak, Real gamma_m, Real gamma_c) {
    Real gamma_peak = std::min(gamma_m, gamma_c);
    Real nu_peak = compute_syn_freq(gamma_peak, B);
    Real ad_idx = adiabatic_idx(Gamma_rel);

    Real kT = (gamma_peak - 1) * (con::me * con::c2) * (ad_idx - 1);
    // 2kT(nu_a/c)^2 = I_peak*(nu_a/nu_peak)^(1/3)
    Real nu_a = fast_pow(I_syn_peak * con::c2 / (std::cbrt(nu_peak) * 2 * kT), 0.6);

    // nu_peak is not real peak, peak at nu_a; kT = (gamma_a-1) * me *c^2*(ad_idx-1), I_syn = I_peak;
    // strong absorption
    if (nu_a > nu_peak) {
        nu_a = fast_pow(
            I_syn_peak / (2 * con::me * (ad_idx - 1) * std::sqrt((4 * con::pi * con::me * con::c / (3 * con::e)) / B)),
            0.4);
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

/**
 * <!-- ************************************************************************************** -->
 * @brief Determines the peak Lorentz factor directly from a SynElectrons object.
 * @details Convenient wrapper around the three-parameter version.
 * @param e Synchrotron electron object
 * @return Peak Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real compute_gamma_peak(SynElectrons const& e) { return compute_gamma_peak(e.gamma_a, e.gamma_m, e.gamma_c); }

void update_electrons_4Y(SynElectronGrid& e, Shock const& shock) {
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
                Real p = e(i, j, k).p;
                auto& Ys = e(i, j, k).Ys;
                auto& electron = e(i, j, k);

                electron.gamma_M = compute_syn_gamma_M(B, Ys, p);  // Update maximum electron Lorentz factor
                if (k <= k_inj) {
                    electron.gamma_c = compute_gamma_c(t_com, B, Ys, p);  // Update cooling electron Lorentz factor
                } else {  // no shocked electron injection, just adiabatic cooling
                    electron.gamma_c = e(i, j, k_inj).gamma_c * electron.gamma_m / e(i, j, k_inj).gamma_m;
                    electron.gamma_M = electron.gamma_c;
                }
                electron.gamma_a =
                    compute_syn_gamma_a(Gamma_rel, B, electron.I_nu_peak, electron.gamma_m, electron.gamma_c);
                electron.regime = determine_regime(electron.gamma_a, electron.gamma_c, electron.gamma_m);
                electron.Y_c = InverseComptonY::compute_Y_tilt_at_gamma(Ys, electron.gamma_c, p);
            }
        }
    }
}

SynElectronGrid generate_syn_electrons(Shock const& shock, Real p, Real xi) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    SynElectronGrid electrons({phi_size, theta_size, t_size});

    constexpr Real gamma_cyclotron = 3;

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
                Real Sigma = shock.column_num_den(i, j, k);

                auto& e = electrons(i, j, k);

                e.gamma_M = compute_syn_gamma_M(B, electrons(i, j, k).Ys, p);
                e.gamma_m = compute_syn_gamma_m(Gamma_rel, e.gamma_M, shock.eps_e, p, xi);
                // Fraction of synchrotron electrons; the rest are cyclotron
                Real f = 1.;
                if (1 < e.gamma_m && e.gamma_m < gamma_cyclotron) {
                    f = std::min(fast_pow((gamma_cyclotron - 1) / (e.gamma_m - 1), 1 - p), 1_r);
                    e.gamma_m = gamma_cyclotron;
                }
                e.column_num_den = Sigma * f * xi;
                e.I_nu_peak = compute_syn_peak_power(B, p) * e.column_num_den / (4 * con::pi);
                if (k <= k_inj) {
                    e.gamma_c = compute_gamma_c(t_com, B, electrons(i, j, k).Ys, p);
                } else {  // no shocked electron injection, just adiabatic cooling
                    e.gamma_c = electrons(i, j, k_inj).gamma_c * e.gamma_m / electrons(i, j, k_inj).gamma_m;
                    e.gamma_M = e.gamma_c;
                }
                e.gamma_a = compute_syn_gamma_a(Gamma_rel, B, e.I_nu_peak, e.gamma_m, e.gamma_c);
                e.regime = determine_regime(e.gamma_a, e.gamma_c, e.gamma_m);
                e.p = p;
            }
        }
    }
    return electrons;
}

void generate_syn_electrons(SynElectronGrid& electrons, Shock const& shock, Real p, Real xi) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    electrons.resize({phi_size, theta_size, t_size});
    constexpr Real gamma_syn_limit = 3;

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
                Real Sigma = shock.column_num_den(i, j, k);

                auto& e = electrons(i, j, k);

                e.gamma_M = compute_syn_gamma_M(B, electrons(i, j, k).Ys, p);
                e.gamma_m = compute_syn_gamma_m(Gamma_rel, e.gamma_M, shock.eps_e, p, xi);
                // Fraction of synchrotron electrons; the rest are cyclotron
                Real f = 1.;
                if (1 < e.gamma_m && e.gamma_m < gamma_syn_limit) {
                    f = std::min(fast_pow((gamma_syn_limit - 1) / (e.gamma_m - 1), 1 - p), 1_r);
                    e.gamma_m = gamma_syn_limit;
                }
                e.column_num_den = Sigma * f * xi;
                e.I_nu_peak = compute_syn_peak_power(B, p) * e.column_num_den / (4 * con::pi);
                if (k <= k_inj) {
                    e.gamma_c = compute_gamma_c(t_com, B, electrons(i, j, k).Ys, p);
                } else {  // no shocked electron injection, just adiabatic cooling
                    e.gamma_c = electrons(i, j, k_inj).gamma_c * e.gamma_m / electrons(i, j, k_inj).gamma_m;
                    e.gamma_M = e.gamma_c;
                }
                e.gamma_a = compute_syn_gamma_a(Gamma_rel, B, e.I_nu_peak, e.gamma_m, e.gamma_c);
                e.regime = determine_regime(e.gamma_a, e.gamma_c, e.gamma_m);
                e.p = p;
            }
        }
    }
}

SynPhotonGrid generate_syn_photons(Shock const& shock, SynElectronGrid const& e) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    SynPhotonGrid ph({phi_size, theta_size, t_size});

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < t_size; ++k) {
                ph(i, j, k).e = &(e(i, j, k));
                if (shock.required(i, j, k) == 0) {
                    continue;
                }
                Real B = shock.B(i, j, k);

                ph(i, j, k).nu_M = compute_syn_freq(e(i, j, k).gamma_M, B);
                ph(i, j, k).nu_m = compute_syn_freq(e(i, j, k).gamma_m, B);
                ph(i, j, k).nu_c = compute_syn_freq(e(i, j, k).gamma_c, B);
                ph(i, j, k).nu_a = compute_syn_freq(e(i, j, k).gamma_a, B);

                ph(i, j, k).log2_I_nu_peak = fast_log2(e(i, j, k).I_nu_peak);
                ph(i, j, k).log2_nu_m = fast_log2(ph(i, j, k).nu_m);
                ph(i, j, k).log2_nu_c = fast_log2(ph(i, j, k).nu_c);
                ph(i, j, k).log2_nu_a = fast_log2(ph(i, j, k).nu_a);
                ph(i, j, k).log2_nu_M = fast_log2(ph(i, j, k).nu_M);
                ph(i, j, k).update_constant();
            }
        }
    }
    return ph;
}

void generate_syn_photons(SynPhotonGrid& ph, Shock const& shock, SynElectronGrid const& e) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    ph.resize({phi_size, theta_size, t_size});
    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < t_size; ++k) {
                ph(i, j, k).e = &(e(i, j, k));
                if (shock.required(i, j, k) == 0) {
                    continue;
                }
                Real B = shock.B(i, j, k);

                ph(i, j, k).nu_M = compute_syn_freq(e(i, j, k).gamma_M, B);
                ph(i, j, k).nu_m = compute_syn_freq(e(i, j, k).gamma_m, B);
                ph(i, j, k).nu_c = compute_syn_freq(e(i, j, k).gamma_c, B);
                ph(i, j, k).nu_a = compute_syn_freq(e(i, j, k).gamma_a, B);

                ph(i, j, k).log2_I_nu_peak = fast_log2(e(i, j, k).I_nu_peak);
                ph(i, j, k).log2_nu_m = fast_log2(ph(i, j, k).nu_m);
                ph(i, j, k).log2_nu_c = fast_log2(ph(i, j, k).nu_c);
                ph(i, j, k).log2_nu_a = fast_log2(ph(i, j, k).nu_a);
                ph(i, j, k).log2_nu_M = fast_log2(ph(i, j, k).nu_M);
                ph(i, j, k).update_constant();
            }
        }
    }
}
