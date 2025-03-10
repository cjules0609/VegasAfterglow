//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _UTILITIES_H_
#define _UTILITIES_H_
#include "macros.h"
#include "mesh.h"
/********************************************************************************************************************
 * FUNCTION: Basic Math Functions                                                                                   *
 * DESCRIPTION: Inline functions for specific power calculations, a step function, and unit conversion.             *
 ********************************************************************************************************************/
inline Real pow52(Real a) { return std::sqrt(a * a * a * a * a); }
inline Real pow43(Real a) { return std::cbrt(a * a * a * a); }
inline Real pow23(Real a) { return std::cbrt(a * a); }
inline Real stepFunc(Real x) { return x > 0 ? 1 : 0; }
inline Real eVtoHz(Real eV) { return eV / con::h; }

/********************************************************************************************************************
 * FUNCTION: Fast Math & Interpolation Prototypes                                                                   *
 * DESCRIPTION: Prototypes for fast power, logarithm, and various interpolation functions.                          *
 ********************************************************************************************************************/
Real fastPow(Real a, Real b);
Real fastLog(Real x);
Real interp(Real x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);
Real interpEqSpaced(Real x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);
Real loglogInterp(Real x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);
Real loglogInterpEqSpaced(Real x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);

/********************************************************************************************************************
 * FUNCTION: Root Finding (Bisection Method)                                                                        *
 * DESCRIPTION: Template function to find the root of a function using the bisection method.                        *
 ********************************************************************************************************************/
template <typename Fun>
auto rootBisection(Fun f, decltype(f(0)) low, decltype(f(0)) high, decltype(f(0)) eps = 1e-6) -> decltype(f(0)) {
    using Scalar = decltype(f(0));
    for (; (high - low) > std::abs((high + low) * 0.5) * eps;) {
        Scalar mid = 0.5 * (high + low);
        if (f(mid) * f(high) > 0)
            high = mid;
        else
            low = mid;
    }
    return 0.5 * (high + low);
}

/********************************************************************************************************************
 * FUNCTION: Utility Templates                                                                                      *
 * DESCRIPTION: Template functions for computing the minimum and maximum of provided values.                        *
 ********************************************************************************************************************/
template <typename T>
T min(T value) {
    return value;
}

template <typename T, typename... Args>
T min(T first, Args... args) {
    return std::min(first, std::min(args...));
}

template <typename T>
T max(T value) {
    return value;
}

template <typename T, typename... Args>
T max(T first, Args... args) {
    return std::max(first, std::max(args...));
}

/********************************************************************************************************************
 * FUNCTION: Fast Exponential and Logarithm Functions                                                               *
 * DESCRIPTION: Inline functions that provide fast approximations of exponential and logarithm functions using      *
 *              alternative methods when EXTREME_SPEED is defined.                                                  *
 ********************************************************************************************************************/
inline Real fastExp(Real x) {
#ifdef EXTREME_SPEED
    // if (std::isnan(x)) return std::numeric_limits<Real>::quiet_NaN();
    // if (x == std::numeric_limits<Real>::infinity()) return std::numeric_limits<Real>::infinity();
    // if (x == -std::numeric_limits<Real>::infinity()) return 0.0;

    constexpr Real ln2 = 0.6931471805599453;
    constexpr Real inv_ln2 = 1.4426950408889634;

    Real y = x * inv_ln2;
    int64_t k = static_cast<int64_t>(y + (y >= 0 ? 0.5 : -0.5));
    Real r = x - k * ln2;

    // Real p = 1.0 + r * (1.0 + r * (0.5 + r * (0.166666666666666 + r * 0.041666666666666664)));

    Real p = 1.0 + r * (1.0 + r * (0.5 + r * (0.166666666666666)));

    return std::ldexp(p, k);
#else
    return std::exp(x);
#endif
}

inline double fastLog(double x) {
#ifdef EXTREME_SPEED
    if (x <= 0.) return -std::numeric_limits<double>::infinity();
    if (std::isnan(x)) return std::numeric_limits<double>::quiet_NaN();
    if (x == std::numeric_limits<double>::infinity()) return std::numeric_limits<double>::infinity();

    uint64_t bits;
    std::memcpy(&bits, &x, sizeof(x));
    int64_t exponent = ((bits >> 52) & 0x7FF) - 1023;
    bits = (bits & 0x000FFFFFFFFFFFFFULL) | 0x3FF0000000000000ULL;
    double f;
    std::memcpy(&f, &bits, sizeof(f));
    double p = -1.49278 + (2.11263 + (-0.729104 + 0.10969 * f) * f) * f;
    return p + 0.6931471806 * exponent;
#else
    return std::log(x);
#endif
}

inline Real fastPow(Real a, Real b) { return fastExp(b * fastLog(a)); }

#endif  // _UTILITIES_H_
