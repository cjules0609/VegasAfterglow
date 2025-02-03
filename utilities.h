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
inline double pow52(double a) { return std::sqrt(a * a * a * a * a); }
inline double pow43(double a) { return std::cbrt(a * a * a * a); }
inline double pow23(double a) { return std::cbrt(a * a); }
inline double stepFunc(double x) { return x > 0 ? 1 : 0; }
inline double eVtoHz(double eV) { return eV / con::h; }

double fastPow(double a, double b);
double fastLog(double x);
double interp(double x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);
double interpEqSpaced(double x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);
double loglogInterp(double x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);
double loglogInterpEqSpaced(double x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);
// Function to find the root of a function using the bisection method
template <typename Fun>
auto rootBisection(Fun f, decltype(f(0)) low, decltype(f(0)) high, decltype(f(0)) eps = 1e-6) -> decltype(f(0)) {
    using Scalar = decltype(f(0));
    for (; std::fabs((high - low)) > std::fabs(high) * eps;) {
        Scalar mid = 0.5 * (high + low);
        if (f(mid) * f(high) > 0)
            high = mid;
        else
            low = mid;
    }
    return 0.5 * (high + low);
}

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

inline double fastExp(double x) {
#ifdef EXTREME_SPEED
    if (std::isnan(x)) return std::numeric_limits<double>::quiet_NaN();
    if (x == std::numeric_limits<double>::infinity()) return std::numeric_limits<double>::infinity();
    if (x == -std::numeric_limits<double>::infinity()) return 0.0;

    constexpr double ln2 = 0.6931471805599453;
    constexpr double inv_ln2 = 1.4426950408889634;

    double y = x * inv_ln2;
    int64_t k = static_cast<int64_t>(y + (y >= 0 ? 0.5 : -0.5));
    double r = x - k * ln2;

    // double p = 1.0 + r * (1.0 + r * (0.5 + r * (0.166666666666666 + r * 0.041666666666666664)));

    double p = 1.0 + r * (1.0 + r * (0.5 + r * (0.166666666666666)));

    return std::ldexp(p, k);
#else
    return std::exp(x);
#endif
}

inline double fastLog(double x) {
#ifdef EXTREME_SPEED
    if (x <= 0.0) return -std::numeric_limits<double>::infinity();
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

inline double fastPow(double a, double b) { return fastExp(b * fastLog(a)); }

#endif  // _UTILITIES_H_
