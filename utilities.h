#ifndef _UTILITIES_H_
#define _UTILITIES_H_
#include "macros.h"
#include "mesh.h"

inline double pow52(double a) { return std::sqrt(a * a * a * a * a); }
inline double pow43(double a) { return std::cbrt(a * a * a * a); }
inline double pow23(double a) { return std::cbrt(a * a); }
inline double stepFunc(double x) { return x > 0 ? 1 : 0; }
inline double eVtoHz(double eV) { return eV / con::h; }

double fastExp(double x);
double fastPow(double a, double b);
double fastLog(double x);
double interp(double x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);
double interpEqSpaced(double x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);
double loglogInterp(double x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);
double loglogInterpEqSpaced(double x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);
Array adaptiveThetaSpace(size_t n, Profile const& gamma);
Array adaptiveThetaSpace(size_t n, Profile const& gamma, double edge);

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

inline double boxcar(double x, double a, double b) {
    if (a <= x && x <= b) {
        return 1;
    } else {
        return 0;
    }
}
#endif  // _UTILITIES_H_
