#ifndef _UTILITIES_H_
#define _UTILITIES_H_
#include "macros.h"
#include "mesh.h"

// Function to find the root of a function using the bisection method
template <typename Fun>
auto root_bisection(Fun f, decltype(f(0)) low, decltype(f(0)) high, decltype(f(0)) eps = 1e-6) -> decltype(f(0)) {
    using Scalar = decltype(f(0));
    for (; fabs((high - low)) > fabs(high) * eps;) {
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
    return std::min(first, min(args...));
}

template <typename T>
T max(T value) {
    return value;
}

template <typename T, typename... Args>
T max(T first, Args... args) {
    return std::max(first, max(args...));
}

double exp_fast(double a);

// Linear interpolation function
double interp(double x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);

double interp_eq_spaced(double x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);

double loglog_interp(double x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);

double loglog_interp_eq_spaced(double x0, Array const& x, Array const& y, bool lo_extrap = false,
                               bool hi_extrap = false);

inline double step_func(double x) { return x > 0 ? 1 : 0; }

inline double eVtoHz(double eV) { return eV / con::h; }

Array adaptive_theta_space(size_t n, Profile const& gamma);

Array adaptive_theta_space(size_t n, Profile const& gamma, double edge);

#endif  // _UTILITIES_H_
