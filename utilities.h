#ifndef _UTILIT_
#define _UTILIT_
#include "macros.h"
#include "mesh.h"
template <typename Fun>
auto root_bisection(Fun f, decltype(f(0)) low, decltype(f(0)) high, decltype(f(0)) eps = 1e-6) -> decltype(f(0)) {
    using Scalar = decltype(f(0));

    for (; fabs((high - low)) > fabs(high) * eps;) {
        Scalar mid = 0.5 * (high + low);
        if (f(mid) > 0)
            high = mid;
        else
            low = mid;
    }
    return 0.5 * (high + low);
}

double exp_fast(double a);
double logscale_interp_extrap(double x0, Array const& x, Array const& y);
double logscale_interp_extrap_eq_spaced(double x0, Array const& x, Array const& y);
Array adaptive_theta_space(size_t n, Profile const& gamma);
double interp(double x0, Array const& x, Array const& y);
double interp_extrap(double x0, Array const& x, Array const& y);
double interp_extrap_eq_spaced(double x0, Array const& x, Array const& y);
inline double step_func(double x) { return x > 0 ? 1 : 0; }
inline double eVtoHz(double eV) { return eV / con::h; };
#endif