#ifndef _UTILIT_
#define _UTILIT_
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

double interp(double x0, Array const& x, Array const& y);
inline double stepfun(double x) { return x > 0 ? 1 : 0; }
#endif