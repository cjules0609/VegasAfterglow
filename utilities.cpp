#include "utilities.h"

#include <cmath>
#include <iostream>

#include "macros.h"
double interp(double x0, Array const& x, Array const& y) {
    if (x0 < x[0]) {
        return y[0];
    } else if (x0 > x.back()) {
        return y.back();
    } else {
        auto it = std::lower_bound(x.begin(), x.end(), x0);
        size_t i = it - x.begin();
        if (x0 == x[i]) {
            return y[i];
        } else {
            double a = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
            return y[i - 1] + a * (x0 - x[i - 1]);
        }
    }
}

double logscale_interp_extrap(double xi, Array const& x, Array const& y) {
    double logxi = log(xi);
    if (xi < x[0]) {
        if (y[0] == 0 || y[1] == 0) {
            return 0;
        }
        double y1 = log(y[1]);
        double y0 = log(y[0]);
        double x1 = log(x[1]);
        double x0 = log(x[0]);
        double a = (y1 - y0) / (x1 - x0);
        return exp(y0 + a * (logxi - x0));
    } else if (xi > x.back()) {
        if (y[y.size() - 1] == 0 || y[y.size() - 2] == 0) {
            return 0;
        }
        double y1 = log(y[y.size() - 2]);
        double y0 = log(y[y.size() - 1]);
        double x1 = log(x[x.size() - 2]);
        double x0 = log(x[x.size() - 1]);
        double a = (y1 - y0) / (x1 - x0);
        return exp(y0 + a * (logxi - x0));
    } else {
        auto it = std::lower_bound(x.begin(), x.end(), xi);
        size_t i = it - x.begin();
        if (xi == x[i]) {
            return y[i];
        } else {
            if (y[i] == 0 || y[i - 1] == 0) {
                return 0;
            }
            double y1 = log(y[i]);
            double y0 = log(y[i - 1]);
            double x1 = log(x[i]);
            double x0 = log(x[i - 1]);
            double a = (y1 - y0) / (x1 - x0);
            return exp(y0 + a * (logxi - x0));
        }
    }
}

double interp_extrap(double x0, Array const& x, Array const& y) {
    if (x0 < x[0]) {
        return (y[1] - y[0]) / (x[1] - x[0]) * (x0 - x[0]) + y[0];
    } else if (x0 > x.back()) {
        size_t i = y.size() - 2;
        size_t j = y.size() - 1;
        return (y[j] - y[i]) / (x[j] - x[i]) * (x0 - x[j]) + y[j];
    } else {
        auto it = std::lower_bound(x.begin(), x.end(), x0);
        size_t i = it - x.begin();
        if (x0 == x[i]) {
            return y[i];
        } else {
            double a = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
            return y[i - 1] + a * (x0 - x[i - 1]);
        }
    }
}

double logscale_interp_extrap_eq_spaced(double xi, Array const& x, Array const& y) {
    double logxi = log(xi);
    if (xi < x[0]) {
        if (y[0] == 0 || y[1] == 0) {
            return 0;
        }
        double y1 = log(y[1]);
        double y0 = log(y[0]);
        double x1 = log(x[1]);
        double x0 = log(x[0]);
        double a = (y1 - y0) / (x1 - x0);
        return exp(y0 + a * (logxi - x0));
    } else if (xi > x.back()) {
        if (y[y.size() - 1] == 0 || y[y.size() - 2] == 0) {
            return 0;
        }
        double y1 = log(y[y.size() - 2]);
        double y0 = log(y[y.size() - 1]);
        double x1 = log(x[x.size() - 2]);
        double x0 = log(x[x.size() - 1]);
        double a = (y1 - y0) / (x1 - x0);
        return exp(y0 + a * (logxi - x0));
    } else {
        double x0_ = log(x[0]);
        double x1_ = log(x[1]);
        double dx = x1_ - x0_;
        size_t i = (logxi - x0_) / dx + 1;
        if (xi == x[i]) {
            return y[i];
        } else {
            if (y[i] == 0 || y[i - 1] == 0) {
                return 0;
            }
            double y1 = log(y[i]);
            double y0 = log(y[i - 1]);
            double x0 = log(x[i - 1]);
            double a = (y1 - y0) / dx;
            return exp(y0 + a * (logxi - x0));
        }
    }
}

double interp_extrap_eq_spaced(double x0, Array const& x, Array const& y) {
    if (x0 < x[0]) {
        return (y[1] - y[0]) / (x[1] - x[0]) * (x0 - x[0]) + y[0];
    } else if (x0 > x.back()) {
        size_t i = y.size() - 2;
        size_t j = y.size() - 1;
        return (y[j] - y[i]) / (x[j] - x[i]) * (x0 - x[j]) + y[j];
    } else {
        size_t i = (x0 - x[0]) / (x[1] - x[0]) + 1;
        if (x0 == x[i]) {
            return y[i];
        } else {
            double a = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
            return y[i - 1] + a * (x0 - x[i - 1]);
        }
    }
}

double exp_fast(double a) {
    union {
        double d;
        long long x;
    } u;
    u.x = (long long)(6497320848556798LL * a + 0x3fef127e83d16f12LL);
    return u.d;
}

double jet_edge(Profile const& gamma) {
    if (gamma(con::pi / 2) > 1) {
        return con::pi / 2;
    }
    double low = 0;
    double hi = con::pi / 2;
    double eps = 1e-9;
    for (; hi - low > eps;) {
        double mid = 0.5 * (low + hi);
        if (gamma(mid) > 1) {
            low = mid;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (low + hi);
}

Array adaptive_theta_space(size_t n, Profile const& gamma) {
    if (n == 1) {
        return Array(0, con::pi / 2);
    }
    double edge = jet_edge(gamma);

    double dx = con::pi / (n - 1);
    Array space(n, 0);
    for (size_t i = 0; i < n - 1; ++i) {
        double x = i * dx;
        space[i + 1] = space[i] + 1 / gamma(x);
    }
    double rescale = edge / (space[n - 1]);
    for (size_t i = 0; i < n; ++i) {
        space[i] = (space[i]) * rescale;
    }
    return space;
}
