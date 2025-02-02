#include "utilities.h"

#include <cmath>
#include <iostream>

#include "macros.h"

double pointInterp(double x0, double x1, double y0, double y1, double xi) {
    if (x0 == x1) return y0;
    double slope = (y1 - y0) / (x1 - x0);
    return y0 + slope * (xi - x0);
}

double pointLoglogInterp(double x0, double x1, double y0, double y1, double xi) {
    if (y0 == 0 || y1 == 0) return 0;
    if (x0 == x1) return y0;
    double log_x0 = std::log(x0);
    double log_x1 = std::log(x1);
    double log_y0 = std::log(y0);
    double log_y1 = std::log(y1);
    double slope = (log_y1 - log_y0) / (log_x1 - log_x0);
    return std::exp(log_y0 + slope * (std::log(xi) - log_x0));
}

double interp(double xi, Array const& x, Array const& y, bool lo_extrap, bool hi_extrap) {
    if (x.size() < 2 || y.size() < 2 || x.size() != y.size()) {
        std::cout << "incorrect array size for interpolation!\n";
        return 0;
    }
    auto x_back = x[x.size() - 1];
    auto y_back = y[y.size() - 1];

    if (xi < x[0]) {
        return (!lo_extrap || x[0] == xi) ? y[0] : pointInterp(x[0], x[1], y[0], y[1], xi);
    } else if (xi > x_back) {
        return (!hi_extrap || x_back == xi) ? y_back
                                            : pointInterp(x[x.size() - 2], x_back, y[y.size() - 2], y_back, xi);
    } else {
        auto it = std::lower_bound(x.begin(), x.end(), xi);
        size_t idx = it - x.begin();
        if (*it == xi) return y[idx];  // Exact match
        return pointInterp(x[idx - 1], x[idx], y[idx - 1], y[idx], xi);
    }
}

double interpEqSpaced(double xi, Array const& x, Array const& y, bool lo_extrap, bool hi_extrap) {
    if (x.size() < 2 || y.size() < 2 || x.size() != y.size()) {
        std::cout << "incorrect array size for interpolation!\n";
        return 0;
    }

    auto x_back = x[x.size() - 1];
    auto y_back = y[y.size() - 1];

    if (xi <= x[0])
        return (!lo_extrap || x[0] == xi) ? y[0] : pointInterp(x[0], x[1], y[0], y[1], xi);
    else if (xi >= x_back)
        return (!hi_extrap || x_back == xi) ? y_back
                                            : pointInterp(x[x.size() - 2], x_back, y[y.size() - 2], y_back, xi);
    else {
        double dx = x[1] - x[0];
        size_t idx = static_cast<size_t>((xi - x[0]) / dx + 1);
        if (xi == x[idx]) return y[idx];
        return pointInterp(x[idx - 1], x[idx], y[idx - 1], y[idx], xi);
    }
}

double loglogInterp(double xi, const Array& x, const Array& y, bool lo_extrap, bool hi_extrap) {
    if (x.size() < 2 || y.size() < 2 || x.size() != y.size()) {
        std::cout << "incorrect array size for interpolation!\n";
        return 0;
    }
    auto x_back = x[x.size() - 1];
    auto y_back = y[y.size() - 1];

    if (xi <= x[0]) {
        return (!lo_extrap || x[0] == xi) ? y[0] : pointLoglogInterp(x[0], x[1], y[0], y[1], xi);
    } else if (xi >= x_back) {
        return (!hi_extrap || x_back == xi) ? y_back
                                            : pointLoglogInterp(x[x.size() - 2], x_back, y[y.size() - 2], y_back, xi);
    } else {
        auto it = std::lower_bound(x.begin(), x.end(), xi);
        size_t idx = it - x.begin();
        if (*it == xi) return y[idx];  // Exact match
        return pointLoglogInterp(x[idx - 1], x[idx], y[idx - 1], y[idx], xi);
    }
}

double loglogInterpEqSpaced(double xi, const Array& x, const Array& y, bool lo_extrap, bool hi_extrap) {
    if (x.size() < 2 || y.size() < 2 || x.size() != y.size()) {
        std::cout << "incorrect array size for interpolation!\n";
        return 0;
    }
    auto x_back = x[x.size() - 1];
    auto y_back = y[y.size() - 1];

    if (xi <= x[0]) {
        // std::cout << "here!" << (!lo_extrap || x[0] == xi) ? y[0] : point_loglog_interp(x[0], x[1], y[0], y[1], xi);
        return (!lo_extrap || x[0] == xi) ? y[0] : pointLoglogInterp(x[0], x[1], y[0], y[1], xi);
    } else if (xi >= x_back) {
        return (!hi_extrap || x_back == xi) ? y_back
                                            : pointLoglogInterp(x[x.size() - 2], x_back, y[y.size() - 2], y_back, xi);
    } else {
        double log_x0 = std::log(x[0]);
        double dx = std::log(x[1]) - log_x0;
        size_t idx = static_cast<size_t>((std::log(xi) - log_x0) / dx + 1);

        if (xi == x[idx]) return y[idx];  // Exact match
        return pointLoglogInterp(x[idx - 1], x[idx], y[idx - 1], y[idx], xi);
    }
}

double fastLog(double a) {
#ifdef EXTREME_SPEED
    /*uint64_t bx = std::bit_cast<uint64_t>(a);
    uint64_t ex = bx >> 52;
    int64_t t = static_cast<int64_t>(ex) - 1023;
    bx = 0x3FF0000000000000ULL | (bx & 0x000FFFFFFFFFFFFFULL);
    double x = std::bit_cast<double>(bx);*/
    // if (a <= 0) return std::numeric_limits<double>::quiet_NaN();

    uint64_t bx = *(uint64_t*)&a;
    uint64_t ex = bx >> 52;
    int64_t t = static_cast<int64_t>(ex) - 1023;
    bx = 0x3FF0000000000000ULL | (bx & 0x000FFFFFFFFFFFFFULL);
    double x = *(double*)&bx;  // C-style cast back to double

    return -1.49278 + (2.11263 + (-0.729104 + 0.10969 * x) * x) * x + 0.6931471806 * t;
#else
    return std::log(a);
#endif
}

double fastExp(double a) {
#ifdef EXTREME_SPEED
    static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__, "Little endian required!");
    // if (a < -708.39641853226408) return 0;
    // if (a > 709.782712893384) return std::numeric_limits<double>::infinity();
    union {
        double d;
        long long x;
    } u;
    u.x = (long long)(6497320848556798LL * a + 0x3fef127e83d16f12LL);
    return u.d;
#else
    return std::exp(a);
#endif
}

double fastPow(double a, double b) {
#ifdef EXTREME_SPEED
    static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__, "Little endian required!");
    union {
        double d;
        int x[2];
    } u = {a};
    u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
    u.x[0] = 0;
    return u.d;
#else
    return std::pow(a, b);
#endif
}

double jetEdge(UnaryFunc const& gamma) {
    if (gamma(con::pi / 2) > 1) {
        return con::pi / 2;
    }
    double low = 0;
    double hi = con::pi / 2;
    double eps = 1e-6;
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

Array adaptiveThetaSpace(size_t n, UnaryFunc const& gamma) {
    if (n == 1) {
        return linspace(0, 2 * con::pi, 2);
        ;
    }
    double edge = jetEdge(gamma);

    return adaptiveThetaSpace(n, gamma, edge);
}

Array adaptiveThetaSpace(size_t n, UnaryFunc const& gamma, double edge) {
    if (n == 1) {
        return linspace(0, 2 * con::pi, 2);
        ;
    }

    double dx = con::pi / (n - 1);
    Array space = zeros(n);
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