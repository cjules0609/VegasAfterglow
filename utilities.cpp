#include "utilities.h"

#include <cmath>
#include <iostream>

#include "macros.h"

double point_interp(double x0, double x1, double y0, double y1, double xi) {
    if (x0 == x1) return y0;
    double slope = (y1 - y0) / (x1 - x0);
    return y0 + slope * (xi - x0);
}

double point_loglog_interp(double x0, double x1, double y0, double y1, double xi) {
    if (y0 == 0 || y1 == 0) return 0;
    if (x0 == x1) return y0;
    double log_x0 = log(x0);
    double log_x1 = log(x1);
    double log_y0 = log(y0);
    double log_y1 = log(y1);
    double slope = (log_y1 - log_y0) / (log_x1 - log_x0);
    return exp(log_y0 + slope * (log(xi) - log_x0));
}

double interp(double xi, Array const& x, Array const& y, bool lo_extrap, bool hi_extrap) {
    if (x.size() < 2 || y.size() < 2 || x.size() != y.size()) {
        std::cout << "incorrect array size for interpolation!\n";
        return 0;
    }
    auto x_back = x[x.size() - 1];
    auto y_back = y[y.size() - 1];

    if (xi < x[0]) {
        return (!lo_extrap || x[0] == xi) ? y[0] : point_interp(x[0], x[1], y[0], y[1], xi);
    } else if (xi > x_back) {
        return (!hi_extrap || x_back == xi) ? y_back
                                            : point_interp(x[x.size() - 2], x_back, y[y.size() - 2], y_back, xi);
    } else {
        auto it = std::lower_bound(x.begin(), x.end(), xi);
        size_t idx = it - x.begin();
        if (*it == xi) return y[idx];  // Exact match
        return point_interp(x[idx - 1], x[idx], y[idx - 1], y[idx], xi);
    }
}

double interp_eq_spaced(double xi, Array const& x, Array const& y, bool lo_extrap, bool hi_extrap) {
    if (x.size() < 2 || y.size() < 2 || x.size() != y.size()) {
        std::cout << "incorrect array size for interpolation!\n";
        return 0;
    }

    auto x_back = x[x.size() - 1];
    auto y_back = y[y.size() - 1];

    if (xi <= x[0])
        return (!lo_extrap || x[0] == xi) ? y[0] : point_interp(x[0], x[1], y[0], y[1], xi);
    else if (xi >= x_back)
        return (!hi_extrap || x_back == xi) ? y_back
                                            : point_interp(x[x.size() - 2], x_back, y[y.size() - 2], y_back, xi);
    else {
        double dx = x[1] - x[0];
        size_t idx = static_cast<size_t>((xi - x[0]) / dx + 1);
        if (xi == x[idx]) return y[idx];
        return point_interp(x[idx - 1], x[idx], y[idx - 1], y[idx], xi);
    }
}

double loglog_interp(double xi, const Array& x, const Array& y, bool lo_extrap, bool hi_extrap) {
    if (x.size() < 2 || y.size() < 2 || x.size() != y.size()) {
        std::cout << "incorrect array size for interpolation!\n";
        return 0;
    }
    auto x_back = x[x.size() - 1];
    auto y_back = y[y.size() - 1];

    if (xi <= x[0]) {
        return (!lo_extrap || x[0] == xi) ? y[0] : point_loglog_interp(x[0], x[1], y[0], y[1], xi);
    } else if (xi >= x_back) {
        return (!hi_extrap || x_back == xi) ? y_back
                                            : point_loglog_interp(x[x.size() - 2], x_back, y[y.size() - 2], y_back, xi);
    } else {
        auto it = std::lower_bound(x.begin(), x.end(), xi);
        size_t idx = it - x.begin();
        if (*it == xi) return y[idx];  // Exact match
        return point_loglog_interp(x[idx - 1], x[idx], y[idx - 1], y[idx], xi);
    }
}

double loglog_interp_eq_spaced(double xi, const Array& x, const Array& y, bool lo_extrap, bool hi_extrap) {
    if (x.size() < 2 || y.size() < 2 || x.size() != y.size()) {
        std::cout << "incorrect array size for interpolation!\n";
        return 0;
    }
    auto x_back = x[x.size() - 1];
    auto y_back = y[y.size() - 1];

    if (xi <= x[0]) {
        // std::cout << "here!" << (!lo_extrap || x[0] == xi) ? y[0] : point_loglog_interp(x[0], x[1], y[0], y[1], xi);
        return (!lo_extrap || x[0] == xi) ? y[0] : point_loglog_interp(x[0], x[1], y[0], y[1], xi);
    } else if (xi >= x_back) {
        return (!hi_extrap || x_back == xi) ? y_back
                                            : point_loglog_interp(x[x.size() - 2], x_back, y[y.size() - 2], y_back, xi);
    } else {
        double log_x0 = log(x[0]);
        double dx = log(x[1]) - log_x0;
        size_t idx = static_cast<size_t>((log(xi) - log_x0) / dx + 1);

        if (xi == x[idx]) return y[idx];  // Exact match
        return point_loglog_interp(x[idx - 1], x[idx], y[idx - 1], y[idx], xi);
    }
}

double fast_exp(double x) {
    /*constexpr double a = (1ll << 52) / 0.6931471805599453;
    constexpr double b = (1ll << 52) * (1023 - 0.04367744890362246);
    x = a * x + b;

    // Remove these lines if bounds checking is not needed
    // constexpr double c = (1ll << 52);
    // constexpr double d = (1ll << 52) * 2047;
    // if (x < c || x > d)
    //    x = (x < c) ? 0.0 : d;

    // With C++20 one can use std::bit_cast instead
    uint64_t n = static_cast<uint64_t>(x);
    memcpy(&x, &n, 8);

    return x;*/
    return std::exp(x);
}

/*inline double exp_fast(double a) {
    union {
        double d;
        long long x;
    } u;
    u.x = (long long)(6497320848556798LL * a + 0x3fef127e83d16f12LL);
    return u.d;
}*/

double fast_pow(double a, double b) {
    uint64_t bits = std::bit_cast<uint64_t>(a);

    uint64_t exponent = static_cast<uint64_t>(b * ((static_cast<int64_t>((bits >> 52) & 0x7FF) - 1023)) + 1023);

    bits = (exponent << 52);

    return std::bit_cast<double>(bits);

    // return std::pow(a, b);
}
/*
inline double fast_pow(double a, double b) {
    union {
        double d;
        int x[2];
    } u = {a};
    u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
    u.x[0] = 0;
    return u.d;
    // return pow(a, b);
}*/

double jet_edge(Profile const& gamma) {
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

Array adaptive_theta_space(size_t n, Profile const& gamma) {
    if (n == 1) {
        return linspace(0, 2 * con::pi, 2);
        ;
    }
    double edge = jet_edge(gamma);

    return adaptive_theta_space(n, gamma, edge);
}

Array adaptive_theta_space(size_t n, Profile const& gamma, double edge) {
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