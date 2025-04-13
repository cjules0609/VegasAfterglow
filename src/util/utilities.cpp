//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "utilities.h"

#include <cmath>
#include <iostream>
#include <numeric>

#include "macros.h"

/********************************************************************************************************************
 * FUNCTION: Point Interpolation Utility Functions
 * DESCRIPTION: Provides basic point-wise interpolation functions used by the higher-level interpolation routines.
 ********************************************************************************************************************/
Real pointInterp(Real x0, Real x1, Real y0, Real y1, Real xi) {
    if (x0 == x1) return y0;
    Real slope = (y1 - y0) / (x1 - x0);
    return y0 + slope * (xi - x0);
}

Real pointLoglogInterp(Real x0, Real x1, Real y0, Real y1, Real xi) {
    if (y0 == 0 || y1 == 0) return 0;
    if (x0 == x1) return y0;
    Real log_x0 = std::log(x0);
    Real log_x1 = std::log(x1);
    Real log_y0 = std::log(y0);
    Real log_y1 = std::log(y1);
    Real slope = (log_y1 - log_y0) / (log_x1 - log_x0);
    return std::exp(log_y0 + slope * (std::log(xi) - log_x0));
}

/********************************************************************************************************************
 * FUNCTION: Linear Interpolation Functions
 * DESCRIPTION: Implements linear interpolation (with and without the assumption of equally spaced x values).
 ********************************************************************************************************************/
Real interp(Real xi, Array const& x, Array const& y, bool lo_extrap, bool hi_extrap) {
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

Real interpEqSpaced(Real xi, Array const& x, Array const& y, bool lo_extrap, bool hi_extrap) {
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
        Real dx = x[1] - x[0];
        size_t idx = static_cast<size_t>((xi - x[0]) / dx + 1);
        if (xi == x[idx]) return y[idx];
        return pointInterp(x[idx - 1], x[idx], y[idx - 1], y[idx], xi);
    }
}

/********************************************************************************************************************
 * FUNCTION: Log-Log Interpolation Functions
 * DESCRIPTION: Implements log–log interpolation (with and without the assumption of equally spaced x values in log
 *space).
 ********************************************************************************************************************/
Real loglogInterp(Real xi, const Array& x, const Array& y, bool lo_extrap, bool hi_extrap) {
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

Real loglogInterpEqSpaced(Real xi, const Array& x, const Array& y, bool lo_extrap, bool hi_extrap) {
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
        Real log_x0 = std::log(x[0]);
        Real dx = std::log(x[1]) - log_x0;
        size_t idx = static_cast<size_t>((std::log(xi) - log_x0) / dx + 1);

        if (xi == x[idx]) return y[idx];  // Exact match
        return pointLoglogInterp(x[idx - 1], x[idx], y[idx - 1], y[idx], xi);
    }
}