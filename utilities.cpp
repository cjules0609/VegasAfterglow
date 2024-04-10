#include "utilities.h"

#include <cmath>
#include <iostream>
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

double interp_log(double xi, Array const& x, Array const& y) {
    if (xi < x[0]) {
        return y[0];
    } else if (xi > x.back()) {
        return y.back();
    } else {
        auto it = std::lower_bound(x.begin(), x.end(), xi);
        size_t i = it - x.begin();
        if (xi == x[i]) {
            return y[i];
        } else {
            if (y[i] == 0 || y[i - 1] == 0) {
                return 0;
            }
            double lgxi = log10(xi);
            double y1 = log10(y[i]);
            double y0 = log10(y[i - 1]);
            double x1 = log10(x[i]);
            double x0 = log10(x[i - 1]);
            double a = (y1 - y0) / (x1 - x0);
            return pow(10, y0 + a * (lgxi - x0));
        }
    }
}

double interp_log_extra_lo(double xi, Array const& x, Array const& y) {
    if (xi < x[0]) {
        if (y[0] == 0 || y[1] == 0) {
            return 0;
        }
        double lgxi = log10(xi);
        double y1 = log10(y[1]);
        double y0 = log10(y[0]);
        double x1 = log10(x[1]);
        double x0 = log10(x[0]);
        double a = (y1 - y0) / (x1 - x0);
        return pow(10, y0 + a * (lgxi - x0));
    } else if (xi > x.back()) {
        return y.back();
    } else {
        auto it = std::lower_bound(x.begin(), x.end(), xi);
        size_t i = it - x.begin();
        if (xi == x[i]) {
            return y[i];
        } else {
            if (y[i] == 0 || y[i - 1] == 0) {
                return 0;
            }
            double lgxi = log10(xi);
            double y1 = log10(y[i]);
            double y0 = log10(y[i - 1]);
            double x1 = log10(x[i]);
            double x0 = log10(x[i - 1]);
            double a = (y1 - y0) / (x1 - x0);
            return pow(10, y0 + a * (lgxi - x0));
        }
    }
}

double interp_log_extra_hi(double xi, Array const& x, Array const& y) {
    if (xi < x[0]) {
        return y[0];
    } else if (xi > x.back()) {
        if (y[y.size() - 1] == 0 || y[y.size() - 2] == 0) {
            return 0;
        }
        double lgxi = log10(xi);
        double y1 = log10(y[y.size() - 2]);
        double y0 = log10(y[y.size() - 1]);
        double x1 = log10(x[x.size() - 2]);
        double x0 = log10(x[x.size() - 1]);
        double a = (y1 - y0) / (x1 - x0);
        return pow(10, y0 + a * (lgxi - x0));
    } else {
        auto it = std::lower_bound(x.begin(), x.end(), xi);
        size_t i = it - x.begin();
        if (xi == x[i]) {
            return y[i];
        } else {
            if (y[i] == 0 || y[i - 1] == 0) {
                return 0;
            }
            double lgxi = log10(xi);
            double y1 = log10(y[i]);
            double y0 = log10(y[i - 1]);
            double x1 = log10(x[i]);
            double x0 = log10(x[i - 1]);
            double a = (y1 - y0) / (x1 - x0);
            return pow(10, y0 + a * (lgxi - x0));
        }
    }
}

double interp_log_extra_both(double xi, Array const& x, Array const& y) {
    double lgxi = log10(xi);
    if (xi < x[0]) {
        if (y[0] == 0 || y[1] == 0) {
            return 0;
        }
        double y1 = log10(y[1]);
        double y0 = log10(y[0]);
        double x1 = log10(x[1]);
        double x0 = log10(x[0]);
        double a = (y1 - y0) / (x1 - x0);
        return pow(10, y0 + a * (lgxi - x0));
    } else if (xi > x.back()) {
        if (y[y.size() - 1] == 0 || y[y.size() - 2] == 0) {
            return 0;
        }
        double y1 = log10(y[y.size() - 2]);
        double y0 = log10(y[y.size() - 1]);
        double x1 = log10(x[x.size() - 2]);
        double x0 = log10(x[x.size() - 1]);
        double a = (y1 - y0) / (x1 - x0);
        return pow(10, y0 + a * (lgxi - x0));
    } else {
        auto it = std::lower_bound(x.begin(), x.end(), xi);
        size_t i = it - x.begin();
        if (xi == x[i]) {
            return y[i];
        } else {
            if (y[i] == 0 || y[i - 1] == 0) {
                return 0;
            }
            double y1 = log10(y[i]);
            double y0 = log10(y[i - 1]);
            double x1 = log10(x[i]);
            double x0 = log10(x[i - 1]);
            double a = (y1 - y0) / (x1 - x0);
            return pow(10, y0 + a * (lgxi - x0));
        }
    }
}

double interp_extra_lo(double x0, Array const& x, Array const& y) {
    if (x0 < x[0]) {
        return (y[1] - y[0]) / (x[1] - x[0]) * (x0 - x[0]) + y[0];
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

double interp_extra_hi(double x0, Array const& x, Array const& y) {
    if (x0 < x[0]) {
        return y[0];
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

double interp_extra_both(double x0, Array const& x, Array const& y) {
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

double exp_fast(double a) {
    union {
        double d;
        long long x;
    } u;
    u.x = (long long)(6497320848556798LL * a + 0x3fef127e83d16f12LL);
    return u.d;
}