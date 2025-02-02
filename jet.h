#ifndef _JET_
#define _JET_

#include "macros.h"
#include "mesh.h"

inline auto return_zero = [](double phi, double theta, double t) -> double { return 0; };

inline auto return_one = [](double phi, double theta, double t) -> double { return 1; };

class Ejecta {
   public:
    TernaryFunc dEdOmega{return_zero};
    TernaryFunc Gamma0{return_one};
    TernaryFunc dLdOmega{return_zero};
    TernaryFunc sigma0{return_zero};
    double duration{0};
};

namespace math {
    inline auto combine(auto f_spatial, auto f_temporal) {
        return [=](double phi, double theta, double t) -> double { return f_spatial(phi, theta) * f_temporal(t); };
    }

    inline auto isotropic(double height) {
        return [=](double phi, double theta, double t = 0) -> double { return height; };
    }

    inline auto tophat(double theta_c, double hight) {
        return [=](double phi, double theta, double t = 0) -> double { return theta < theta_c ? hight : 0; };
    }

    inline auto gaussian(double theta_c, double height) {
        return [=](double phi, double theta, double t = 0) -> double {
            return height * std::exp(-theta * theta / (2 * theta_c * theta_c));
        };
    }

    inline auto powerLaw(double theta_c, double height, double k) {
        return [=](double phi, double theta, double t = 0) -> double {
            return theta < theta_c ? height : height * std::pow(theta / theta_c, -k);
        };
    }

    inline auto constInjection() {
        return [=](double t) -> double { return 1; };
    }

    inline auto constIntegral() {
        return [=](double t) -> double { return t; };
    }

    inline auto stepInjection(double t0) {
        return [=](double t) -> double { return t > t0 ? 1 : 0; };
    }

    inline auto stepIntegral(double t0) {
        return [=](double t) -> double { return t > t0 ? t - t0 : 0; };
    }

    inline auto squareInjection(double t0, double t1) {
        return [=](double t) -> double { return t > t0 && t < t1 ? 1 : 0; };
    }

    inline auto squareIntegral(double t0, double t1) {
        return [=](double t) -> double {
            if (t < t0) {
                return 0;
            } else if (t < t1) {
                return t - t0;
            } else {
                return t1 - t0;
            }
        };
    }

    inline auto powerLawInjection(double t0, double q) {
        return [=](double t) -> double { return std::pow(1 + t / t0, -q); };
    }

    inline auto powerLawIntegral(double t0, double q) {
        return [=](double t) -> double {
            if (std::fabs(q - 1) > 1e-6) {
                return t0 / (1 - q) * (std::pow(1 + t / t0, 1 - q) - 1);
            } else {
                return t0 * std::log(1 + t / t0);
            }
        };
    }

}  // namespace math

auto LiangGhirlanda2010(auto energy_func, double e_max, double gamma_max, double idx);

Ejecta tophatJet(double theta_c, double E_iso, double Gamma0);
Ejecta gaussianJet(double theta_c, double E_iso, double Gamma0, double idx = 1);
Ejecta powerLawJet(double theta_c, double E_iso, double Gamma0, double k, double idx = 1);
std::tuple<double, double> findRadiusRange(double t_min, double t_max, double z, Ejecta const& jet);
#endif