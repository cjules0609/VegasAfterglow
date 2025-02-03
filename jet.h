//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _JET_
#define _JET_

#include "macros.h"
#include "mesh.h"
#include "utilities.h"

namespace func {
    inline constexpr auto zero = [](double phi, double theta, double t) -> double { return 0; };
    inline constexpr auto one = [](double phi, double theta, double t) -> double { return 1; };
}  // namespace func
class Ejecta {
   public:
    TernaryFunc dEdOmega{func::zero};
    TernaryFunc Gamma0{func::one};
    TernaryFunc dLdOmega{func::zero};
    TernaryFunc sigma0{func::zero};
    double duration{0};
};

class TophatJet {
   public:
    TophatJet(double theta_c, double E_iso, double Gamma0)
        : theta_c_(theta_c), dEdOmega_(E_iso / (4 * con::pi)), Gamma0_(Gamma0) {}
    double dEdOmega(double phi, double theta, double t) const;
    double Gamma0(double phi, double theta, double t) const;
    double sigma0(double phi, double theta, double t) const { return 0; };

    double duration{0.02 * con::sec};

   private:
    double theta_c_{0};
    double dEdOmega_{0};
    double Gamma0_{1};
};
class GaussianJet {
   public:
    GaussianJet(double theta_c, double E_iso, double Gamma0, double idx = 1)
        : theta_c_(theta_c), dEdOmega_(E_iso / (4 * con::pi)), Gamma0_(Gamma0), idx_(idx) {}
    double dEdOmega(double phi, double theta, double t) const;
    double Gamma0(double phi, double theta, double t) const;
    double sigma0(double phi, double theta, double t) const { return 0; };

    double duration{0.02 * con::sec};

   private:
    double theta_c_{0};
    double dEdOmega_{0};
    double Gamma0_{1};
    double idx_{0};
};

class PowerLawJet {
   public:
    PowerLawJet(double theta_c, double E_iso, double Gamma0, double k, double idx = 1)
        : theta_c_(theta_c), dEdOmega_(E_iso / (4 * con::pi)), Gamma0_(Gamma0), k_(k), idx_(idx) {}
    double dEdOmega(double phi, double theta, double t) const;
    double Gamma0(double phi, double theta, double t) const;
    double sigma0(double phi, double theta, double t) const { return 0; };

    double duration{0.02 * con::sec};

   private:
    double theta_c_{0};
    double dEdOmega_{0};
    double Gamma0_{1};
    double k_{0};
    double idx_{0};
};

namespace inject {
    inline struct {
        inline double dEdOmega(double phi, double theta, double t) const { return 0; };
        inline double Gamma0(double phi, double theta, double t) const { return 1; };
        inline double dLdOmega(double phi, double theta, double t) const { return 0; };
        inline double sigma0(double phi, double theta, double t) const { return 0; };
    } none;
}  // namespace inject

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
            return height * fastExp(-theta * theta / (2 * theta_c * theta_c));
        };
    }

    inline auto powerLaw(double theta_c, double height, double k) {
        return [=](double phi, double theta, double t = 0) -> double {
            return theta < theta_c ? height : height * fastPow(theta / theta_c, -k);
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
        return [=](double t) -> double { return fastPow(1 + t / t0, -q); };
    }

    inline auto powerLawIntegral(double t0, double q) {
        return [=](double t) -> double {
            if (std::fabs(q - 1) > 1e-6) {
                return t0 / (1 - q) * (fastPow(1 + t / t0, 1 - q) - 1);
            } else {
                return t0 * fastLog(1 + t / t0);
            }
        };
    }

}  // namespace math

auto LiangGhirlanda2010(auto energy_func, double e_max, double gamma_max, double idx);
#endif