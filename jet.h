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

/********************************************************************************************************************
 * NAMESPACE: func
 * DESCRIPTION: Contains inline constexpr lambda functions that return constant values.
 ********************************************************************************************************************/
namespace func {
    // Always returns 0 regardless of the input.
    inline constexpr auto zero = [](double phi, double theta, double t) -> double { return 0; };
    // Always returns 1 regardless of the input.
    inline constexpr auto one = [](double phi, double theta, double t) -> double { return 1; };
}  // namespace func

/********************************************************************************************************************
 * CLASS: Ejecta
 * DESCRIPTION: Represents ejecta properties for a simulation. It uses ternary functions (TernaryFunc) to
 *              describe various quantities (dEdOmega, Gamma0, dLdOmega, sigma0) as functions of phi, theta, and time.
 ********************************************************************************************************************/
class Ejecta {
   public:
    TernaryFunc dEdOmega{func::zero};  // Differential energy per unit solid angle (default: zero)
    TernaryFunc Gamma0{func::one};     // Initial Lorentz factor (default: one)
    TernaryFunc dLdOmega{func::zero};  // Differential luminosity per unit solid angle (default: zero)
    TernaryFunc sigma0{func::zero};    // Magnetization parameter (default: zero)
    double duration{0};                // Duration of the ejecta

    // (Additional member functions could be defined as needed.)
};

/********************************************************************************************************************
 * CLASS: TophatJet
 * DESCRIPTION: Represents a tophat jet model. The jet has a uniform structure within a cone defined by theta_c.
 *              The energy per unit solid angle and the initial Lorentz factor are constant inside the cone.
 ********************************************************************************************************************/
class TophatJet {
   public:
    // Constructor: Sets the core angle (theta_c), isotropic equivalent energy (E_iso), and initial Lorentz factor
    // (Gamma0).
    TophatJet(double theta_c, double E_iso, double Gamma0)
        : theta_c_(theta_c), dEdOmega_(E_iso / (4 * con::pi)), Gamma0_(Gamma0) {}
    // Returns the energy per unit solid angle at (phi, theta, t). (For a tophat jet, this is constant within the core.)
    double dEdOmega(double phi, double theta, double t) const;
    // Returns the initial Lorentz factor at (phi, theta, t). (Constant for a tophat jet.)
    double Gamma0(double phi, double theta, double t) const;
    // Returns the magnetization parameter, which is zero for a tophat jet.
    double sigma0(double phi, double theta, double t) const { return 0; };

    double duration{0.02 * con::sec};  // Jet duration

   private:
    double theta_c_{0};   // Core opening angle (radians)
    double dEdOmega_{0};  // Energy per unit solid angle (computed from E_iso)
    double Gamma0_{1};    // Initial Lorentz factor
};

/********************************************************************************************************************
 * CLASS: GaussianJet
 * DESCRIPTION: Represents a Gaussian jet model. The jet's properties vary with angle according to a Gaussian profile.
 ********************************************************************************************************************/
class GaussianJet {
   public:
    // Constructor: Sets the core angle (theta_c), isotropic equivalent energy (E_iso), initial Lorentz factor (Gamma0),
    // and the Gaussian index (idx).
    GaussianJet(double theta_c, double E_iso, double Gamma0, double idx = 1)
        : theta_c_(theta_c), dEdOmega_(E_iso / (4 * con::pi)), Gamma0_(Gamma0), idx_(idx) {}
    double dEdOmega(double phi, double theta, double t) const;
    double Gamma0(double phi, double theta, double t) const;
    double sigma0(double phi, double theta, double t) const { return 0; };

    double duration{0.02 * con::sec};

   private:
    double theta_c_{0};   // Core opening angle
    double dEdOmega_{0};  // Energy per unit solid angle (computed from E_iso)
    double Gamma0_{1};    // Initial Lorentz factor
    double idx_{0};       // Gaussian index parameter
};

/********************************************************************************************************************
 * CLASS: PowerLawJet
 * DESCRIPTION: Represents a power-law jet model. Outside the core angle theta_c, the jet properties (e.g., energy per
 *              unit solid angle) fall off as a power law with index k.
 ********************************************************************************************************************/
class PowerLawJet {
   public:
    // Constructor: Sets the core angle (theta_c), isotropic equivalent energy (E_iso), initial Lorentz factor (Gamma0),
    // power-law index (k), and an optional index parameter (idx).
    PowerLawJet(double theta_c, double E_iso, double Gamma0, double k, double idx = 1)
        : theta_c_(theta_c), dEdOmega_(E_iso / (4 * con::pi)), Gamma0_(Gamma0), k_(k), idx_(idx) {}
    double dEdOmega(double phi, double theta, double t) const;
    double Gamma0(double phi, double theta, double t) const;
    double sigma0(double phi, double theta, double t) const { return 0; };

    double duration{0.02 * con::sec};

   private:
    double theta_c_{0};   // Core opening angle
    double dEdOmega_{0};  // Energy per unit solid angle
    double Gamma0_{1};    // Initial Lorentz factor
    double k_{0};         // Power-law index for the angular profile
    double idx_{0};       // Additional index parameter
};

/********************************************************************************************************************
 * NAMESPACE: inject
 * DESCRIPTION: Contains injector models for energy and Lorentz factor injection.
 ********************************************************************************************************************/
namespace inject {
    // A dummy injector that returns zero for energy and luminosity injection and one for the Lorentz factor.
    inline struct {
        inline double dEdOmega(double phi, double theta, double t) const { return 0; };
        inline double Gamma0(double phi, double theta, double t) const { return 1; };
        inline double dLdOmega(double phi, double theta, double t) const { return 0; };
        inline double sigma0(double phi, double theta, double t) const { return 0; };
    } none;
}  // namespace inject

/********************************************************************************************************************
 * NAMESPACE: math
 * DESCRIPTION: Provides a collection of mathematical helper functions for combining functions,
 *              constructing various injection profiles (e.g., tophat, Gaussian, power-law), and computing integrals.
 ********************************************************************************************************************/
namespace math {
    // Combines a spatial function and a temporal function into one function of (phi, theta, t).
    inline auto combine(auto f_spatial, auto f_temporal) {
        return [=](double phi, double theta, double t) -> double { return f_spatial(phi, theta) * f_temporal(t); };
    }

    // Returns a constant (isotropic) function with a fixed height.
    inline auto isotropic(double height) {
        return [=](double phi, double theta, double t = 0) -> double { return height; };
    }

    // Returns a tophat function: it is constant (height) within the core angle theta_c and zero outside.
    inline auto tophat(double theta_c, double hight) {
        return [=](double phi, double theta, double t = 0) -> double { return theta < theta_c ? hight : 0; };
    }

    // Returns a Gaussian profile function for jet properties.
    inline auto gaussian(double theta_c, double height) {
        return [=](double phi, double theta, double t = 0) -> double {
            return height * fastExp(-theta * theta / (2 * theta_c * theta_c));
        };
    }

    // Returns a power-law profile function for jet properties.
    inline auto powerLaw(double theta_c, double height, double k) {
        return [=](double phi, double theta, double t = 0) -> double {
            return theta < theta_c ? height : height * fastPow(theta / theta_c, -k);
        };
    }

    // Constant injection: returns 1 regardless of time.
    inline auto constInjection() {
        return [=](double t) -> double { return 1; };
    }

    // Constant integral: returns t (linearly increasing with time).
    inline auto constIntegral() {
        return [=](double t) -> double { return t; };
    }

    // Step injection: returns 1 if t is greater than t0, else 0.
    inline auto stepInjection(double t0) {
        return [=](double t) -> double { return t > t0 ? 1 : 0; };
    }

    // Step integral: returns t-t0 if t is greater than t0, else 0.
    inline auto stepIntegral(double t0) {
        return [=](double t) -> double { return t > t0 ? t - t0 : 0; };
    }

    // Square injection: returns 1 if t is between t0 and t1, else 0.
    inline auto squareInjection(double t0, double t1) {
        return [=](double t) -> double { return t > t0 && t < t1 ? 1 : 0; };
    }

    // Square integral: integrates the square injection function over time.
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

    // Power-law injection: returns a decaying function with power-law index q.
    inline auto powerLawInjection(double t0, double q) {
        return [=](double t) -> double { return fastPow(1 + t / t0, -q); };
    }

    // Power-law integral: integrates the power-law injection over time, handling the q = 1 case separately.
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

/********************************************************************************************************************
 * FUNCTION: LiangGhirlanda2010
 * DESCRIPTION: Declaration of a function (or function template) that combines an energy function with specified
 *              parameters (e_max, gamma_max, idx) following the model by Liang & Ghirlanda (2010). The precise
 *              implementation details are not provided here.
 ********************************************************************************************************************
 */
auto LiangGhirlanda2010(auto energy_func, double e_max, double gamma_max, double idx);
#endif