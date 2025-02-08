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
    inline constexpr auto zero = [](Real phi, Real theta, Real t) -> Real { return 0; };
    // Always returns 1 regardless of the input.
    inline constexpr auto one = [](Real phi, Real theta, Real t) -> Real { return 1; };
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
    Real duration{0};                  // Duration of the ejecta

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
    TophatJet(Real theta_c, Real E_iso, Real Gamma0)
        : theta_c_(theta_c), dEdOmega_(E_iso / (4 * con::pi)), Gamma0_(Gamma0) {}
    // Returns the energy per unit solid angle at (phi, theta, t). (For a tophat jet, this is constant within the core.)
    Real dEdOmega(Real phi, Real theta, Real t) const;
    // Returns the initial Lorentz factor at (phi, theta, t). (Constant for a tophat jet.)
    Real Gamma0(Real phi, Real theta, Real t) const;
    // Returns the magnetization parameter, which is zero for a tophat jet.
    Real sigma0(Real phi, Real theta, Real t) const { return 0; };

    Real duration{0.02 * con::sec};  // Jet duration

   private:
    Real theta_c_{0};   // Core opening angle (radians)
    Real dEdOmega_{0};  // Energy per unit solid angle (computed from E_iso)
    Real Gamma0_{1};    // Initial Lorentz factor
};

/********************************************************************************************************************
 * CLASS: GaussianJet
 * DESCRIPTION: Represents a Gaussian jet model. The jet's properties vary with angle according to a Gaussian profile.
 ********************************************************************************************************************/
class GaussianJet {
   public:
    // Constructor: Sets the core angle (theta_c), isotropic equivalent energy (E_iso), initial Lorentz factor (Gamma0),
    // and the Gaussian index (idx).
    GaussianJet(Real theta_c, Real E_iso, Real Gamma0, Real idx = 1)
        : two_theta_c_sq(2 * theta_c * theta_c), dEdOmega_(E_iso / (4 * con::pi)), Gamma0_(Gamma0), idx_(idx) {}
    Real dEdOmega(Real phi, Real theta, Real t) const;
    Real Gamma0(Real phi, Real theta, Real t) const;
    Real sigma0(Real phi, Real theta, Real t) const { return 0; };

    Real duration{0.02 * con::sec};

   private:
    Real two_theta_c_sq{0};  // Core opening angle
    Real dEdOmega_{0};       // Energy per unit solid angle (computed from E_iso)
    Real Gamma0_{1};         // Initial Lorentz factor
    Real idx_{0};            // Gaussian index parameter
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
    PowerLawJet(Real theta_c, Real E_iso, Real Gamma0, Real k, Real idx = 1)
        : theta_c_(theta_c), dEdOmega_(E_iso / (4 * con::pi)), Gamma0_(Gamma0), k_(k), idx_(idx) {}
    Real dEdOmega(Real phi, Real theta, Real t) const;
    Real Gamma0(Real phi, Real theta, Real t) const;
    Real sigma0(Real phi, Real theta, Real t) const { return 0; };

    Real duration{0.02 * con::sec};

   private:
    Real theta_c_{0};   // Core opening angle
    Real dEdOmega_{0};  // Energy per unit solid angle
    Real Gamma0_{1};    // Initial Lorentz factor
    Real k_{0};         // Power-law index for the angular profile
    Real idx_{0};       // Additional index parameter
};

/********************************************************************************************************************
 * NAMESPACE: inject
 * DESCRIPTION: Contains injector models for energy and Lorentz factor injection.
 ********************************************************************************************************************/
namespace inject {
    // A dummy injector that returns zero for energy and luminosity injection and one for the Lorentz factor.
    inline struct {
        inline Real dEdOmega(Real phi, Real theta, Real t) const { return 0; };
        inline Real Gamma0(Real phi, Real theta, Real t) const { return 1; };
        inline Real dLdOmega(Real phi, Real theta, Real t) const { return 0; };
        inline Real sigma0(Real phi, Real theta, Real t) const { return 0; };
    } none;
}  // namespace inject

/********************************************************************************************************************
 * NAMESPACE: math
 * DESCRIPTION: Provides a collection of mathematical helper functions for combining functions,
 *              constructing various injection profiles (e.g., tophat, Gaussian, power-law), and computing integrals.
 ********************************************************************************************************************/
namespace math {
    // Combines a spatial function and a temporal function into one function of (phi, theta, t).
    template <typename F1, typename F2>
    inline auto combine(F1 f_spatial, F2 f_temporal) {
        return [=](Real phi, Real theta, Real t) -> Real { return f_spatial(phi, theta) * f_temporal(t); };
    }

    // Returns a constant (isotropic) function with a fixed height.
    inline auto isotropic(Real height) {
        return [=](Real phi, Real theta, Real t = 0) -> Real { return height; };
    }

    // Returns a tophat function: it is constant (height) within the core angle theta_c and zero outside.
    inline auto tophat(Real theta_c, Real hight) {
        return [=](Real phi, Real theta, Real t = 0) -> Real { return theta < theta_c ? hight : 0; };
    }

    // Returns a Gaussian profile function for jet properties.
    inline auto gaussian(Real theta_c, Real height) {
        return [=](Real phi, Real theta, Real t = 0) -> Real {
            return height * fastExp(-theta * theta / (2 * theta_c * theta_c));
        };
    }

    // Returns a power-law profile function for jet properties.
    inline auto powerLaw(Real theta_c, Real height, Real k) {
        return [=](Real phi, Real theta, Real t = 0) -> Real {
            return theta < theta_c ? height : height * fastPow(theta / theta_c, -k);
        };
    }

    // Constant injection: returns 1 regardless of time.
    inline auto constInjection() {
        return [=](Real t) -> Real { return 1; };
    }

    // Constant integral: returns t (linearly increasing with time).
    inline auto constIntegral() {
        return [=](Real t) -> Real { return t; };
    }

    // Step injection: returns 1 if t is greater than t0, else 0.
    inline auto stepInjection(Real t0) {
        return [=](Real t) -> Real { return t > t0 ? 1 : 0; };
    }

    // Step integral: returns t-t0 if t is greater than t0, else 0.
    inline auto stepIntegral(Real t0) {
        return [=](Real t) -> Real { return t > t0 ? t - t0 : 0; };
    }

    // Square injection: returns 1 if t is between t0 and t1, else 0.
    inline auto squareInjection(Real t0, Real t1) {
        return [=](Real t) -> Real { return t > t0 && t < t1 ? 1 : 0; };
    }

    // Square integral: integrates the square injection function over time.
    inline auto squareIntegral(Real t0, Real t1) {
        return [=](Real t) -> Real {
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
    inline auto powerLawInjection(Real t0, Real q) {
        return [=](Real t) -> Real { return fastPow(1 + t / t0, -q); };
    }

    // Power-law integral: integrates the power-law injection over time, handling the q = 1 case separately.
    inline auto powerLawIntegral(Real t0, Real q) {
        return [=](Real t) -> Real {
            if (std::abs(q - 1) > 1e-6) {
                return t0 / (1 - q) * (fastPow(1 + t / t0, 1 - q) - 1);
            } else {
                return t0 * fastLog(1 + t / t0);
            }
        };
    }
}  // namespace math

/********************************************************************************************************************
 * FUNCTION: LiangGhirlanda2010
 * DESCRIPTION: Returns a lambda function that computes a Lorentz factor using the Liang & Ghirlanda (2010)
 *              prescription. The lambda takes an energy function (energy_func) and parameters e_max, gamma_max,
 *              and idx, and returns a function of (phi, theta, t). The returned function calculates:
 *                  e = energy_func(phi, theta, 0)
 *                  u = (e / e_max)^idx * gamma_max
 *                  return sqrt(1 + u^2)
 ********************************************************************************************************************/
template <typename F>
auto LiangGhirlanda2010(F energy_func, Real e_max, Real gamma_max, Real idx) {
    return [=](Real phi, Real theta, Real t = 0) -> Real {
        Real e = energy_func(phi, theta, 0);
        Real u = fastPow(e / e_max, idx) * gamma_max;
        return std::sqrt(1 + u * u);
    };
}
#endif