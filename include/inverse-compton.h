//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <array>
#include <cmath>
#include <vector>

#include "macros.h"
#include "medium.h"
#include "mesh.h"
#include "shock.h"
#include "synchrotron.h"
#include "utilities.h"
/********************************************************************************************************************
 * CONSTANTS AND INLINE FUNCTIONS FOR IC PHOTON CALCULATIONS
 ********************************************************************************************************************/

// A constant used in the IC photon calculation, defined as √2/3.
inline const Real IC_x0 = std::sqrt(2) / 3;

// Computes the Compton scattering cross-section as a function of frequency (nu).
// For x = (h * nu) / (me * c²) <= 1, it returns the Thomson cross-section (con::sigmaT).
// Otherwise, it returns 0.
// (An alternative implementation is commented out.)
inline Real comptonSigma(Real nu) {
    Real x = con::h * nu / (con::me * con::c2);
    if (x <= 1) {
        return con::sigmaT;
    } else {
        return 0;
    }
    /*
    if (x < 1e-2) {
         return con::sigmaT * (1 - 2 * x);
     } else if (x > 1e2) {
         return 3. / 8 * con::sigmaT * (log(2 * x) + 1.0 / 2) / x;
     } else {
         return 0.75 * con::sigmaT *
                ((1 + x) / (x * x * x) * (2 * x * (1 + x) / (1 + 2 * x) - log(1 + 2 * x)) + log(1 + 2 * x) / (2 * x) -
                 (1 + 3 * x) / (1 + 2 * x) / (1 + 2 * x));
     }
    */
}

/********************************************************************************************************************
 * STRUCT: IntegratorGrid
 * DESCRIPTION: Defines a grid for numerical integration in log-space.
 *              Given minimum and maximum values for x and y, it computes logarithmically spaced bins (x_bin and y_bin)
 *              and then determines center values (x and y) from those bins.
 ********************************************************************************************************************/
struct IntegratorGrid {
    // Constructor: Initializes the grid with given x and y boundaries.
    IntegratorGrid(Real x_min, Real x_max, Real y_min, Real y_max)
        : x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max) {
        logspace(std::log10(x_min), std::log10(x_max), x_bin);  // Generate logarithmically spaced bin edges for x.
        logspace(std::log10(y_min), std::log10(y_max), y_bin);  // Generate logarithmically spaced bin edges for y.
        boundaryToCenter(x_bin, x);                             // Compute center values for x.
        boundaryToCenter(y_bin, y);                             // Compute center values for y.
    }

    Real x_min;                                        // Minimum x-value.
    Real x_max;                                        // Maximum x-value.
    Real y_min;                                        // Minimum y-value.
    Real y_max;                                        // Maximum y-value.
    static constexpr size_t num{80};                   // Number of bins.
    std::array<Real, num + 1> x_bin{0};                // Bin edges for x.
    std::array<Real, num + 1> y_bin{0};                // Bin edges for y.
    std::array<Real, num> x{0};                        // Center values for x.
    std::array<Real, num> y{0};                        // Center values for y.
    std::array<Real, num> j_syn{0};                    // Synchrotron intensity at each x center.
    std::array<Real, num> ns{0};                       // Number density at each y center.
    std::array<std::array<Real, num>, num> I0{{{0}}};  // 2D array to store computed intermediate values.
};

/********************************************************************************************************************
 * STRUCT: ICPhoton
 * DESCRIPTION: Represents a single inverse Compton (IC) photon.
 *              Contains methods to compute the photon intensity I_nu and to generate an IC photon spectrum based
 *              on electron and synchrotron photon properties.
 ********************************************************************************************************************/
struct ICPhoton {
   public:
    ICPhoton() = default;

    // Resolution of the computed IC spectrum.
    static constexpr size_t spectrum_resol{50};

    // Returns the photon intensity at frequency nu.
    Real I_nu(Real nu) const;

    // Generates the IC photon spectrum from the given electron and photon data.
    // This template member function uses the properties of the electrons (e) and synchrotron photons (ph) to:
    //   - Determine minimum electron Lorentz factor and minimum synchrotron frequency.
    //   - Define integration limits for the synchrotron frequency (nu0) and electron Lorentz factor (gamma).
    //   - Fill in an IntegratorGrid with computed synchrotron intensity and electron column density.
    //   - Compute a 2D array I0 representing differential contributions.
    //   - Finally, integrate over the grid to populate the IC photon spectrum (j_nu_).
    template <typename Electrons, typename Photons>
    void gen(Electrons const& e, Photons const& ph) {
        // Real gamma_e_min = min(e.gamma_m, e.gamma_c, e.gamma_a);
        Real nu_ph_min = min(ph.nu_m, ph.nu_c, ph.nu_a);

        Real nu0_max = ph.nu_M * 10;
        Real nu0_min = nu_ph_min / 1e5;

        Real gamma_min = std::min(std::min(e.gamma_m, e.gamma_c), e.gamma_a);
        Real gamma_max = e.gamma_M * 10;

        // Construct an integration grid in nu0 and gamma.
        IntegratorGrid grid(nu0_min, nu0_max, gamma_min, gamma_max);

        // For each bin in nu0, compute the synchrotron intensity.
        for (size_t i = 0; i < grid.num; i++) {
            grid.j_syn[i] = ph.I_nu(grid.x[i]);
            grid.ns[i] = e.columnNumDen(grid.y[i]);
        }

        // For each (nu0, gamma) pair, compute differential contributions and fill in I0.
        for (size_t i = 0; i < grid.num; ++i) {
            Real nu0_ = grid.x[i];
            Real dnu = grid.x_bin[i + 1] - grid.x_bin[i];
            for (size_t j = 0; j < grid.num; ++j) {
                Real gamma_ = grid.y[j];
                Real dgamma = grid.y_bin[j + 1] - grid.y_bin[j];
                Real dS = std::fabs(dnu * dgamma);
                Real f = 4 * gamma_ * gamma_ * nu0_ * nu0_;
                grid.I0[i][j] = grid.ns[j] * grid.j_syn[i] * dS / f * comptonSigma(nu0_ / gamma_);
            }
        }

        // Compute integration limits for the IC spectrum.
        Real nu_min = 4 * gamma_min * gamma_min * nu0_min;
        Real nu_max = 4 * gamma_max * gamma_max * nu0_max;

        // Generate the IC frequency grid and allocate the output array.
        nu_IC_ = xt::logspace(std::log10(nu_min), std::log10(nu_max), spectrum_resol);
        j_nu_ = xt::zeros<Real>({spectrum_resol});

        // Integrate over the grid to compute the final IC photon spectrum.
        for (size_t k = 0; k < nu_IC_.size(); ++k) {
            for (size_t i = 0; i < grid.num; ++i) {
                Real nu0_ = grid.x[i];
                for (size_t j = 0; j < grid.num; ++j) {
                    Real gamma_ = grid.y[j];
                    if (nu0_ <= nu_IC_[k] && nu_IC_[k] <= 4 * gamma_ * gamma_ * nu0_ * IC_x0) {
                        j_nu_[k] += grid.I0[i][j] * nu_IC_[k];
                    }
                }
            }
        }
    };

   private:
    Array j_nu_;   // IC photon spectrum intensity array.
    Array nu_IC_;  // Frequency grid for the IC photon spectrum.
};

/********************************************************************************************************************
 * TYPE ALIAS: ICPhotonGrid
 * DESCRIPTION: Defines a 3D grid (using xt::xtensor) for storing ICPhoton objects.
 ********************************************************************************************************************/
using ICPhotonGrid = xt::xtensor<ICPhoton, 3>;

/********************************************************************************************************************
 * FUNCTION PROTOTYPES: IC Photon and Electron Cooling Functions
 * DESCRIPTION: These functions create and generate IC photon grids, and apply electron cooling mechanisms.
 ********************************************************************************************************************/
ICPhotonGrid genICPhotons(SynElectronGrid const& electron, SynPhotonGrid const& photon);
void eCoolingThomson(SynElectronGrid& electron, SynPhotonGrid const& photon, Shock const& shock);
void eCoolingKleinNishina(SynElectronGrid& electron, SynPhotonGrid const& photon, Shock const& shock);
