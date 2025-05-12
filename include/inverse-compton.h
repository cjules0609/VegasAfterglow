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

/**
 * <!-- ************************************************************************************** -->
 * @defgroup IC_Calculation Inverse Compton Calculation Constants and Functions
 * @brief Constants and inline functions for inverse Compton photon calculations
 * <!-- ************************************************************************************** -->
 */

/// A constant used in the IC photon calculation, defined as √2/3.
inline const Real IC_x0 = std::sqrt(2) / 3;

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the Compton scattering cross-section as a function of frequency (nu).
 * @param nu The frequency at which to compute the cross-section
 * @return For x = (h * nu) / (me * c²) <= 1, returns the Thomson cross-section (con::sigmaT), otherwise 0
 * <!-- ************************************************************************************** -->
 */
inline Real compton_sigma(Real nu);

/**
 * <!-- ************************************************************************************** -->
 * @struct IntegratorGrid
 * @brief Defines a grid for numerical integration in log-space.
 * @details Given minimum and maximum values for x and y, it computes logarithmically spaced bins (x_bin and y_bin)
 *          and then determines center values (x and y) from those bins.
 * <!-- ************************************************************************************** -->
 */
struct IntegratorGrid {
    /**
     * @brief Constructor: Initializes the grid with given x and y boundaries.
     * @param x_min Minimum x-value
     * @param x_max Maximum x-value
     * @param y_min Minimum y-value
     * @param y_max Maximum y-value
     */
    IntegratorGrid(Real x_min, Real x_max, Real y_min, Real y_max)
        : x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max) {
        logspace(std::log10(x_min), std::log10(x_max), x_bin);  // Generate logarithmically spaced bin edges for x.
        logspace(std::log10(y_min), std::log10(y_max), y_bin);  // Generate logarithmically spaced bin edges for y.
        boundary_to_center(x_bin, x);                           // Compute center values for x.
        boundary_to_center(y_bin, y);                           // Compute center values for y.
    }

    Real x_min;                                        ///< Minimum x-value.
    Real x_max;                                        ///< Maximum x-value.
    Real y_min;                                        ///< Minimum y-value.
    Real y_max;                                        ///< Maximum y-value.
    static constexpr size_t num{32};                   ///< Number of bins.
    std::array<Real, num + 1> x_bin{0};                ///< Bin edges for x.
    std::array<Real, num + 1> y_bin{0};                ///< Bin edges for y.
    std::array<Real, num> x{0};                        ///< Center values for x.
    std::array<Real, num> y{0};                        ///< Center values for y.
    std::array<Real, num> j_syn{0};                    ///< Synchrotron intensity at each x center.
    std::array<Real, num> ns{0};                       ///< Number density at each y center.
    std::array<std::array<Real, num>, num> I0{{{0}}};  ///< 2D array to store computed intermediate values.
};

/**
 * <!-- ************************************************************************************** -->
 * @struct ICPhoton
 * @brief Represents a single inverse Compton (IC) photon.
 * @details Contains methods to compute the photon intensity I_nu and to generate an IC photon spectrum based
 *          on electron and synchrotron photon properties.
 * <!-- ************************************************************************************** -->
 */
struct ICPhoton {
   public:
    /// Default constructor
    ICPhoton() = default;

    /// Resolution of the computed IC spectrum.
    static constexpr size_t spectrum_resol{50};

    /**
     * <!-- ************************************************************************************** -->
     * @brief Returns the photon intensity at frequency nu.
     * @param nu The frequency at which to compute the intensity
     * @return The intensity at the given frequency
     * <!-- ************************************************************************************** -->
     */
    Real compute_I_nu(Real nu) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the base-2 logarithm of the synchrotron intensity at a given frequency.
     * @param log2_nu The base-2 logarithm of the frequency
     * @return The base-2 logarithm of the synchrotron intensity at the given frequency
     * <!-- ************************************************************************************** -->
     */
    Real compute_log2_I_nu(Real log2_nu) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Generates the IC photon spectrum from the given electron and photon data.
     * @tparam Electrons Type of the electron distribution
     * @tparam Photons Type of the photon distribution
     * @param e The electron distribution
     * @param ph The photon distribution
     * @details This template member function uses the properties of the electrons (e) and synchrotron photons (ph) to:
     *          - Determine minimum electron Lorentz factor and minimum synchrotron frequency.
     *          - Define integration limits for the synchrotron frequency (nu0) and electron Lorentz factor (gamma).
     *          - Fill in an IntegratorGrid with computed synchrotron intensity and electron column density.
     *          - Compute a 2D array I0 representing differential contributions.
     *          - Finally, integrate over the grid to populate the IC photon spectrum (j_nu_).
     * <!-- ************************************************************************************** -->
     */
    template <typename Electrons, typename Photons>
    void gen(Electrons const& e, Photons const& ph, bool KN = true) {
        // Real gamma_min = min(e.gamma_m, e.gamma_c, e.gamma_a);
        Real nu_ph_min = min(ph.nu_m, ph.nu_c, ph.nu_a);

        Real nu0_max = ph.nu_M * 5;
        Real nu0_min = nu_ph_min / 1e5;

        Real gamma_min = std::min(std::min(e.gamma_m, e.gamma_c), e.gamma_a);
        Real gamma_max = e.gamma_M * 5;

        // Construct an integration grid in nu0 and gamma.
        IntegratorGrid grid(nu0_min, nu0_max, gamma_min, gamma_max);

        // For each bin in nu0, compute the synchrotron intensity.
        for (size_t i = 0; i < grid.num; i++) {
            grid.j_syn[i] = ph.compute_I_nu(grid.x[i]);
            grid.ns[i] = e.compute_column_num_den(grid.y[i]);
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
                if (KN) {
                    grid.I0[i][j] = grid.ns[j] * grid.j_syn[i] * dS / f * compton_sigma(nu0_ / gamma_);
                } else {
                    grid.I0[i][j] = grid.ns[j] * grid.j_syn[i] * dS / f * con::sigmaT;
                }
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
                    if (nu_IC_[k] <= 4 * gamma_ * gamma_ * nu0_ * IC_x0) {
                        j_nu_[k] += grid.I0[i][j] * nu_IC_[k];
                    }
                }
            }
        }
        // Compute the base-2 logarithm of the IC photon spectrum intensity and frequency grid.
        log2_j_nu_ = xt::log2(j_nu_);
        log2_nu_IC_ = xt::log2(nu_IC_);
    };

   private:
    Array j_nu_;        ///< IC photon spectrum intensity array.
    Array nu_IC_;       ///< Frequency grid for the IC photon spectrum.
    Array log2_j_nu_;   ///< Base-2 logarithm of the IC photon spectrum intensity array.
    Array log2_nu_IC_;  ///< Base-2 logarithm of the frequency grid for the IC photon spectrum.
};

/// @typedef ICPhotonGrid
/// @brief Defines a 3D grid (using xt::xtensor) for storing ICPhoton objects.
using ICPhotonGrid = xt::xtensor<ICPhoton, 3>;

/**
 * <!-- ************************************************************************************** -->
 * @defgroup IC_Functions IC Photon and Electron Cooling Functions
 * @brief Functions to create and generate IC photon grids, and apply electron cooling mechanisms.
 * <!-- ************************************************************************************** -->
 */

/**
 * <!-- ************************************************************************************** -->
 * @brief Creates and generates an IC photon grid from electron and photon distributions
 * @param electron The electron grid
 * @param photon The photon grid
 * @return A 3D grid of IC photons
 * <!-- ************************************************************************************** -->
 */
ICPhotonGrid gen_IC_photons(SynElectronGrid const& electron, SynPhotonGrid const& photon, bool KN = true);

/**
 * <!-- ************************************************************************************** -->
 * @brief Applies Thomson cooling to electrons based on photon distribution
 * @param electron The electron grid to be modified
 * @param photon The photon grid
 * @param shock The shock properties
 * <!-- ************************************************************************************** -->
 */
void Thomson_cooling(SynElectronGrid& electron, SynPhotonGrid& photon, Shock const& shock);

/**
 * <!-- ************************************************************************************** -->
 * @brief Applies Klein-Nishina cooling to electrons based on photon distribution
 * @param electron The electron grid to be modified
 * @param photon The photon grid
 * @param shock The shock properties
 * <!-- ************************************************************************************** -->
 */
void KN_cooling(SynElectronGrid& electron, SynPhotonGrid& photon, Shock const& shock);
