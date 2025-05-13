//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <array>
#include <cmath>
#include <tuple>
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
inline Real compton_cross_section(Real nu);

/**
 * <!-- ************************************************************************************** -->
 * @struct IntegratorGrid
 * @brief Defines a grid for numerical integration in log-space. Stack struct to avoid memory allocation.
 * @details Given minimum and maximum values for nu and gamma, it computes logarithmically spaced bins (nu_edge and
 * gamma_edge) and then determines center values (nu and gamma) from those bins.
 * <!-- ************************************************************************************** -->
 */
struct IntegratorGrid {
    template <std::size_t N>
    using StackArray = xt::xtensor_fixed<Real, xt::xshape<N>>;

    template <std::size_t N, std::size_t M>
    using StackMesh = xt::xtensor_fixed<Real, xt::xshape<N, M>>;
    /**
     * @brief Constructor: Initializes the grid with given nu and gamma boundaries.
     * @param nu_min Minimum nu-value
     * @param nu_max Maximum nu-value
     * @param gamma_min Minimum gamma-value
     * @param gamma_max Maximum gamma-value
     */
    IntegratorGrid(Real nu_min, Real nu_max, Real gamma_min, Real gamma_max) {
        // Generate logarithmically spaced bin edges for nu.
        nu_edge = xt::logspace(std::log10(nu_min), std::log10(nu_max), num + 1);
        // Generate logarithmically spaced bin edges for gamma.
        gamma_edge = xt::logspace(std::log10(gamma_min), std::log10(gamma_max), num + 1);
        boundary_to_center(nu_edge, nu);        // Compute center values for nu.
        boundary_to_center(gamma_edge, gamma);  // Compute center values for gamma.
    }

    static constexpr size_t num{64};    ///< Number of bins.
    StackArray<num + 1> nu_edge{0};     ///< Bin edges for nu.
    StackArray<num + 1> gamma_edge{0};  ///< Bin edges for gamma.
    StackArray<num> nu{0};              ///< Center values for nu.
    StackArray<num> gamma{0};           ///< Center values for gamma.
    StackArray<num> I_nu_syn{0};        ///< Synchrotron intensity at each nu center.
    StackArray<num> column_den{0};      ///< Number density at each gamma center.
    StackMesh<num, num> I0{{{0}}};      ///< 2D array to store computed intermediate values.
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
    static constexpr size_t spectrum_resol{64};

    /**
     * <!-- ************************************************************************************** -->
     * @brief Returns the photon intensity at frequency nu.
     * @param nu The frequency at which to compute the intensity
     * @return The intensity at the given frequency
     * <!-- ************************************************************************************** -->
     */
    Real compute_I_nu(Real nu) const noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the base-2 logarithm of the synchrotron intensity at a given frequency.
     * @param log2_nu The base-2 logarithm of the frequency
     * @return The base-2 logarithm of the synchrotron intensity at the given frequency
     * <!-- ************************************************************************************** -->
     */
    Real compute_log2_I_nu(Real log2_nu) const noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Generates the IC photon spectrum from the given electron and photon data.
     * @tparam Electrons Type of the electron distribution
     * @tparam Photons Type of the photon distribution
     * @param electrons The electron distribution
     * @param photons The photon distribution
     * @param KN Whether to use the Klein-Nishina corrected cross-section
     * @details This template member function uses the properties of the electrons (e) and synchrotron photons (ph) to:
     *          - Determine minimum electron Lorentz factor and minimum synchrotron frequency.
     *          - Define integration limits for the synchrotron frequency (nu0) and electron Lorentz factor (gamma).
     *          - Fill in an IntegratorGrid with computed synchrotron intensity and electron column density.
     *          - Compute a 2D array I0 representing differential contributions.
     *          - Finally, integrate over the grid to populate the IC photon spectrum (I_nu_IC_).
     * <!-- ************************************************************************************** -->
     */
    template <typename Electrons, typename Photons>
    void compute_IC_spectrum(Electrons const& electrons, Photons const& photons, bool KN = true) noexcept;

   public:
    Array I_nu_IC_;     ///< IC photon spectrum intensity array.
    Array nu_IC_;       ///< Frequency grid for the IC photon spectrum.
    Array log2_I_nu_;   ///< Base-2 logarithm of the IC photon spectrum intensity array.
    Array log2_nu_IC_;  ///< Base-2 logarithm of the frequency grid for the IC photon spectrum.

    /**
     * <!-- ************************************************************************************** -->
     * @brief Get the integration bounds for the IC photon spectrum.
     * @param electrons The electron distribution
     * @param photons The photon distribution
     * @tparam Electrons Type of the electron distribution
     * @tparam Photons Type of the photon distribution
     * @details This template member function computes the minimum and maximum values for the synchrotron frequency
     * (nu0) and electron Lorentz factor (gamma)
     * @return A tuple containing the minimum and maximum values for nu0, gamma, nu_min, and nu_max
     * <!-- ************************************************************************************** -->
     */
    template <typename Electrons, typename Photons>
    std::tuple<Real, Real, Real, Real> get_integration_bounds(Electrons const& electrons,
                                                              Photons const& photons) noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Fill the input spectrum data for the IC photon spectrum.
     * @param grid The integrator grid
     * @param electrons The electron distribution
     * @param photons The photon distribution
     * @tparam Electrons Type of the electron distribution
     * @tparam Photons Type of the photon distribution
     * <!-- ************************************************************************************** -->
     */
    template <typename Electrons, typename Photons>
    void fill_input_spectrum(IntegratorGrid& grid, Electrons const& electrons, Photons const& photons) noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Fill the integration grid for the IC photon spectrum.
     * @param grid The integrator grid
     * @param KN Whether to use the Klein-Nishina cross-section
     * <!-- ************************************************************************************** -->
     */
    void fill_integration_grid(IntegratorGrid& grid, bool KN) noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Integrate the IC photon spectrum.
     * @param grid The integrator grid
     * <!-- ************************************************************************************** -->
     */
    void integrate_IC_spectrum(IntegratorGrid const& grid) noexcept;
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
ICPhotonGrid generate_IC_photons(SynElectronGrid const& electron, SynPhotonGrid const& photon, bool KN = true) noexcept;

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

//========================================================================================================
//                                  template function implementation
//========================================================================================================
template <typename Electrons, typename Photons>
void ICPhoton::compute_IC_spectrum(Electrons const& electrons, Photons const& photons, bool KN) noexcept {
    // Get integration boundaries
    auto [nu0_min, nu0_max, gamma_min, gamma_max] = get_integration_bounds(electrons, photons);

    // Construct an integration grid in nu0 and gamma
    IntegratorGrid grid(nu0_min, nu0_max, gamma_min, gamma_max);

    fill_input_spectrum(grid, electrons, photons);

    // Compute the differential contributions
    fill_integration_grid(grid, KN);

    // Perform final integration
    integrate_IC_spectrum(grid);

    // Convert to log space for interpolation
    log2_I_nu_ = xt::log2(I_nu_IC_);
    log2_nu_IC_ = xt::log2(nu_IC_);
}

template <typename Electrons, typename Photons>
std::tuple<Real, Real, Real, Real> ICPhoton::get_integration_bounds(Electrons const& electrons,
                                                                    Photons const& photons) noexcept {
    Real nu0_min = min(photons.nu_m, photons.nu_c, photons.nu_a) / 1e5;
    Real nu0_max = photons.nu_M * 5;

    Real gamma_min = min(electrons.gamma_m, electrons.gamma_c, electrons.gamma_a);
    Real gamma_max = electrons.gamma_M * 5;

    return std::make_tuple(nu0_min, nu0_max, gamma_min, gamma_max);
}

template <typename Electrons, typename Photons>
void ICPhoton::fill_input_spectrum(IntegratorGrid& grid, Electrons const& electrons, Photons const& photons) noexcept {
    // For each bin in nu0, compute the synchrotron intensity and column number density
    for (size_t i = 0; i < grid.num; i++) {
        grid.I_nu_syn(i) = photons.compute_I_nu(grid.nu(i));
        grid.column_den(i) = electrons.compute_column_num_den(grid.gamma(i));
    }
}
