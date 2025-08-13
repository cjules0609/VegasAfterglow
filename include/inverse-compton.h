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
#include "mesh.h"
#include "shock.h"
#include "utilities.h"

template <std::size_t N>
using StackArray = xt::xtensor_fixed<Real, xt::xshape<N>>;

template <std::size_t N, std::size_t M>
using StackMesh = xt::xtensor_fixed<Real, xt::xshape<N, M>>;

/**
 * <!-- ************************************************************************************** -->
 * @struct InverseComptonY
 * @brief Handles Inverse Compton Y parameter calculations and related threshold values.
 * <!-- ************************************************************************************** -->
 */
struct InverseComptonY {
    /**
     * <!-- ************************************************************************************** -->
     * @brief Initializes an InverseComptonY object with frequency thresholds, magnetic field and Y parameter.
     * @details Computes characteristic gamma values and corresponding frequencies, then determines cooling regime.
     * @param nu_m Characteristic frequency for minimum Lorentz factor
     * @param nu_c Characteristic frequency for cooling Lorentz factor
     * @param B Magnetic field strength
     * @param Y_T Thomson Y parameter
     * <!-- ************************************************************************************** -->
     */
    InverseComptonY(Real nu_m, Real nu_c, Real B, Real Y_T) noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Simple constructor that initializes with only the Thomson Y parameter for special cases.
     * @param Y_T Thomson Y parameter
     * <!-- ************************************************************************************** -->
     */
    InverseComptonY(Real Y_T) noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Default constructor that initializes all member variables to zero.
     * <!-- ************************************************************************************** -->
     */
    InverseComptonY() noexcept;

    // Member variables
    Real nu_hat_m{0};     ///< Frequency threshold for minimum electrons
    Real nu_hat_c{0};     ///< Frequency threshold for cooling electrons
    Real gamma_hat_m{0};  ///< Lorentz factor threshold for minimum energy electrons
    Real gamma_hat_c{0};  ///< Lorentz factor threshold for cooling electrons
    Real Y_T{0};          ///< Thomson scattering Y parameter
    size_t regime{0};     ///< Indicator for the operating regime (1=fast IC cooling, 2=slow IC cooling, 3=special case)

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the effective Y parameter for a given frequency and spectral index.
     * @details Different scaling relations apply depending on the cooling regime and frequency range.
     * @param nu Frequency at which to compute the Y parameter
     * @param p Spectral index of electron distribution
     * @return The effective Y parameter at the given frequency
     * <!-- ************************************************************************************** -->
     */
    Real compute_val_at_nu(Real nu, Real p) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the effective Y parameter for a given Lorentz factor and spectral index.
     * @details Different scaling relations apply depending on the cooling regime and gamma value.
     * @param gamma Electron Lorentz factor
     * @param p Spectral index of electron distribution
     * @return The effective Y parameter at the given gamma
     * <!-- ************************************************************************************** -->
     */
    Real compute_val_at_gamma(Real gamma, Real p) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Returns the Thomson Y parameter from the provided InverseComptonY object.
     * @details Previously supported summing Y parameters from multiple objects.
     * @param Ys InverseComptonY object
     * @return The Thomson Y parameter
     * <!-- ************************************************************************************** -->
     */
    static Real compute_Y_Thompson(InverseComptonY const& Ys);  ///< Returns Y_T parameter

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the effective Y parameter at a specific Lorentz factor and spectral index.
     * @details Previously supported summing contributions from multiple InverseComptonY objects.
     * @param Ys InverseComptonY object
     * @param gamma Electron Lorentz factor
     * @param p Spectral index of electron distribution
     * @return The effective Y parameter at the given gamma
     * <!-- ************************************************************************************** -->
     */
    static Real compute_Y_tilt_at_gamma(InverseComptonY const& Ys, Real gamma, Real p);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the effective Y parameter at a specific frequency and spectral index.
     * @details Previously supported summing contributions from multiple InverseComptonY objects.
     * @param Ys InverseComptonY object
     * @param nu Frequency at which to compute the Y parameter
     * @param p Spectral index of electron distribution
     * @return The effective Y parameter at the given frequency
     * <!-- ************************************************************************************** -->
     */
    static Real compute_Y_tilt_at_nu(InverseComptonY const& Ys, Real nu, Real p);
};

/**
 * <!-- ************************************************************************************** -->
 * @defgroup IC_Calculation Inverse Compton Calculation Constants and Functions
 * @brief Constants and inline functions for inverse Compton photon calculations
 * <!-- ************************************************************************************** -->
 */

/// A constant used in the IC photon calculation, defined as √2/3.
inline constexpr Real IC_x0 = 0.47140452079103166;

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the Compton scattering cross-section as a function of frequency (nu).
 * @param nu The frequency at which to compute the cross-section
 * @return Compton cross-section
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

    static constexpr size_t num{128};   ///< Number of bins.
    StackArray<num + 1> nu_edge{0};     ///< Bin edges for nu.
    StackArray<num + 1> gamma_edge{0};  ///< Bin edges for gamma.
    StackArray<num> nu{0};              ///< Center values for nu.
    StackArray<num> gamma{0};           ///< Center values for gamma.
    StackArray<num> I_nu{0};            ///< Specific intensity at each nu center.
    StackArray<num> column_num_den{0};  ///< column Number density at each gamma center.
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
     * @brief Returns the photon specific intensity.
     * @param nu The frequency at which to compute the specific intensity
     * @return The specific intensity at the given frequency
     * <!-- ************************************************************************************** -->
     */
    Real compute_I_nu(Real nu) const noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the base-2 logarithm of the photon specific intensity at a given frequency.
     * @param log2_nu The base-2 logarithm of the frequency
     * @return The base-2 logarithm of the photon specific intensity at the given frequency
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

   private:
    StackArray<spectrum_resol> P_nu_IC_{0};     ///< IC photon spectrum specific power array.
    StackArray<spectrum_resol> nu_IC_{0};       ///< Frequency grid for the IC photon spectrum.
    StackArray<spectrum_resol> log2_P_nu_{0};   ///< Base-2 logarithm of the IC photon spectrum specific power array.
    StackArray<spectrum_resol> log2_nu_IC_{0};  ///< Base-2 logarithm of the frequency grid for the IC photon spectrum.

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

    /**
     * <!-- ************************************************************************************** -->
     * @brief Remove the zero tail of the IC photon spectrum.
     * <!-- ************************************************************************************** -->
     */
    void remove_zero_tail() noexcept;
};

class SynPhotons;
class SynElectrons;

/// Type alias for 3D grid of synchrotron photons
using SynPhotonGrid = xt::xtensor<SynPhotons, 3>;

/// Type alias for 3D grid of synchrotron electrons
using SynElectronGrid = xt::xtensor<SynElectrons, 3>;

/// Defines a 3D grid (using xt::xtensor) for storing ICPhoton objects.
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
    // IntegratorGrid grid(nu0_min, nu0_max, gamma_min, gamma_max);
    IntegratorGrid grid(nu0_min, nu0_max, gamma_min, gamma_max);

    fill_input_spectrum(grid, electrons, photons);

    // Compute the differential contributions
    fill_integration_grid(grid, KN);

    // Perform final integration
    integrate_IC_spectrum(grid);

    remove_zero_tail();

    // Convert to log space for interpolation
    log2_P_nu_ = xt::log2(P_nu_IC_);
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
        grid.I_nu(i) = photons.compute_I_nu(grid.nu(i));
        grid.column_num_den(i) = electrons.compute_column_den(grid.gamma(i));
    }
}
