//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/
#pragma once

#include <vector>

#include "medium.h"
#include "mesh.h"
#include "shock.h"

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
 * @struct SynElectrons
 * @brief Represents synchrotron-emitting electrons in the comoving frame along with their energy distribution
 *        and properties.
 * <!-- ************************************************************************************** -->
 */
struct SynElectrons {
    // All values in comoving frame
    Real I_nu_peak{0};       ///< Peak intensity at the characteristic frequency
    Real gamma_m{0};         ///< Minimum electron Lorentz factor
    Real gamma_c{0};         ///< Cooling electron Lorentz factor
    Real gamma_a{0};         ///< Self-absorption Lorentz factor
    Real gamma_M{0};         ///< Maximum electron Lorentz factor
    Real p{2.3};             ///< Power-law index for the electron energy distribution
    Real column_num_den{0};  ///< Normalized column number density
    Real Y_c{0};             ///< Inverse Compton Y parameter at cooling frequency
    size_t regime{0};        ///< Regime indicator (1-6, determines spectral shape)
    InverseComptonY Ys;      ///< InverseComptonY parameters for this electron population

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the electron column number density at a specific Lorentz factor.
     * @details Includes corrections for inverse Compton cooling effects above the cooling Lorentz factor.
     * @param gamma Electron Lorentz factor
     * @return Column number density at the specified Lorentz factor
     * <!-- ************************************************************************************** -->
     */
    Real compute_column_num_den(Real gamma) const;

   private:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the electron energy spectrum at a given Lorentz factor.
     * @details Different spectral forms apply based on the current regime and relative to
     *          characteristic Lorentz factors (gamma_a, gamma_c, gamma_m, gamma_M).
     * @param gamma Electron Lorentz factor
     * @return The normalized electron energy spectrum value
     * <!-- ************************************************************************************** -->
     */
    inline Real compute_gamma_spectrum(Real gamma) const;
};

/**
 * <!-- ************************************************************************************** -->
 * @struct SynPhotons
 * @brief Represents synchrotron photons in the comoving frame and provides spectral functions.
 * <!-- ************************************************************************************** -->
 */
struct SynPhotons {
    // All values in comoving frame
    Real nu_m{0};  ///< Characteristic frequency corresponding to gamma_m
    Real nu_c{0};  ///< Cooling frequency corresponding to gamma_c
    Real nu_a{0};  ///< Self-absorption frequency
    Real nu_M{0};  ///< Maximum photon frequency

    Real log2_I_nu_peak{0};          ///< Log2 of peak intensity (for computational efficiency)
    Real log2_nu_m{0};               ///< Log2 of nu_m
    Real log2_nu_c{0};               ///< Log2 of nu_c
    Real log2_nu_a{0};               ///< Log2 of nu_a
    Real log2_nu_M{0};               ///< Log2 of nu_M
    const SynElectrons* e{nullptr};  ///< Pointer to the associated SynElectrons

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the synchrotron intensity at a given frequency.
     * @details Includes inverse Compton corrections for frequencies above the cooling frequency.
     * @param nu Frequency at which to compute the intensity
     * @return The synchrotron intensity at the specified frequency
     * <!-- ************************************************************************************** -->
     */
    Real compute_I_nu(Real nu) const;  ///< Linear intensity

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the base-2 logarithm of synchrotron intensity at a given frequency.
     * @details Optimized for numerical computation by using logarithmic arithmetic.
     * @param log2_nu Base-2 logarithm of the frequency
     * @return Base-2 logarithm of synchrotron intensity
     * <!-- ************************************************************************************** -->
     */
    Real compute_log2_I_nu(Real log2_nu) const;  ///<  Log2 intensity (for computational efficiency)

    /**
     * <!-- ************************************************************************************** -->
     * @brief Updates cached calculation constants used for efficiently computing synchrotron spectra.
     * @details Constants vary based on the electron regime (1-6) and involve different power laws.
     * <!-- ************************************************************************************** -->
     */
    void update_constant();

   private:
    // Cached calculation constants for spectral computations
    // Optimized calculation constants
    Real C1_{0};  ///< Cached spectral coefficient 1
    Real C2_{0};  ///< Cached spectral coefficient 2
    Real C3_{0};  ///< Cached spectral coefficient 3

    // Log2 of calculation constants for faster computation
    Real log2_C1_{0};
    Real log2_C2_{0};
    Real log2_C3_{0};
    Real log2_C4_{0};

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the synchrotron spectrum at a given frequency based on the electron regime.
     * @details Implements the broken power-law with exponential cutoff formulae for different regimes.
     * @param nu The frequency at which to compute the spectrum
     * @return The normalized synchrotron spectrum value
     * <!-- ************************************************************************************** -->
     */
    inline Real compute_spectrum(Real nu) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the base-2 logarithm of synchrotron spectrum at a given frequency.
     * @details Uses logarithmic arithmetic for numerical stability in different spectral regimes.
     * @param log2_nu Base-2 logarithm of the frequency
     * @return Base-2 logarithm of the synchrotron spectrum
     * <!-- ************************************************************************************** -->
     */
    inline Real compute_log2_spectrum(Real log2_nu) const;
};

/**
 * <!-- ************************************************************************************** -->
 * @defgroup SynchrotronGrids Synchrotron Grid Type Aliases
 * @brief Defines multi-dimensional grid types for Synchrotron Photons and Electrons.
 * <!-- ************************************************************************************** -->
 */

/// Type alias for 3D grid of synchrotron photons
using SynPhotonGrid = xt::xtensor<SynPhotons, 3>;
/// Type alias for 3D grid of synchrotron electrons
using SynElectronGrid = xt::xtensor<SynElectrons, 3>;

/**
 * <!-- ************************************************************************************** -->
 * @defgroup SynchrotronFunctions Synchrotron Grid Creation and Generation
 * @brief Functions to create and generate grids for Synchrotron electrons and photons.
 * <!-- ************************************************************************************** -->
 */

/**
 * <!-- ************************************************************************************** -->
 * @brief Creates and returns a new electron grid based on shock parameters
 * @details Initializes all electron properties including Lorentz factors, column densities,
 *          and peak intensities for each grid cell.
 * @param shock The shock object containing physical properties
 * @param p Power-law index for the electron energy distribution
 * @param xi Energy partitioning parameter (default: 1)
 * @return A new grid of synchrotron electrons
 * <!-- ************************************************************************************** -->
 */
SynElectronGrid generate_syn_electrons(Shock const& shock, Real p, Real xi = 1);

/**
 * <!-- ************************************************************************************** -->
 * @brief Populates an existing electron grid with values based on shock parameters
 * @details Modifies a grid supplied by the caller rather than creating a new one.
 * @param electrons The electron grid to populate
 * @param shock The shock object containing physical properties
 * @param p Power-law index for the electron energy distribution
 * @param xi Energy partitioning parameter (default: 1)
 * <!-- ************************************************************************************** -->
 */
void generate_syn_electrons(SynElectronGrid& electrons, Shock const& shock, Real p, Real xi = 1);

/**
 * <!-- ************************************************************************************** -->
 * @brief Creates and returns a new photon grid based on shock and electron grid
 * @details Computes characteristic frequencies and updates calculation constants for each grid cell.
 *          Returns the populated photon grid.
 * @param shock The shock object containing physical properties
 * @param electrons The electron grid providing energy distribution information
 * @return A new grid of synchrotron photons
 * <!-- ************************************************************************************** -->
 */
SynPhotonGrid generate_syn_photons(Shock const& shock, SynElectronGrid const& electrons);

/**
 * <!-- ************************************************************************************** -->
 * @brief Populates an existing photon grid with values based on shock and electron grid
 * @param photons The photon grid to populate
 * @param shock The shock object containing physical properties
 * @param electrons The electron grid providing energy distribution information
 * <!-- ************************************************************************************** -->
 */
void generate_syn_photons(SynPhotonGrid& photons, Shock const& shock, SynElectronGrid const& electrons);

/**
 * <!-- ************************************************************************************** -->
 * @defgroup SynchrotronUpdates Synchrotron Update and Parameter Calculation
 * @brief Functions for updating electron grids and calculating synchrotron parameters.
 * <!-- ************************************************************************************** -->
 */

/**
 * <!-- ************************************************************************************** -->
 * @brief Updates electron properties throughout the grid based on shock parameters and IC Y values.
 * @details Recalculates gamma_M, gamma_c, gamma_a, regime, and Y_c parameters for each grid cell.
 *          Handles both freshly-shocked and adiabatic cooling regions.
 * @param e Synchrotron electron grid to update
 * @param shock The shock object containing evolution data
 * <!-- ************************************************************************************** -->
 */
void update_electrons_4Y(SynElectronGrid& e, Shock const& shock);

/**
 * <!-- ************************************************************************************** -->
 * @brief Calculates cooling Lorentz factor based on comoving time, magnetic field, and IC parameters
 * @details Accounts for synchrotron and inverse Compton cooling using an iterative approach
 *          to handle the Lorentz factor dependent IC cooling.
 * @param t_com Comoving time
 * @param B Magnetic field
 * @param Ys Inverse Compton Y parameters
 * @param p Power-law index for electron energy distribution
 * @return The cooling Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real compute_gamma_c(Real t_com, Real B, InverseComptonY const& Ys, Real p);

/**
 * <!-- ************************************************************************************** -->
 * @brief Determines the electron Lorentz factor at which the number density peaks.
 * @details Based on the relative ordering of absorption, minimum, and cooling Lorentz factors.
 * @param gamma_a Absorption Lorentz factor
 * @param gamma_m Minimum electron Lorentz factor
 * @param gamma_c Cooling electron Lorentz factor
 * @return Peak Lorentz factor
 * <!-- ************************************************************************************** -->
 */
Real compute_gamma_peak(Real gamma_a, Real gamma_m, Real gamma_c);

/**
 * <!-- ************************************************************************************** -->
 * @brief Calculates synchrotron frequency for a given Lorentz factor and magnetic field
 * @param gamma Electron Lorentz factor
 * @param B Magnetic field
 * @return The synchrotron frequency
 * <!-- ************************************************************************************** -->
 */
Real compute_syn_freq(Real gamma, Real B);

