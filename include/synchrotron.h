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

/********************************************************************************************************************
 * STRUCT: InverseComptonY
 * DESCRIPTION: Handles Inverse Compton Y parameter calculations and related threshold values.
 ********************************************************************************************************************/
struct InverseComptonY {
    // Constructors
    InverseComptonY(Real nu_m, Real nu_c, Real B,
                    Real Y_T) noexcept;  // Constructor with frequency thresholds and magnetic field
    InverseComptonY(Real Y_T) noexcept;  // Constructor with only Y_T parameter
    InverseComptonY() noexcept;          // Default constructor

    // Member variables
    Real nu_hat_m{0};     // Frequency threshold for minimum electrons
    Real nu_hat_c{0};     // Frequency threshold for cooling electrons
    Real gamma_hat_m{0};  // Lorentz factor threshold for minimum energy electrons
    Real gamma_hat_c{0};  // Lorentz factor threshold for cooling electrons
    Real Y_T{0};          // Thomson scattering Y parameter
    size_t regime{0};     // Indicator for the operating regime (1=fast IC cooling, 2=slow IC cooling, 3=special case)

    // Member functions
    Real compute_val_at_nu(Real nu, Real p) const;        // Computes effective Y parameter at frequency nu
    Real compute_val_at_gamma(Real gamma, Real p) const;  // Computes effective Y parameter at Lorentz factor gamma

    // Static member functions for combined Y parameters
    static Real compute_Y_Thompson(InverseComptonY const& Ys);  // Returns Y_T parameter
    static Real compute_Y_tilt_at_gamma(InverseComptonY const& Ys, Real gamma,
                                        Real p);                                   // Y parameter at specific gamma
    static Real compute_Y_tilt_at_nu(InverseComptonY const& Ys, Real nu, Real p);  // Y parameter at specific frequency
};

/********************************************************************************************************************
 * STRUCT: SynElectrons
 * DESCRIPTION: Represents synchrotron-emitting electrons in the comoving frame along with their energy distribution
 *              and properties.
 ********************************************************************************************************************/
struct SynElectrons {
    // All values in comoving frame
    Real I_nu_peak{0};       // Peak intensity at the characteristic frequency
    Real gamma_m{0};         // Minimum electron Lorentz factor
    Real gamma_c{0};         // Cooling electron Lorentz factor
    Real gamma_a{0};         // Self-absorption Lorentz factor
    Real gamma_M{0};         // Maximum electron Lorentz factor
    Real p{2.3};             // Power-law index for the electron energy distribution
    Real column_num_den{0};  // Normalized column number density
    Real Y_c{0};             // Inverse Compton Y parameter at cooling frequency
    size_t regime{0};        // Regime indicator (1-6, determines spectral shape)
    InverseComptonY Ys;      // InverseComptonY parameters for this electron population

    // Calculates the column number density for a given electron Lorentz factor
    Real compute_column_num_den(Real gamma) const;

   private:
    // Computes the electron energy spectrum at a given Lorentz factor
    inline Real compute_gamma_spectrum(Real gamma) const;
};

/********************************************************************************************************************
 * STRUCT: SynPhotons
 * DESCRIPTION: Represents synchrotron photons in the comoving frame and provides spectral functions.
 ********************************************************************************************************************/
struct SynPhotons {
    // All values in comoving frame
    Real nu_m{0};  // Characteristic frequency corresponding to gamma_m
    Real nu_c{0};  // Cooling frequency corresponding to gamma_c
    Real nu_a{0};  // Self-absorption frequency
    Real nu_M{0};  // Maximum photon frequency

    Real log2_I_nu_peak{0};          // Log2 of peak intensity (for computational efficiency)
    Real log2_nu_m{0};               // Log2 of nu_m
    Real log2_nu_c{0};               // Log2 of nu_c
    Real log2_nu_a{0};               // Log2 of nu_a
    Real log2_nu_M{0};               // Log2 of nu_M
    const SynElectrons* e{nullptr};  // Pointer to the associated SynElectrons

    // Returns the intensity at a given frequency nu
    Real compute_I_nu(Real nu) const;            // Linear intensity
    Real compute_log2_I_nu(Real log2_nu) const;  // Log2 intensity (for computational efficiency)

    // Updates internal constants used in the spectral calculations
    void update_constant();

   private:
    // Cached calculation constants for spectral computations
    // Optimized calculation constants
    Real C1_{0};  // Cached spectral coefficient 1
    Real C2_{0};  // Cached spectral coefficient 2
    Real C3_{0};  // Cached spectral coefficient 3

    // Log2 of calculation constants for faster computation
    Real log2_C1_{0};  //
    Real log2_C2_{0};  //
    Real log2_C3_{0};  //
    Real log2_C4_{0};  //

    // Computes the photon spectrum at a given frequency nu
    inline Real compute_spectrum(Real nu) const;            // Linear spectrum
    inline Real compute_log2_spectrum(Real log2_nu) const;  // Log2 spectrum
};

/********************************************************************************************************************
 * TYPE ALIASES
 * DESCRIPTION: Defines multi-dimensional grid types for Synchrotron Photons and Electrons.
 ********************************************************************************************************************/
using SynPhotonGrid = xt::xtensor<SynPhotons, 3>;      // 3D grid of synchrotron photons
using SynElectronGrid = xt::xtensor<SynElectrons, 3>;  // 3D grid of synchrotron electrons

/********************************************************************************************************************
 * FUNCTION PROTOTYPES: Synchrotron Grid Creation and Generation
 * DESCRIPTION: Functions to create and generate grids for Synchrotron electrons and photons.
 ********************************************************************************************************************/
// Creates and returns a new electron grid based on shock parameters
SynElectronGrid generate_syn_electrons(Shock const& shock, Real p, Real xi = 1);

// Populates an existing electron grid with values based on shock parameters
void generate_syn_electrons(SynElectronGrid& electrons, Shock const& shock, Real p, Real xi = 1);

// Creates and returns a new photon grid based on shock and electron grid
SynPhotonGrid generate_syn_photons(Shock const& shock, SynElectronGrid const& electrons);

// Populates an existing photon grid with values based on shock and electron grid
void generate_syn_photons(SynPhotonGrid& photons, Shock const& shock, SynElectronGrid const& electrons);

/********************************************************************************************************************
 * FUNCTION PROTOTYPES: Synchrotron Update and Parameter Calculation
 * DESCRIPTION: Functions for updating electron grids and calculating synchrotron parameters.
 ********************************************************************************************************************/
// Updates electron grid with new Y parameter values
void update_electrons_4Y(SynElectronGrid& e, Shock const& shock);

// Calculates cooling Lorentz factor based on comoving time, magnetic field, and IC parameters
Real compute_gamma_c(Real t_com, Real B, InverseComptonY const& Ys, Real p);

// Determines the peak Lorentz factor based on absorption, minimum, and cooling factors
Real compute_gamma_peak(Real gamma_a, Real gamma_m, Real gamma_c);

// Calculates synchrotron frequency for a given Lorentz factor and magnetic field
Real compute_syn_freq(Real gamma, Real B);
