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
    InverseComptonY(Real nu_m, Real nu_c, Real B, Real Y_T);
    InverseComptonY(Real Y_T);
    InverseComptonY();

    // Member variables
    Real nu_hat_m{0};     // Frequency threshold for minimum electrons
    Real nu_hat_c{0};     // Frequency threshold for cooling electrons
    Real gamma_hat_m{0};  // Lorentz factor threshold for minimum energy electrons
    Real gamma_hat_c{0};  // Lorentz factor threshold for cooling electrons
    Real Y_T{0};          // Thomson scattering Y parameter
    size_t regime{0};     // Indicator for the operating regime

    // Member functions
    Real as_nu(Real nu, Real p) const;        // Computes based on frequency and power-law index
    Real as_gamma(Real gamma, Real p) const;  // Computes based on Lorentz factor and power-law index

    // Static member functions for combined Y parameters
    static Real Y_Thompson(InverseComptonY const& Ys);
    static Real Y_tilt_gamma(InverseComptonY const& Ys, Real gamma, Real p);
    static Real Y_tilt_nu(InverseComptonY const& Ys, Real nu, Real p);
};

/********************************************************************************************************************
 * STRUCT: SynElectrons
 * DESCRIPTION: Represents synchrotron electrons in the comoving frame along with their energy distribution and
 *properties.
 ********************************************************************************************************************/
struct SynElectrons {
    // all in comoving frame
    Real I_nu_peak{0};       // Peak intensity at the characteristic frequency
    Real gamma_m{0};         // Minimum electron Lorentz factor
    Real gamma_c{0};         // Cooling electron Lorentz factor
    Real gamma_a{0};         // Self-absorption Lorentz factor
    Real gamma_M{0};         // Maximum electron Lorentz factor
    Real p{2.3};             // Power-law index for the electron energy distribution
    Real column_num_den{0};  // Normalized column number density
    Real Y_c{0};             // Compton Y parameter for electrons
    size_t regime{0};        // Regime indicator
    InverseComptonY Ys;      // InverseComptonY parameters

    // Calculates the column number density for a given electron Lorentz factor
    Real columnNumDen(Real gamma) const;

   private:
    // the electron energy spectrum at a given Lorentz factor
    inline Real gammaSpectrum(Real gamma) const;
};

/********************************************************************************************************************
 * STRUCT: SynPhotons
 * DESCRIPTION: Represents synchrotron photons in the comoving frame and provides spectral functions.
 ********************************************************************************************************************/
struct SynPhotons {
    // all in comoving frame
    Real nu_m{0};                    // Characteristic frequency corresponding to gamma_m
    Real nu_c{0};                    // Cooling frequency corresponding to gamma_c
    Real nu_a{0};                    // Self-absorption frequency
    Real nu_M{0};                    // Maximum photon frequency
    const SynElectrons* e{nullptr};  // Pointer to the associated SynElectrons

    // Returns the intensity at a given frequency nu
    Real I_nu(Real nu) const;
    // Updates internal constants used in the spectral calculations
    void updateConstant();

   private:
    /*Real a_m_1_3{0};     // (nu_a / nu_m)^(1/3)
    Real c_m_1_2{0};     // (nu_c / nu_m)^(1/2)
    Real m_a_pa4_2{0};   // (nu_m / nu_a)^((p+4)/2)
    Real a_m_mpa1_2{0};  // (nu_a / nu_m)^((-p+1)/2)
    Real a_c_1_3{0};     // (nu_a / nu_c)^(1/3)
    Real a_m_1_2{0};     // (nu_a / nu_m)^(1/2)
    Real R4{0};          // R coefficient for case 4 (Bing Zhang's Book, page 199)
    Real R6{0};          // R coefficient for case 6 (Bing Zhang's Book, page 200)*/
    Real C1_{0};
    Real C2_{0};
    Real C3_{0};

    // Computes the photon spectrum at a given frequency nu
    inline Real spectrum(Real nu) const;
};

/********************************************************************************************************************
 * TYPE ALIASES
 * DESCRIPTION: Defines multi-dimensional grid types for Synchrotron Photons and Electrons.
 ********************************************************************************************************************/
using SynPhotonGrid = boost::multi_array<SynPhotons, 3>;
using SynElectronGrid = boost::multi_array<SynElectrons, 3>;

/********************************************************************************************************************
 * FUNCTION PROTOTYPES: Synchrotron Grid Creation and Generation
 * DESCRIPTION: Functions to create and generate grids for Synchrotron electrons and photons.
 ********************************************************************************************************************/
SynElectronGrid createSynElectronGrid(size_t phi_size, size_t theta_size, size_t t_size);
SynElectronGrid genSynElectrons(Shock const& shock, Real p, Real xi = 1);

SynPhotonGrid createSynPhotonGrid(size_t phi_size, size_t theta_size, size_t t_size);
SynPhotonGrid genSynPhotons(Shock const& shock, SynElectronGrid const& electrons);

/********************************************************************************************************************
 * FUNCTION PROTOTYPES: Synchrotron Update and Parameter Calculation
 * DESCRIPTION: Functions for updating electron grids and calculating synchrotron parameters.
 ********************************************************************************************************************/
void updateElectrons4Y(SynElectronGrid& e, Shock const& shock);
Real syn_gamma_c(Real t_com, Real B, InverseComptonY const& Ys, Real p);
Real syn_gamma_N_peak(Real gamma_a, Real gamma_m, Real gamma_c);
Real syn_nu(Real gamma, Real B);
