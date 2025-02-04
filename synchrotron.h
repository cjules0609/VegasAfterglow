//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/
#ifndef _SYNCHROTRON_
#define _SYNCHROTRON_

/********************************************************************************************************************
 * INCLUDES
 * DESCRIPTION: Standard and project-specific headers required for synchrotron calculations.
 ********************************************************************************************************************/
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
    InverseComptonY(double nu_m, double nu_c, double B, double Y_T);
    InverseComptonY(double Y_T);
    InverseComptonY();

    // Member variables
    double nu_hat_m{0};     // Frequency threshold for minimum electrons
    double nu_hat_c{0};     // Frequency threshold for cooling electrons
    double gamma_hat_m{0};  // Lorentz factor threshold for minimum energy electrons
    double gamma_hat_c{0};  // Lorentz factor threshold for cooling electrons
    double Y_T{0};          // Thomson scattering Y parameter
    size_t regime{0};       // Indicator for the operating regime

    // Member functions
    double as_nu(double nu, double p) const;        // Computes based on frequency and power-law index
    double as_gamma(double gamma, double p) const;  // Computes based on Lorentz factor and power-law index

    // Static member functions for combined Y parameters
    static double Y_Thompson(std::vector<InverseComptonY> const& Ys);
    static double Y_tilt_gamma(std::vector<InverseComptonY> const& Ys, double gamma, double p);
    static double Y_tilt_nu(std::vector<InverseComptonY> const& Ys, double nu, double p);
};

/********************************************************************************************************************
 * STRUCT: SynElectrons
 * DESCRIPTION: Represents synchrotron electrons in the comoving frame along with their energy distribution and
 *properties.
 ********************************************************************************************************************/
struct SynElectrons {
    // all in comoving frame
    double I_nu_peak{0};              // Peak intensity at the characteristic frequency
    double gamma_m{0};                // Minimum electron Lorentz factor
    double gamma_c{0};                // Cooling electron Lorentz factor
    double gamma_a{0};                // Self-absorption Lorentz factor
    double gamma_M{0};                // Maximum electron Lorentz factor
    double p{2.3};                    // Power-law index for the electron energy distribution
    double column_num_den{0};         // Normalized column number density
    double gamma_N_peak;              // Lorentz factor at which the column density peaks
    double Y_c{0};                    // Compton Y parameter for electrons
    size_t regime{0};                 // Regime indicator
    std::vector<InverseComptonY> Ys;  // Collection of InverseComptonY parameters

    // Calculates the column number density for a given electron Lorentz factor
    double columnNumDen(double gamma) const;

   private:
    // Computes the electron energy spectrum at a given Lorentz factor
    inline double gammaSpectrum(double gamma) const;
};

/********************************************************************************************************************
 * STRUCT: SynPhotons
 * DESCRIPTION: Represents synchrotron photons in the comoving frame and provides spectral functions.
 ********************************************************************************************************************/
struct SynPhotons {
    // all in comoving frame
    double nu_m{0};                  // Characteristic frequency corresponding to gamma_m
    double nu_c{0};                  // Cooling frequency corresponding to gamma_c
    double nu_a{0};                  // Self-absorption frequency
    double nu_M{0};                  // Maximum photon frequency
    const SynElectrons* e{nullptr};  // Pointer to the associated SynElectrons

    // Returns the intensity at a given frequency nu
    double I_nu(double nu) const;
    // Updates internal constants used in the spectral calculations
    void updateConstant();

   private:
    double a_m_1_3{0};     // (nu_a / nu_m)^(1/3)
    double c_m_1_2{0};     // (nu_c / nu_m)^(1/2)
    double m_a_pa4_2{0};   // (nu_m / nu_a)^((p+4)/2)
    double a_m_mpa1_2{0};  // (nu_a / nu_m)^((-p+1)/2)
    double a_c_1_3{0};     // (nu_a / nu_c)^(1/3)
    double a_m_1_2{0};     // (nu_a / nu_m)^(1/2)
    double R4{0};          // R coefficient for case 4 (Bing Zhang's Book, page 199)
    double R6{0};          // R coefficient for case 6 (Bing Zhang's Book, page 200)

    // Computes the photon spectrum at a given frequency nu
    inline double spectrum(double nu) const;
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
SynElectronGrid createSynElectronGrid(size_t phi_size, size_t theta_size, size_t r_size);
SynElectronGrid genSynElectrons(Shock const& shock, double p, double xi = 1);

SynPhotonGrid createSynPhotonGrid(size_t phi_size, size_t theta_size, size_t r_size);
SynPhotonGrid genSynPhotons(Shock const& shock, SynElectronGrid const& electrons);

/********************************************************************************************************************
 * FUNCTION PROTOTYPES: Synchrotron Update and Parameter Calculation
 * DESCRIPTION: Functions for updating electron grids and calculating synchrotron parameters.
 ********************************************************************************************************************/
void updateElectrons4Y(SynElectronGrid& e, Shock const& shock);
double syn_gamma_c(double t_com, double B, std::vector<InverseComptonY> const& Ys, double p);
double syn_gamma_N_peak(double gamma_a, double gamma_m, double gamma_c);
double syn_nu(double gamma, double B);

#endif