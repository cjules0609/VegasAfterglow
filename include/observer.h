//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/
#pragma once

#include <iostream>
#include <stdexcept>
#include <thread>
#include <vector>

#include "afterglow.h"
#include "macros.h"
#include "mesh.h"

/********************************************************************************************************************
 * CLASS: LogScaleInterp
 * DESCRIPTION: Provides logarithmic scale interpolation for various shock-related quantities. This class is used
 *              to interpolate physical quantities (radius, intensity, Doppler factor) in logarithmic space, which
 *              is particularly useful for GRB afterglow calculations where quantities often span many orders of
 *              magnitude. The interpolation is performed using the logarithm of the observation time as the
 *              independent variable.
 ********************************************************************************************************************/
class LogScaleInterp {
   public:
    size_t jet_3d{0};  // Flag indicating if the jet is non-axis-symmetric (non-zero if true)

    // Interpolates the luminosity using the observation time (t_obs) in logarithmic space
    // Returns the interpolated luminosity value
    Real interpLuminosity(Real t_obs) const noexcept;

    // Validates and sets the interpolation boundaries using the provided grids
    // Returns true if both lower and upper boundaries are valid for interpolation
    template <typename... PhotonGrid>
    bool validateInterpBoundary(size_t i, size_t j, size_t k, Real z, MeshGrid3d const& surface,
                                MeshGrid3d const& t_obs, MeshGrid3d const& doppler, Real nu_obs,
                                PhotonGrid const&... photons) noexcept;

   private:
    Real slope{0};

    Real t_obs_lo{0};  // Lower boundary of observation time
    Real log_I_lo{0};  // Lower boundary of specific intensity
    Real log_I_hi{0};  // Upper boundary of specific intensity
    Real surface_lo{0};

    Real nu_last{0};
    size_t idx_hi{0};  // Index for the upper boundary in the grid
    size_t idx_i{0};   // Index for phi coordinate
    size_t idx_j{0};   // Index for theta coordinate
};

/********************************************************************************************************************
 * CLASS: Observer
 * DESCRIPTION: Represents an observer in the GRB afterglow simulation. This class handles the calculation of
 *              observed quantities such as specific flux, integrated flux, and spectra. It accounts for relativistic
 *              effects (Doppler boosting), cosmological effects (redshift), and geometric effects (solid angle).
 *              The observer can be placed at any viewing angle relative to the jet axis.
 ********************************************************************************************************************/
class Observer {
   public:
    // Default constructor
    Observer() = default;

    // Grids for storing simulation data
    MeshGrid3d t_obs_grid;  // Grid of observation times
    MeshGrid3d doppler;     // Grid of Doppler factors
    MeshGrid3d surface;     // Grid of effective emission surface L_obs = I^\prime * S

    // Physical parameters
    Real lumi_dist{1};  // Luminosity distance
    Real z{0};          // Redshift

    // Main observation function that sets up the observer's view of the simulation
    void observe(Coord const& coord, Shock const& shock, Real luminosity_dist, Real redshift);

    void observeAt(Array const& t_obs, Coord const& coord, Shock& shock, Real luminosity_dist, Real redshift);

    // Computes the specific flux at a single observed frequency
    template <typename... PhotonGrid>
    Array specificFlux(Array const& t_obs, Real nu_obs, PhotonGrid const&... photons);

    // Computes the specific flux at multiple observed frequencies
    template <typename... PhotonGrid>
    MeshGrid specificFlux(Array const& t_obs, Array const& nu_obs, PhotonGrid const&... photons);

    // Computes the integrated flux over a frequency band
    template <typename... PhotonGrid>
    Array flux(Array const& t_obs, Array const& band_freq, PhotonGrid const&... photons);

    // Computes the spectrum at multiple frequencies for a single observation time
    template <typename... PhotonGrid>
    MeshGrid spectrum(Array const& freqs, Array const& t_obs, PhotonGrid const&... photons);

    // Computes the spectrum at multiple frequencies for a single observation time
    template <typename... PhotonGrid>
    Array spectrum(Array const& freqs, Real t_obs, PhotonGrid const&... photons);

    // Updates the required grid points for observation
    void updateRequired(MaskGrid& required, Array const& t_obs);

   private:
    LogScaleInterp interp;   // Helper class for logarithmic interpolation
    size_t eff_phi_size{1};  // Effective number of phi grid points
    size_t theta_size{0};    // Number of theta grid points
    size_t t_size{0};        // Number of time grid points

    void buildObsTimeGrid(Coord const& coord, Shock const& shock, Real luminosity_dist, Real redshift);

    // Calculates the observation time grid based on Lorentz factor and engine time
    void calcObsTimeGrid(Coord const& coord, MeshGrid3d const& Gamma, MeshGrid3d const& r_grid);

    // Calculates the effective emission surface for each grid point
    void calcEmissionSurface(Coord const& coord, Shock const& shock);
};

/********************************************************************************************************************
 * TEMPLATE METHOD: LogScaleInterp::validateInterpBoundary
 * DESCRIPTION: Attempts to set the lower and upper boundary values for logarithmic interpolation.
 *              It updates the internal boundary members for:
 *                - Logarithmic observation time (t_obs_lo, log_t_ratio)
 *                - Logarithmic luminosity (L_lo, L_hi)
 *              The boundaries are set using data from the provided grids and the photon grids (via parameter pack).
 *              Returns true if both lower and upper boundaries are finite such that interpolation can proceed.
 ********************************************************************************************************************/
template <typename... PhotonGrid>
bool LogScaleInterp::validateInterpBoundary(size_t i, size_t j, size_t k_lo, Real z, MeshGrid3d const& surface,
                                            MeshGrid3d const& t_obs, MeshGrid3d const& doppler, Real nu_obs,
                                            PhotonGrid const&... photons) noexcept {
    t_obs_lo = t_obs(i, j, k_lo);
    Real log_t_ratio = fastLog2(t_obs(i, j, k_lo + 1) / t_obs_lo);

    if (!std::isfinite(log_t_ratio) || log_t_ratio == 0) [[unlikely]] {
        return false;
    }

    size_t phi_idx = i * jet_3d;
    surface_lo = surface(i, j, k_lo);
    Real log_S_ratio = fastLog2(surface(i, j, k_lo + 1) / surface_lo);

    // continuing from previous boundary, shift the high boundary to lower. Calling .I_nu()/.log_I_nu() could be
    // expensive.
    if (idx_hi != 0 && k_lo == idx_hi && idx_i == i && idx_j == j && nu_last == nu_obs) {
        log_I_lo = log_I_hi;
    } else {
        Real nu = (1 + z) * nu_obs / doppler(i, j, k_lo);
        log_I_lo = (photons(phi_idx, j, k_lo).log2_I_nu(nu) + ...);
    }

    Real nu = (1 + z) * nu_obs / doppler(i, j, k_lo + 1);
    log_I_hi = (photons(phi_idx, j, k_lo + 1).log2_I_nu(nu) + ...);

    slope = (log_I_hi - log_I_lo + log_S_ratio) / log_t_ratio;

    if (!std::isfinite(slope)) [[unlikely]] {
        return false;
    }
    nu_last = nu_obs;
    idx_hi = k_lo + 1;
    idx_i = i;
    idx_j = j;
    return true;
}

/********************************************************************************************************************
 * TEMPLATE METHOD: Observer::specificFlux (multi-frequency overload)
 * DESCRIPTION: Returns the specific flux (as a MeshGrid) for multiple observed frequencies (nu_obs) by computing
 *              the specific flux for each frequency and assembling the results into a grid.
 ********************************************************************************************************************/
inline size_t iterate_to(Real value, Array const& arr, size_t it) {
    while (it < arr.size() && arr(it) < value) {
        it++;
    }
    return it;
}

template <typename... PhotonGrid>
MeshGrid Observer::specificFlux(Array const& t_obs, Array const& nu_obs, PhotonGrid const&... photons) {
    size_t t_obs_size = t_obs.size();
    size_t nu_size = nu_obs.size();

    MeshGrid F_nu({nu_size, t_obs_size}, 0);

    // Loop over effective phi and theta grid points.
    for (size_t i = 0; i < eff_phi_size; i++) {
        for (size_t j = 0; j < theta_size; j++) {
            // Skip observation times that are below the grid's start time
            size_t t_idx = iterate_to(t_obs_grid(i, j, 0), t_obs, 0);

            // Interpolate for observation times within the grid.
            for (size_t k = 0; k < t_size - 1 && t_idx < t_obs_size; k++) {
                Real const t_lo = t_obs_grid(i, j, k);
                Real const t_hi = t_obs_grid(i, j, k + 1);

                if (t_hi < t_obs(t_idx)) {
                    continue;
                }

                if (t_lo <= t_obs(t_idx) && t_obs(t_idx) < t_hi) {
                    for (size_t l = 0; l < nu_size; l++) {
                        if (interp.validateInterpBoundary(i, j, k, z, surface, t_obs_grid, doppler, nu_obs[l],
                                                          photons...)) {
                            for (size_t tt = t_idx; tt < t_obs_size && t_obs(tt) < t_hi; tt++) {
                                F_nu(l, tt) += interp.interpLuminosity(t_obs(tt));
                            }
                        }
                    }
                    t_idx = iterate_to(t_hi, t_obs, t_idx);
                }
            }
        }
    }

    // Normalize the flux by the factor (1+z)/(lumi_dist^2).
    Real const coef = (1 + z) / (lumi_dist * lumi_dist);
    F_nu *= coef;

    return F_nu;
}

/********************************************************************************************************************
 * TEMPLATE METHOD: Observer::specificFlux (single-frequency overload)
 * DESCRIPTION: Returns the specific flux (as an Array) for a single observed frequency (nu_obs) by computing the
 *              specific flux over the observation times.
 ********************************************************************************************************************/
template <typename... PhotonGrid>
Array Observer::specificFlux(Array const& t_obs, Real nu_obs, PhotonGrid const&... photons) {
    return xt::view(specificFlux(t_obs, Array({nu_obs}), photons...), 0);
}

template <typename... PhotonGrid>
Array Observer::spectrum(Array const& freqs, Real t_obs, PhotonGrid const&... photons) {
    return xt::view(spectrum(freqs, Array({t_obs}), photons...), 0);
}

template <typename... PhotonGrid>
MeshGrid Observer::spectrum(Array const& freqs, Array const& t_obs, PhotonGrid const&... photons) {
    return xt::transpose(specificFlux(t_obs, freqs, photons...));
}

/********************************************************************************************************************
 * TEMPLATE METHOD: Observer::flux
 * DESCRIPTION: Computes the integrated flux over a frequency band specified by band_freq.
 *              It converts band boundaries to center frequencies, computes the specific flux at each frequency,
 *              and integrates (sums) the flux contributions weighted by the frequency bin widths.
 ********************************************************************************************************************/
template <typename... PhotonGrid>
Array Observer::flux(Array const& t_obs, Array const& band_freq, PhotonGrid const&... photons) {
    Array nu_obs = boundaryToCenterLog(band_freq);
    MeshGrid F_nu = specificFlux(t_obs, nu_obs, photons...);
    Array flux({t_obs.size()}, 0);
    for (size_t i = 0; i < nu_obs.size(); ++i) {
        Real dnu = band_freq(i + 1) - band_freq(i);
        for (size_t j = 0; j < flux.size(); ++j) {
            flux(j) += dnu * F_nu(i, j);
        }
    }
    return flux;
}
