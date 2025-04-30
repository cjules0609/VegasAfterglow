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
    MeshGrid3d time;     // Grid of observation times
    MeshGrid3d doppler;  // Grid of Doppler factors
    MeshGrid3d surface;  // Grid of effective emission surface L_obs = I^\prime * S

    // Physical parameters
    Real lumi_dist{1};   // Luminosity distance
    Real one_plus_z{0};  // Redshift z+1

    // Main observation function that sets up the observer's view of the simulation
    void observe(Coord const& coord, Shock const& shock, Real luminosity_dist, Real redshift);

    // Configures observer to observe at specific time points
    void observe_at(Array const& t_obs, Coord const& coord, Shock& shock, Real luminosity_dist, Real redshift);

    // Computes the specific flux at a single observed frequency
    template <typename... PhotonGrid>
    Array specific_flux(Array const& t_obs, Real nu_obs, PhotonGrid const&... photons);

    // Computes the specific flux at multiple observed frequencies
    template <typename... PhotonGrid>
    MeshGrid specific_flux(Array const& t_obs, Array const& nu_obs, PhotonGrid const&... photons);

    // Computes the integrated flux over a frequency band
    template <typename... PhotonGrid>
    Array flux(Array const& t_obs, Array const& band_freq, PhotonGrid const&... photons);

    // Computes the spectrum at multiple observation times
    template <typename... PhotonGrid>
    MeshGrid spectra(Array const& freqs, Array const& t_obs, PhotonGrid const&... photons);

    // Computes the spectrum at single observation time
    template <typename... PhotonGrid>
    Array spectrum(Array const& freqs, Real t_obs, PhotonGrid const&... photons);

    // Updates the required grid points for observation
    void update_required(MaskGrid& required, Array const& t_obs);

   private:
    // for observer gird
    size_t jet_3d{0};        // Flag indicating if the jet is non-axis-symmetric (non-zero if true)
    size_t eff_phi_size{1};  // Effective number of phi grid points
    size_t theta_size{0};    // Number of theta grid points
    size_t t_size{0};        // Number of time grid points

    // Builds the observation time, doppler, and surface grids
    void build_time_grid(Coord const& coord, Shock const& shock, Real luminosity_dist, Real redshift);

    // Calculates the observation time grid based on Lorentz factor and engine time
    void calc_t_obs(Coord const& coord, Shock const& shock);

    // Calculates the effective emission surface for each grid point
    void calc_emission_surface(Coord const& coord, Shock const& shock);

    struct InterpState {
        Real slope{0};      // Slope for logarithmic interpolation
        Real log_I_lo{0};   // Lower boundary of specific intensity (log2 scale)
        Real log_I_hi{0};   // Upper boundary of specific intensity (log2 scale)
        size_t last_hi{0};  // Index for the upper boundary in the grid
    };

    // Interpolates the luminosity using the observation time (t_obs) in logarithmic space
    // Returns the interpolated luminosity value
    Real interpolate(InterpState const& state, size_t i, size_t j, size_t k, Real t_obs) const noexcept;

    // Validates and sets the interpolation boundaries.
    // Returns true if both lower and upper boundaries are valid for interpolation
    template <typename... PhotonGrid>
    bool set_boundaries(InterpState& state, size_t i, size_t j, size_t k, Real nu_obs,
                        PhotonGrid const&... photons) noexcept;
};

/********************************************************************************************************************
 * TEMPLATE METHOD: LogScaleInterp::set_boundaries
 * DESCRIPTION: Attempts to set the lower and upper boundary values for logarithmic interpolation.
 *              It updates the internal boundary members for:
 *                - Logarithmic observation time (t_obs_lo, log_t_ratio)
 *                - Logarithmic luminosity (L_lo, L_hi)
 *              The boundaries are set using data from the provided grids and the photon grids (via parameter pack).
 *              Returns true if both lower and upper boundaries are finite such that interpolation can proceed.
 ********************************************************************************************************************/
template <typename... PhotonGrid>
bool Observer::set_boundaries(InterpState& state, size_t i, size_t j, size_t k, Real nu_obs,
                              PhotonGrid const&... photons) noexcept {
    Real log_t_ratio = fast_log2(time(i, j, k + 1) / time(i, j, k));

    if (!std::isfinite(log_t_ratio)) [[unlikely]] {
        return false;
    }

    size_t eff_i = i * jet_3d;
    Real log_S_ratio = fast_log2(surface(i, j, k + 1) / surface(i, j, k));

    // continuing from previous boundary, shift the high boundary to lower.
    // Calling .I_nu()/.log_I_nu() could be expensive.
    if (state.last_hi != 0 && k == state.last_hi) {
        state.log_I_lo = state.log_I_hi;
    } else {
        Real nu = one_plus_z * nu_obs / doppler(i, j, k);
        state.log_I_lo = (photons(eff_i, j, k).compute_log2_I_nu(nu) + ...);
    }

    Real nu = one_plus_z * nu_obs / doppler(i, j, k + 1);
    state.log_I_hi = (photons(eff_i, j, k + 1).compute_log2_I_nu(nu) + ...);

    state.slope = (state.log_I_hi - state.log_I_lo + log_S_ratio) / log_t_ratio;

    if (!std::isfinite(state.slope)) [[unlikely]] {
        return false;
    }

    state.last_hi = k + 1;
    return true;
}

/********************************************************************************************************************
 * FUNCTION: iterate_to
 * DESCRIPTION: Helper function that advances an iterator (it) through an array (arr) until the value at that
 *              position exceeds the target value or the end of the array is reached. This is used for efficiently
 *              finding the appropriate position in a sorted array without binary search.
 ********************************************************************************************************************/
inline void iterate_to(Real value, Array const& arr, size_t& it) noexcept {
    while (it < arr.size() && arr(it) < value) {
        it++;
    }
}

/********************************************************************************************************************
 * TEMPLATE METHOD: Observer::specific_flux (multi-frequency overload)
 * DESCRIPTION: Returns the specific flux (as a MeshGrid) for multiple observed frequencies (nu_obs) by computing
 *              the specific flux for each frequency and assembling the results into a grid. This method accounts for
 *              relativistic beaming and cosmological effects.
 ********************************************************************************************************************/
template <typename... PhotonGrid>
MeshGrid Observer::specific_flux(Array const& t_obs, Array const& nu_obs, PhotonGrid const&... photons) {
    size_t t_obs_size = t_obs.size();
    size_t nu_size = nu_obs.size();

    MeshGrid F_nu({nu_size, t_obs_size}, 0);

    InterpState state;
    for (size_t l = 0; l < nu_size; l++) {
        // Loop over effective phi and theta grid points.
        for (size_t i = 0; i < eff_phi_size; i++) {
            for (size_t j = 0; j < theta_size; j++) {
                // Skip observation times that are below the grid's start time
                size_t t_idx = 0;
                iterate_to(time(i, j, 0), t_obs, t_idx);

                // Interpolate for observation times within the grid.
                for (size_t k = 0; k < t_size - 1 && t_idx < t_obs_size; k++) {
                    Real const t_hi = time(i, j, k + 1);

                    if (t_hi < t_obs(t_idx)) {
                        continue;
                    } else {
                        if (set_boundaries(state, i, j, k, nu_obs[l], photons...)) {
                            for (; t_idx < t_obs_size && t_obs(t_idx) < t_hi; t_idx++) {
                                F_nu(l, t_idx) += interpolate(state, i, j, k, t_obs(t_idx));
                            }
                        } else {
                            iterate_to(t_hi, t_obs, t_idx);
                        }
                    }
                }
            }
        }
    }

    // Normalize the flux by the factor (1+z)/(lumi_dist^2).
    Real const coef = one_plus_z / (lumi_dist * lumi_dist);
    F_nu *= coef;

    return F_nu;
}

/********************************************************************************************************************
 * TEMPLATE METHOD: Observer::specific_flux (single-frequency overload)
 * DESCRIPTION: Returns the specific flux (as an Array) for a single observed frequency (nu_obs) by computing the
 *              specific flux over the observation times.
 ********************************************************************************************************************/
template <typename... PhotonGrid>
Array Observer::specific_flux(Array const& t_obs, Real nu_obs, PhotonGrid const&... photons) {
    return xt::view(specific_flux(t_obs, Array({nu_obs}), photons...), 0);
}

/********************************************************************************************************************
 * TEMPLATE METHOD: Observer::spectrum (single-time overload)
 * DESCRIPTION: Returns the spectrum (as an Array) at a single observation time by computing the specific flux
 *              for each frequency in the given array.
 ********************************************************************************************************************/
template <typename... PhotonGrid>
Array Observer::spectrum(Array const& freqs, Real t_obs, PhotonGrid const&... photons) {
    return xt::view(spectra(freqs, Array({t_obs}), photons...), 0);
}

/********************************************************************************************************************
 * TEMPLATE METHOD: Observer::spectra (multi-time overload)
 * DESCRIPTION: Returns the spectra (as a MeshGrid) for multiple observation times by computing the specific flux
 *              for each frequency and transposing the result to get freq x time format.
 ********************************************************************************************************************/
template <typename... PhotonGrid>
MeshGrid Observer::spectra(Array const& freqs, Array const& t_obs, PhotonGrid const&... photons) {
    return xt::transpose(specific_flux(t_obs, freqs, photons...));
}

/********************************************************************************************************************
 * TEMPLATE METHOD: Observer::flux
 * DESCRIPTION: Computes the integrated flux over a frequency band specified by band_freq.
 *              It converts band boundaries to center frequencies, computes the specific flux at each frequency,
 *              and integrates (sums) the flux contributions weighted by the frequency bin widths.
 ********************************************************************************************************************/
template <typename... PhotonGrid>
Array Observer::flux(Array const& t_obs, Array const& band_freq, PhotonGrid const&... photons) {
    Array nu_obs = boundary_to_center_log(band_freq);
    MeshGrid F_nu = specific_flux(t_obs, nu_obs, photons...);
    Array flux({t_obs.size()}, 0);
    for (size_t i = 0; i < nu_obs.size(); ++i) {
        Real dnu = band_freq(i + 1) - band_freq(i);
        for (size_t j = 0; j < flux.size(); ++j) {
            flux(j) += dnu * F_nu(i, j);
        }
    }
    return flux;
}
