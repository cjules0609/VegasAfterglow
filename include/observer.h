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
 * DESCRIPTION: Provides logarithmic scale interpolation for various shock-related quantities. It includes methods
 *              to interpolate radius, intensity, and Doppler factor based on the logarithm of the observation time.
 *              It also defines a method to attempt to set boundary conditions for the interpolation using provided
 *              grids.
 ********************************************************************************************************************/
class LogScaleInterp {
   public:
    Real z{0};         // Redshift
    size_t jet_3d{0};  // Flag indicating if the jet is non-axis-symmetric (non-zero if true)

    // Interpolates the radius, Intensity, and Doppler factor using  the observation time (t)
    Real interpLuminosity(Real t_obs) const;

    // Tries to set the interpolation boundaries using the provided logarithmic radius array, observation time grid,
    // Doppler grid, observed frequency, and one or more photon grids.
    template <typename... PhotonGrid>
    bool validateInterpBoundary(size_t i, size_t j, size_t k, MeshGrid3d const& dOmega, MeshGrid3d const& r_grid,
                                MeshGrid3d const& t_obs, MeshGrid3d const& doppler, Real nu_obs,
                                PhotonGrid const&... photons);

   private:
    Real log_t_ratio{0};  // Ratio of logarithmic observation time
    Real log_L_ratio{0};

    Real t_obs_lo{0};  // Lower boundary of observation time

    Real L_lo{0};  // Lower boundary of luminosity
    Real L_hi{0};  // Upper boundary of luminosity

    size_t idx_hi{0};  // Index for the upper boundary in the grid
    size_t idx_i{0};
    size_t idx_j{0};
};

/********************************************************************************************************************
 * CLASS: Observer
 * DESCRIPTION: Represents an observer in the shock simulation. The Observer stores observation grids (for
 *              observation time and Doppler factors), as well as parameters such as the observation angle,
 *              luminosity distance, and redshift. It provides methods for "observing" the simulation dynamics
 *              and for computing the specific flux, integrated flux, and spectrum.
 ********************************************************************************************************************/
class Observer {
   public:
    // Constructor: Requires a coordinate reference to initialize the observer.
    template <typename Dynamics>
    Observer(Coord const& coord, Dynamics const& dyn, Real theta_view, Real luminosity_dist, Real redshift);
    Observer() = delete;    // Default constructor is deleted.
    MeshGrid3d t_obs_grid;  // Grid of observation times
    MeshGrid3d doppler;     // Grid of Doppler factors
    Real theta_obs{0};      // Observer's theta angle
    Real lumi_dist{1};      // Luminosity distance
    Real z{0};              // Redshift

    // Observes the provided dynamics (dyn) with the given observation parameters.
    void changeViewingAngle(Real theta_obs);

    // Computes the specific flux at a single observed frequency (nu_obs) for the given observation times,
    // using one or more photon grids.
    template <typename... PhotonGrid>
    Array specificFlux(Array const& t_obs, Real nu_obs, PhotonGrid const&... photons);

    // Computes the specific flux at multiple observed frequencies (nu_obs) for the given observation times.
    template <typename... PhotonGrid>
    MeshGrid specificFlux(Array const& t_obs, Array const& nu_obs, PhotonGrid const&... photons);

    // Computes the integrated flux over a frequency band defined by band_freq.
    template <typename... PhotonGrid>
    Array flux(Array const& t_obs, Array const& band_freq, PhotonGrid const&... photons);

    // Computes the spectrum over a frequency band.
    template <typename... PhotonGrid>
    MeshGrid spectrum(Array const& freqs, Array const& t_obs, PhotonGrid const&... photons);

    template <typename... PhotonGrid>
    Array spectrum(Array const& freqs, Real t_obs, PhotonGrid const&... photons);

   private:
    LogScaleInterp interp;         // Log-scale interpolation helper
    MeshGrid3d dOmega;             // Grid of solid angles
    MeshGrid3d const& r_grid;      // Grid of radius
    MeshGrid3d const& theta_grid;  // Grid of theta
    MeshGrid3d const& Gamma;       // Grid of Lorentz factors
    Coord const& coord;            // Reference to the coordinate object
    size_t eff_phi_size{1};        // Effective number of phi grid points

    // Calculates the observation time grid based on Gamma and engine time array.
    void calcObsTimeGrid();
    // Calculates the solid angle grid.
    void calcSolidAngle();

    // Template helper method to compute specific flux and store the result in a provided iterator (f_nu).
    template <typename View, typename... PhotonGrid>
    void calcSpecificFlux(View& f_nu, Array const& t_obs, Real nu_obs, PhotonGrid const&... photons);

    template <typename View, typename... PhotonGrid>
    void calcSpectrum(View& f_nu, Array const& nu_obs, Real t_obs, PhotonGrid const&... photons);
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
bool LogScaleInterp::validateInterpBoundary(size_t i, size_t j, size_t k_lo, MeshGrid3d const& dOmega,
                                            MeshGrid3d const& r_grid, MeshGrid3d const& t_obs,
                                            MeshGrid3d const& doppler, Real nu_obs, PhotonGrid const&... photons) {
    t_obs_lo = t_obs(i, j, k_lo);
    log_t_ratio = fastLog(t_obs(i, j, k_lo + 1) / t_obs(i, j, k_lo));

    if (!std::isfinite(log_t_ratio) || log_t_ratio == 0) {
        return false;
    }
    size_t phi_idx = i * jet_3d;

    // If continuing from previous boundary, shift the high boundary values to lower. Calling .I_nu() is expensive.
    if (idx_hi != 0 && k_lo == idx_hi && idx_i == i && idx_j == j) {
        L_lo = L_hi;
    } else {
        Real D = doppler(i, j, k_lo);
        Real nu = (1 + z) * nu_obs / D;
        Real I_lo = (photons(phi_idx, j, k_lo).I_nu(nu) + ...);
        Real r = r_grid(phi_idx, j, k_lo);
        Real solid_angle = dOmega(i, j, k_lo);
        L_lo = I_lo * r * r * D * D * D * solid_angle;
    }
    Real D = doppler(i, j, k_lo + 1);
    Real nu = (1 + z) * nu_obs / D;
    Real I_hi = (photons(phi_idx, j, k_lo + 1).I_nu(nu) + ...);
    Real r = r_grid(phi_idx, j, k_lo + 1);
    Real solid_angle = dOmega(i, j, k_lo + 1);
    L_hi = I_hi * r * r * D * D * D * solid_angle;

    log_L_ratio = fastLog(L_hi / L_lo);
    if (!std::isfinite(log_L_ratio)) {
        return false;
    }
    idx_hi = k_lo + 1;
    idx_i = i;
    idx_j = j;
    return true;
}

/********************************************************************************************************************
 * CONSTRUCTOR: Observer::Observer
 * DESCRIPTION: Constructs an Observer object with the given coordinate grid, shocks, observation angle, luminosity
 *              and redshift. It initializes the observation time and Doppler factor grids, as well as the
 *              interpolation object.
 ********************************************************************************************************************/
template <typename Dynamics>
Observer::Observer(Coord const& coord, Dynamics const& shock, Real theta_view, Real luminosity_dist, Real redshift)
    : t_obs_grid({coord.phi.size(), coord.theta.size(), coord.t.size()}, 0),
      doppler({coord.phi.size(), coord.theta.size(), coord.t.size()}, 1),
      theta_obs(theta_view),
      lumi_dist(luminosity_dist),
      z(redshift),
      interp(),
      dOmega({coord.phi.size(), coord.theta.size(), coord.t.size()}, 0),
      r_grid(shock.r),
      theta_grid(shock.theta),
      Gamma(shock.Gamma),
      coord(coord),
      eff_phi_size(1) {
    auto [phi_size, theta_size, t_size] = shock.shape();
    interp.z = redshift;
    // Determine if the jet is 3D (more than one phi value)
    interp.jet_3d = static_cast<size_t>((phi_size > 1));
    // Set effective phi grid size based on the observation angle and jet dimensionality.
    if (theta_view == 0 && interp.jet_3d == 0) {
        eff_phi_size = 1;  // optimize for on-axis observer
    } else {
        eff_phi_size = coord.phi.size();
    }
    // Calculate the solid angle grid and observation time grid.
    calcSolidAngle();
    calcObsTimeGrid();
}

/********************************************************************************************************************
 * TEMPLATE METHOD: Observer::calcSpecificFlux
 * DESCRIPTION: Calculates the specific flux at each observation time for a given observed frequency nu_obs.
 *              For each effective phi and theta grid point, it:
 *                - Interpolates the Doppler factor, radius, and intensity using the LogScaleInterp object.
 *                - Updates the flux array by accumulating contributions from each grid cell multiplied by its
 *                  solid angle.
 *              Finally, the flux is normalized by the factor (1+z)/(lumi_dist^2).
 ********************************************************************************************************************/
template <typename View, typename... PhotonGrid>
void Observer::calcSpecificFlux(View& f_nu, Array const& t_obs, Real nu_obs, const PhotonGrid&... photons) {
    auto [ignore, theta_size, t_size] = coord.shape();
    size_t t_obs_size = t_obs.size();

    // Loop over effective phi and theta grid points.
    for (size_t i = 0; i < eff_phi_size; i++) {
        for (size_t j = 0; j < theta_size; j++) {
            // Skip observation times that are below the grid's start time
            size_t t_idx = 0;
            while (t_idx < t_obs_size && t_obs(t_idx) < t_obs_grid(i, j, 0)) {
                t_idx++;
            }

            // Interpolate for observation times within the grid.
            for (size_t k = 0; k < t_size - 1 && t_idx < t_obs_size; k++) {
                Real const t_lo = t_obs_grid(i, j, k);
                Real const t_hi = t_obs_grid(i, j, k + 1);

                if (t_lo <= t_obs(t_idx) && t_obs(t_idx) < t_hi) {
                    if (!interp.validateInterpBoundary(i, j, k, dOmega, r_grid, t_obs_grid, doppler, nu_obs,
                                                       photons...)) {
                        for (; t_idx < t_obs_size && t_obs(t_idx) < t_hi; t_idx++) {
                        }
                        continue;
                    }
                }

                for (; t_idx < t_obs_size && t_lo <= t_obs(t_idx) && t_obs(t_idx) < t_hi; t_idx++) {
                    f_nu(t_idx) += interp.interpLuminosity(t_obs(t_idx));
                }
            }
        }
    }

    // Normalize the flux by the factor (1+z)/(lumi_dist^2).
    Real const coef = (1 + z) / (lumi_dist * lumi_dist);
    for (size_t i = 0; i < t_obs_size; ++i) {
        f_nu(i) *= coef;
    }
}

/********************************************************************************************************************
 * TEMPLATE METHOD: Observer::calcSpectrum
 * DESCRIPTION: Calculates the specific flux at each observation frequency for a given time t_obs.
 *              For each effective phi and theta grid point, it:
 *                - Interpolates the Doppler factor, radius, and intensity using the LogScaleInterp object.
 *                - Updates the flux array by accumulating contributions from each grid cell multiplied by its
 *                  solid angle.
 *              Finally, the flux is normalized by the factor (1+z)/(lumi_dist^2).
 ********************************************************************************************************************/
template <typename View, typename... PhotonGrid>
void Observer::calcSpectrum(View& f_nu, Array const& nu_obs, Real t_obs, const PhotonGrid&... photons) {
    auto [ignore, theta_size, t_size] = coord.shape();
    size_t nu_size = nu_obs.size();

    // Loop over effective phi and theta grid points.
    for (size_t i = 0; i < eff_phi_size; i++) {
        for (size_t j = 0; j < theta_size; j++) {
            for (size_t k = 0; k < t_size - 1; k++) {
                Real const t_lo = t_obs_grid(i, j, k);
                Real const t_hi = t_obs_grid(i, j, k + 1);

                if (t_lo <= t_obs && t_obs < t_hi) {
                    for (size_t nu_idx = 0; nu_idx < nu_size; ++nu_idx) {
                        if (interp.validateInterpBoundary(i, j, k, dOmega, r_grid, t_obs_grid, doppler, nu_obs[nu_idx],
                                                          photons...)) {
                            f_nu(nu_idx) += interp.interpLuminosity(t_obs);
                        }
                    }
                    break;
                }
            }
        }
    }

    // Normalize the flux by the factor (1+z)/(lumi_dist^2).
    Real const coef = (1 + z) / (lumi_dist * lumi_dist);
    for (size_t i = 0; i < nu_size; ++i) {
        f_nu(i) *= coef;
    }
}

/********************************************************************************************************************
 * TEMPLATE METHOD: Observer::specificFlux (single-frequency overload)
 * DESCRIPTION: Returns the specific flux (as an Array) for a single observed frequency (nu_obs) by computing the
 *              specific flux over the observation times.
 ********************************************************************************************************************/
template <typename... PhotonGrid>
Array Observer::specificFlux(Array const& t_obs, Real nu_obs, PhotonGrid const&... photons) {
    Array F_nu({t_obs.size()}, 0);
    auto view = xt::view(F_nu, xt::all());
    calcSpecificFlux(view, t_obs, nu_obs, photons...);
    return F_nu;
}

/********************************************************************************************************************
 * TEMPLATE METHOD: Observer::specificFlux (multi-frequency overload)
 * DESCRIPTION: Returns the specific flux (as a MeshGrid) for multiple observed frequencies (nu_obs) by computing
 *              the specific flux for each frequency and assembling the results into a grid.
 ********************************************************************************************************************/
template <typename... PhotonGrid>
MeshGrid Observer::specificFlux(Array const& t_obs, Array const& nu_obs, PhotonGrid const&... photons) {
    MeshGrid F_nu({nu_obs.size(), t_obs.size()}, 0);
    for (size_t l = 0; l < nu_obs.size(); ++l) {
        auto view = xt::view(F_nu, l, xt::all());
        calcSpecificFlux(view, t_obs, nu_obs(l), photons...);
    }
    return F_nu;
}

template <typename... PhotonGrid>
Array Observer::spectrum(Array const& freqs, Real t_obs, PhotonGrid const&... photons) {
    Array F_nu({freqs.size()}, 0);
    auto view = xt::view(F_nu, xt::all());
    calcSpectrum(view, freqs, t_obs, photons...);
    return F_nu;
}

template <typename... PhotonGrid>
MeshGrid Observer::spectrum(Array const& freqs, Array const& t_obs, PhotonGrid const&... photons) {
    MeshGrid F_nu({t_obs.size(), freqs.size()}, 0);
    for (size_t l = 0; l < t_obs.size(); ++l) {
        auto view = xt::view(F_nu, l, xt::all());
        calcSpectrum(view, freqs, t_obs(l), photons...);
    }
    return F_nu;
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
