//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _OBSERVER_
#define _OBSERVER_

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

    // Interpolates the radius using the logarithm of the observation time (log_t)
    Real interpRadius(Real log_t) const;
    // Interpolates the intensity using the logarithm of the observation time (log_t)
    Real interpIntensity(Real log_t) const;
    // Interpolates the Doppler factor using the logarithm of the observation time (log_t)
    Real interpDoppler(Real log_t) const;

    // Tries to set the interpolation boundaries using the provided logarithmic radius array, observation time grid,
    // Doppler grid, observed frequency, and one or more photon grids.
    template <typename... PhotonGrid>
    bool trySetBoundary(size_t i, size_t j, size_t k, MeshGrid3d const& r, MeshGrid3d const& t_obs,
                        MeshGrid3d const& doppler, Real nu_obs, PhotonGrid const&... photons);

   private:
    Real log_r_lo{0};  // Lower boundary of logarithmic radius
    Real log_r_hi{0};  // Upper boundary of logarithmic radius

    Real log_t_lo{0};  // Lower boundary of logarithmic observation time
    Real log_t_hi{0};  // Upper boundary of logarithmic observation time

    Real log_d_lo{0};  // Lower boundary of logarithmic Doppler factor
    Real log_d_hi{0};  // Upper boundary of logarithmic Doppler factor

    Real log_I_lo{0};  // Lower boundary of logarithmic intensity
    Real log_I_hi{0};  // Upper boundary of logarithmic intensity

    size_t idx_hi{0};  // Index for the upper boundary in the grid
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
    Observer(Coord const& coord);
    Observer() = delete;  // Default constructor is deleted.

    MeshGrid3d t_obs_grid;  // Grid of observation times
    MeshGrid3d doppler;     // Grid of Doppler factors
    Real theta_obs{0};      // Observer's theta angle
    Real lumi_dist{1};      // Luminosity distance
    Real z{0};              // Redshift

    // Observes the provided dynamics (dyn) with the given observation parameters.
    template <typename Dynamics>
    void observe(Dynamics const& dyn, Real theta_obs, Real lumi_dist, Real z);

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
    MeshGrid spectrum(Array const& t_obs, Array const& band_freq, PhotonGrid const&... photons);

   private:
    MeshGrid dOmega;                    // Grid of solid angles
    MeshGrid3d const* r_grid{nullptr};  // Grid of radius in logarithmic scale
    LogScaleInterp interp;              // Log-scale interpolation helper
    Coord const& coord;                 // Reference to the coordinate object
    size_t eff_phi_size{1};             // Effective number of phi grid points

    // Calculates the observation time grid based on Gamma and engine time array.
    void calcObsTimeGrid(MeshGrid3d const& Gamma, MeshGrid3d const& r);
    // Calculates the solid angle grid.
    void calcSolidAngle();

    // Template helper method to compute specific flux and store the result in a provided iterator (f_nu).
    template <typename Iter, typename... PhotonGrid>
    void calcSpecificFlux(Iter f_nu, Array const& t_obs, Real nu_obs, PhotonGrid const&... photons);
};

/********************************************************************************************************************
 * TEMPLATE METHOD: Observer::observe
 * DESCRIPTION: Updates the Observer object using the provided dynamics (dyn). It sets the observer’s theta, redshift,
 *              and luminosity distance, and updates internal parameters such as the effective phi grid size.
 *              It then calculates the solid angle grid and observation time grid based on the provided dynamics.
 ********************************************************************************************************************/
template <typename Dynamics>
void Observer::observe(Dynamics const& dyn, Real theta_view, Real luminosity_dist, Real redshift) {
    // Extract grid dimensions from the dynamics object.
    auto [phi_size, theta_size, t_size] = dyn.shape();

    theta_obs = theta_view;
    lumi_dist = luminosity_dist;
    z = redshift;
    interp.z = redshift;
    // Determine if the jet is 3D (more than one phi value)
    interp.jet_3d = static_cast<size_t>((phi_size > 1));

    // Set effective phi grid size based on the observation angle and jet dimensionality.
    if (theta_view == 0 && interp.jet_3d == 0) {
        eff_phi_size = 1;
    } else {
        eff_phi_size = coord.phi.size();
    }
    // Calculate the solid angle grid and observation time grid.
    r_grid = &(dyn.r);
    calcSolidAngle();
    calcObsTimeGrid(dyn.Gamma_rel, dyn.r);
}

/********************************************************************************************************************
 * TEMPLATE METHOD: LogScaleInterp::trySetBoundary
 * DESCRIPTION: Attempts to set the lower and upper boundary values for logarithmic interpolation.
 *              It updates the internal boundary members for:
 *                - Logarithmic radius (log_r_lo, log_r_hi)
 *                - Logarithmic observation time (log_t_lo, log_t_hi)
 *                - Logarithmic Doppler factor (log_d_lo, log_d_hi)
 *                - Logarithmic intensity (log_I_lo, log_I_hi)
 *              The boundaries are set using data from the provided grids and the photon grids (via parameter pack).
 *              Returns true if both lower and upper log observation time boundaries are finite.
 ********************************************************************************************************************/
template <typename... PhotonGrid>
bool LogScaleInterp::trySetBoundary(size_t i, size_t j, size_t k_lo, MeshGrid3d const& r, MeshGrid3d const& t_obs,
                                    MeshGrid3d const& doppler, Real nu_obs, PhotonGrid const&... photons) {
    // If continuing from previous boundary, shift the high boundary values to lower. Calling .I_nu() is expensive.
    if (idx_hi != 0 && k_lo == idx_hi) {
        log_r_lo = log_r_hi;
        log_t_lo = log_t_hi;
        log_d_lo = log_d_hi;
        log_I_lo = log_I_hi;
    } else {
        // Otherwise, set lower boundary values using the k_lo index.
        log_r_lo = fastLog(r[i * jet_3d][j][k_lo]);
        log_t_lo = fastLog(t_obs[i][j][k_lo]);
        log_d_lo = fastLog(doppler[i][j][k_lo]);

        Real D = doppler[i][j][k_lo];
        Real nu = (1 + z) * nu_obs / D;

        log_I_lo = fastLog((photons[i * jet_3d][j][k_lo].I_nu(nu) + ...));
    }
    // Set upper boundary values at index k_lo + 1.
    log_r_hi = fastLog(r[i * jet_3d][j][k_lo + 1]);
    log_t_hi = fastLog(t_obs[i][j][k_lo + 1]);
    log_d_hi = fastLog(doppler[i][j][k_lo + 1]);

    Real D = doppler[i][j][k_lo + 1];
    Real nu = (1 + z) * nu_obs / D;

    log_I_hi = fastLog((photons[i * jet_3d][j][k_lo + 1].I_nu(nu) + ...));

    idx_hi = k_lo + 1;
    return std::isfinite(log_t_lo) && std::isfinite(log_t_hi);
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
template <typename Iter, typename... PhotonGrid>
void Observer::calcSpecificFlux(Iter f_nu, Array const& t_obs, Real nu_obs, const PhotonGrid&... photons) {
    auto [phi_size, theta_size, t_size] = coord.shape();
    size_t t_obs_size = t_obs.size();

    // Define a lambda to update the flux for a given observation time index.
    auto update_flux = [&](size_t id, Real solid_angle) {
        Real const log_t = fastLog(t_obs[id]);
        Real const D = interp.interpDoppler(log_t);
        Real const r = interp.interpRadius(log_t);
        Real const I_nu = interp.interpIntensity(log_t);
        // The flux contribution scales as D^3 * I_nu * r^2 * solid_angle.
        f_nu[id] += D * D * D * I_nu * r * r * solid_angle;
    };

    // Loop over effective phi and theta grid points.
    for (size_t i = 0; i < eff_phi_size; i++) {
        for (size_t j = 0; j < theta_size; j++) {
            Real const solid_angle = dOmega[i][j];
            // Attempt to set the initial boundary values; if unsuccessful, skip this grid cell.
            if (!interp.trySetBoundary(i, j, 0, *r_grid, t_obs_grid, doppler, nu_obs, photons...)) {
                continue;
            }

            size_t t_idx = 0;
            // Extrapolation for observation times below the grid (if enabled). Otherwise, skip to the next required
            // cell untill the first observation time is reached.
            for (; t_idx < t_obs_size && t_obs[t_idx] < t_obs_grid[i][j][0]; t_idx++) {
#ifdef EXTRAPOLATE
                update_flux(t_idx, solid_angle);
#endif
            }

            // Interpolate for observation times within the grid.
            for (size_t k = 0; k < t_size - 1 && t_idx < t_obs_size; k++) {
                Real const t_obs_lo = t_obs_grid[i][j][k];
                Real const t_obs_hi = t_obs_grid[i][j][k + 1];

                if (t_obs_lo <= t_obs[t_idx] && t_obs[t_idx] < t_obs_hi) {
                    if (!interp.trySetBoundary(i, j, k, *r_grid, t_obs_grid, doppler, nu_obs, photons...)) {
                        continue;
                    }
                }

                for (; t_idx < t_obs_size && t_obs_lo <= t_obs[t_idx] && t_obs[t_idx] < t_obs_hi; t_idx++) {
                    update_flux(t_idx, solid_angle);
                }
            }
#ifdef EXTRAPOLATE
            //   Extrapolation for observation times above the grid.
            for (; t_idx < t_obs_size; t_idx++) {
                print("extrapolating");
                update_flux(t_idx, solid_angle);
            }
#endif
        }
    }

    // Normalize the flux by the factor (1+z)/(lumi_dist^2).
    Real const coef = (1 + z) / (lumi_dist * lumi_dist);
    for (size_t i = 0; i < t_obs_size; ++i) {
        f_nu[i] *= coef;
    }
}

/********************************************************************************************************************
 * TEMPLATE METHOD: Observer::specificFlux (single-frequency overload)
 * DESCRIPTION: Returns the specific flux (as an Array) for a single observed frequency (nu_obs) by computing the
 *              specific flux over the observation times.
 ********************************************************************************************************************/
template <typename... PhotonGrid>
Array Observer::specificFlux(Array const& t_obs, Real nu_obs, PhotonGrid const&... photons) {
    if (r_grid == nullptr) {
        throw std::runtime_error("Observer::specificFlux: r_grid is not set. Call observe() first.");
    }
    Array F_nu = zeros(t_obs.size());
    calcSpecificFlux(F_nu.data(), t_obs, nu_obs, photons...);
    return F_nu;
}

/********************************************************************************************************************
 * TEMPLATE METHOD: Observer::specificFlux (multi-frequency overload)
 * DESCRIPTION: Returns the specific flux (as a MeshGrid) for multiple observed frequencies (nu_obs) by computing
 *              the specific flux for each frequency and assembling the results into a grid.
 ********************************************************************************************************************/
template <typename... PhotonGrid>
MeshGrid Observer::specificFlux(Array const& t_obs, Array const& nu_obs, PhotonGrid const&... photons) {
    if (r_grid == nullptr) {
        throw std::runtime_error("Observer::specificFlux: r_grid is not set. Call observe() first.");
    }
    MeshGrid F_nu = createGrid(nu_obs.size(), t_obs.size(), 0);
    size_t t_num = t_obs.size();
    for (size_t l = 0; l < nu_obs.size(); ++l) {
        calcSpecificFlux(F_nu.data() + l * t_num, t_obs, nu_obs[l], photons...);
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
    if (r_grid == nullptr) {
        throw std::runtime_error("Observer::specificFlux: r_grid is not set. Call observe() first.");
    }
    Array nu_obs = boundaryToCenterLog(band_freq);
    MeshGrid F_nu = specificFlux(t_obs, nu_obs, photons...);
    Array flux = zeros(t_obs.size());
    for (size_t i = 0; i < F_nu.size(); ++i) {
        Real dnu = band_freq[i + 1] - band_freq[i];
        for (size_t j = 0; j < flux.size(); ++j) {
            flux[j] += dnu * F_nu[i][j];
        }
    }
    return flux;
}
#endif