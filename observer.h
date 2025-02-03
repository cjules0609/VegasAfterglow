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

class LogScaleInterp {
   public:
    double z{0};
    size_t jet_3d{0};

    double interpRadius(double log_t) const;
    double interpIntensity(double log_t) const;
    double interpDoppler(double log_t) const;

    template <typename... PhotonGrid>
    bool trySetBoundary(size_t i, size_t j, size_t k, Array const& log_r, MeshGrid3d const& t_obs,
                        MeshGrid3d const& doppler, double nu_obs, PhotonGrid const&... photons);

   private:
    double log_r_lo{0};
    double log_r_hi{0};

    double log_t_lo{0};
    double log_t_hi{0};

    double log_d_lo{0};
    double log_d_hi{0};

    double log_I_lo{0};
    double log_I_hi{0};

    size_t idx_hi{0};
};

class Observer {
   public:
    Observer(Coord const& coord);
    Observer() = delete;

    MeshGrid3d t_obs_grid;
    MeshGrid3d doppler;
    double theta_obs{0};
    double lumi_dist{1};
    double z{0};

    template <typename Dynamics>
    void observe(Dynamics const& dyn, double theta_obs, double lumi_dist, double z);

    template <typename... PhotonGrid>
    Array specificFlux(Array const& t_obs, double nu_obs, PhotonGrid const&... photons);

    template <typename... PhotonGrid>
    MeshGrid specificFlux(Array const& t_obs, Array const& nu_obs, PhotonGrid const&... photons);

    template <typename... PhotonGrid>
    Array flux(Array const& t_obs, Array const& band_freq, PhotonGrid const&... photons);

    template <typename... PhotonGrid>
    MeshGrid spectrum(Array const& t_obs, Array const& band_freq, PhotonGrid const&... photons);

   private:
    MeshGrid dOmega;
    Array log_r;
    LogScaleInterp interp;
    Coord const& coord;
    size_t eff_phi_size{1};

    double dopplerInterp(double log_t, size_t i, size_t j, size_t k) const;
    double radiusInterp(double log_t, size_t i, size_t j, size_t k) const;
    double intensityInterp(double log_t, size_t i, size_t j, size_t k, double log_I_lo, double log_I_hi) const;
    void calcObsTimeGrid(MeshGrid3d const& Gamma, MeshGrid3d const& t_eng);
    void calcSolidAngle();

    template <typename Iter, typename... PhotonGrid>
    void calcSpecificFlux(Iter f_nu, Array const& t_obs, double nu_obs, PhotonGrid const&... photons);
};

template <typename Dynamics>
void Observer::observe(Dynamics const& dyn, double theta_obs, double lumi_dist, double z) {
    auto [phi_size, theta_size, r_size] = dyn.shape();

    this->theta_obs = theta_obs;
    this->z = z;
    this->lumi_dist = lumi_dist;
    this->interp.z = z;
    this->interp.jet_3d = static_cast<size_t>((phi_size > 1));

    if (theta_obs == 0 && interp.jet_3d == 0) {
        eff_phi_size = 1;
    } else {
        eff_phi_size = coord.phi.size();
    }
    calcSolidAngle();
    calcObsTimeGrid(dyn.Gamma_rel, dyn.t_eng);
}

template <typename... PhotonGrid>
bool LogScaleInterp::trySetBoundary(size_t i, size_t j, size_t k_lo, Array const& log_r, MeshGrid3d const& t_obs,
                                    MeshGrid3d const& doppler, double nu_obs, PhotonGrid const&... photons) {
    if (idx_hi != 0 && k_lo == idx_hi) {  // just move one grid, lower->higher
        log_r_lo = log_r_hi;
        log_t_lo = log_t_hi;
        log_d_lo = log_d_hi;
        log_I_lo = log_I_hi;
    } else {
        log_r_lo = log_r[k_lo];
        log_t_lo = fastLog(t_obs[i][j][k_lo]);
        log_d_lo = fastLog(doppler[i][j][k_lo]);

        double D = doppler[i][j][k_lo];
        double nu = (1 + z) * nu_obs / D;

        log_I_lo = fastLog((photons[i * jet_3d][j][k_lo].I_nu(nu) + ...));
    }
    log_r_hi = log_r[k_lo + 1];
    log_t_hi = fastLog(t_obs[i][j][k_lo + 1]);
    log_d_hi = fastLog(doppler[i][j][k_lo + 1]);

    double D = doppler[i][j][k_lo + 1];
    double nu = (1 + z) * nu_obs / D;

    log_I_hi = fastLog((photons[i * jet_3d][j][k_lo + 1].I_nu(nu) + ...));

    idx_hi = k_lo + 1;
    return std::isfinite(log_t_lo) && std::isfinite(log_t_hi);
}

template <typename Iter, typename... PhotonGrid>
void Observer::calcSpecificFlux(Iter f_nu, Array const& t_obs, double nu_obs, const PhotonGrid&... photons) {
    auto [phi_size, theta_size, r_size] = coord.shape();

    size_t t_size = t_obs.size();

    auto update_flux = [&](size_t id, double solid_angle) {
        double const log_t = fastLog(t_obs[id]);
        double const D = interp.interpDoppler(log_t);
        double const r = interp.interpRadius(log_t);
        double const I_nu = interp.interpIntensity(log_t);
        f_nu[id] += D * D * D * I_nu * r * r * solid_angle;
    };

    for (size_t i = 0; i < eff_phi_size; i++) {
        for (size_t j = 0; j < theta_size; j++) {
            double const solid_angle = dOmega[i][j];
            if (!interp.trySetBoundary(i, j, 0, log_r, t_obs_grid, doppler, nu_obs, photons...)) {
                continue;
            }

            size_t t_idx = 0;
            // extrapolation for EAT outside the grid
            for (; t_idx < t_size && t_obs[t_idx] < t_obs_grid[i][j][0]; t_idx++) {
#ifdef EXTRAPOLATE
                update_flux(t_idx, solid_angle);
#endif
            }

            // interpolation for EAT inside the grid
            for (size_t k = 0; k < r_size - 1 && t_idx < t_size; k++) {
                double const t_obs_lo = t_obs_grid[i][j][k];
                double const t_obs_hi = t_obs_grid[i][j][k + 1];

                if (t_obs_lo <= t_obs[t_idx] && t_obs[t_idx] < t_obs_hi) {
                    if (!interp.trySetBoundary(i, j, k, log_r, t_obs_grid, doppler, nu_obs, photons...)) {
                        continue;
                    }
                }

                for (; t_idx < t_size && t_obs_lo <= t_obs[t_idx] && t_obs[t_idx] < t_obs_hi; t_idx++) {
                    update_flux(t_idx, solid_angle);
                }
            }
#ifdef EXTRAPOLATE
            // extrapolation for EAT outside the grid
            for (; t_idx < t_size; t_idx++) {
                update_flux(t_idx, solid_angle);
            }
#endif
        }
    }

    // normalize the flux
    double const coef = (1 + z) / (lumi_dist * lumi_dist);
    for (size_t i = 0; i < t_size; ++i) {
        f_nu[i] *= coef;
    }
}

template <typename... PhotonGrid>
Array Observer::specificFlux(Array const& t_obs, double nu_obs, PhotonGrid const&... photons) {
    Array F_nu = zeros(t_obs.size());
    calcSpecificFlux(F_nu.data(), t_obs, nu_obs, photons...);
    return F_nu;
}

template <typename... PhotonGrid>
MeshGrid Observer::specificFlux(Array const& t_obs, Array const& nu_obs, PhotonGrid const&... photons) {
    MeshGrid F_nu = createGrid(nu_obs.size(), t_obs.size(), 0);
    size_t t_num = t_obs.size();
    for (size_t l = 0; l < nu_obs.size(); ++l) {
        calcSpecificFlux(F_nu.data() + l * t_num, t_obs, nu_obs[l], photons...);
    }
    return F_nu;
}

template <typename... PhotonGrid>
Array Observer::flux(Array const& t_obs, Array const& band_freq, PhotonGrid const&... photons) {
    Array nu_obs = boundaryToCenterLog(band_freq);
    MeshGrid F_nu = specificFlux(t_obs, nu_obs, photons...);
    Array flux = zeros(t_obs.size());
    for (size_t i = 0; i < F_nu.size(); ++i) {
        double dnu = band_freq[i + 1] - band_freq[i];
        for (size_t j = 0; j < flux.size(); ++j) {
            flux[j] += dnu * F_nu[i][j];
        }
    }
    return flux;
}
#endif