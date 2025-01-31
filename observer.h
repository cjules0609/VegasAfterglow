#ifndef _OBSERVER_
#define _OBSERVER_

#include <iostream>
#include <stdexcept>
#include <thread>
#include <vector>

#include "afterglow.h"
#include "macros.h"
#include "mesh.h"
class Observer {
   public:
    Observer(Coord const& coord);
    Observer() = delete;

    MeshGrid3d t_obs_grid;
    MeshGrid3d doppler;
    Array dphi;
    double theta_obs{0};
    double lumi_dist{1};
    double z{0};

    template <typename Dynamics>
    void observe(Dynamics const& dyn, double theta_obs, double lumi_dist, double z);

    template <typename... PhotonMesh>
    Array specificFlux(Array const& t_obs, double nu_obs, PhotonMesh const&... photons);

    template <typename... PhotonMesh>
    MeshGrid specificFlux(Array const& t_obs, Array const& nu_obs, PhotonMesh const&... photons);

    template <typename... PhotonMesh>
    Array flux(Array const& t_obs, Array const& band_pass_freq, PhotonMesh const&... photons);

    template <typename... PhotonMesh>
    MeshGrid spectrum(Array const& t_obs, Array const& band_pass_freq, PhotonMesh const&... photons);

   private:
    MeshGrid3d log_t_obs;    // for log scale interpolation
    MeshGrid3d log_doppler;  // for log scale interpolation
    Array log_r;             // for log scale interpolation
    Coord const& coord;
    size_t effective_phi_size{1};

    double dopplerInterp(double log_t, size_t i, size_t j, size_t k) const;
    double radiusInterp(double log_t, size_t i, size_t j, size_t k) const;
    double intensityInterp(double log_t, size_t i, size_t j, size_t k, double log_I_lo, double log_I_hi) const;
    void calcObsTimeGrid(MeshGrid const& Gamma, MeshGrid const& t_eng);

    template <typename Iter, typename... PhotonMesh>
    void calcSpecificFlux(Iter f_nu, Array const& t_obs, double nu_obs, PhotonMesh const&... photons) const;

    template <typename... PhotonMesh>
    double getSpecificIntensity(size_t i, size_t j, size_t k, double nu, PhotonMesh const&... photons) const;
};

#include "observer.tpp"

template <typename Dynamics>
void Observer::observe(Dynamics const& dyn, double theta_obs, double lumi_dist, double z) {
    this->theta_obs = theta_obs;
    this->z = z;
    this->lumi_dist = lumi_dist;
    if (theta_obs == 0) {
        effective_phi_size = 1;
        std::fill(dphi.begin(), dphi.end(), 2 * con::pi);
    } else {
        effective_phi_size = coord.phi.size();
        for (size_t i = 0; i < coord.phi.size(); ++i) {
            dphi[i] = coord.phi_b[i + 1] - coord.phi_b[i];
        }
    }
    calcObsTimeGrid(dyn.Gamma, dyn.t_eng);
}

template <typename... PhotonMesh>
inline double Observer::getSpecificIntensity(size_t i, size_t j, size_t k, double nu_obs,
                                             PhotonMesh const&... photons) const {
    double D = doppler[i][j][k];
    double nu = (1 + z) * nu_obs / D;
    return (photons[j][k].I_nu(nu) + ...);
}

template <typename Iter, typename... PhotonMesh>
void Observer::calcSpecificFlux(Iter f_nu, Array const& t_obs, double nu_obs, const PhotonMesh&... photons) const {
    size_t const theta_size = coord.theta.size();
    size_t const r_size = coord.r.size();
    size_t const t_size = t_obs.size();

    for (size_t i = 0; i < effective_phi_size; i++) {
        for (size_t j = 0; j < theta_size; j++) {
            double const dOmega = std::fabs(std::cos(coord.theta_b[j + 1]) - std::cos(coord.theta_b[j])) * dphi[i];

            double log_I_lo = 0;
            double log_I_hi = 0;
            for (size_t k = 0, t_idx = 0, k_old = 0; k < r_size - 1 && t_idx < t_size; k++) {
                double const t_obs_lo = t_obs_grid[i][j][k];
                double const t_obs_hi = t_obs_grid[i][j][k + 1];

                if (k == 0) {
                    while (t_idx < t_size - 1 && t_obs[t_idx] < t_obs_lo) ++t_idx;
                }

                if (t_obs_lo <= t_obs[t_idx] && t_obs[t_idx] < t_obs_hi) {
                    if (k_old == k && k != 0) {
                        log_I_lo = log_I_hi;
                    } else {
                        log_I_lo = fastLog(getSpecificIntensity(i, j, k, nu_obs, photons...));
                    }
                    log_I_hi = fastLog(getSpecificIntensity(i, j, k + 1, nu_obs, photons...));
                    k_old = k + 1;
                }

#pragma omp simd
                for (; t_obs_lo <= t_obs[t_idx] && t_obs[t_idx] < t_obs_hi && t_idx < t_size; t_idx++) {
                    double const log_t = fastLog(t_obs[t_idx]);
                    double const D = dopplerInterp(log_t, i, j, k);
                    double const r = radiusInterp(log_t, i, j, k);
                    double I_nu = intensityInterp(log_t, i, j, k, log_I_lo, log_I_hi);

                    f_nu[t_idx] += D * D * D * I_nu * r * r * dOmega;
                }
            }
        }
    }

    // normalize the flux
    double coef = (1 + z) / (lumi_dist * lumi_dist);
    for (size_t i = 0; i < t_size; ++i) {
        f_nu[i] *= coef;
    }
}

template <typename... PhotonMesh>
Array Observer::specificFlux(Array const& t_obs, double nu_obs, PhotonMesh const&... photons) {
    Array F_nu = zeros(t_obs.size());
    calcSpecificFlux(F_nu.data(), t_obs, nu_obs, photons...);
    return F_nu;
}

template <typename... PhotonMesh>
MeshGrid Observer::specificFlux(Array const& t_obs, Array const& nu_obs, PhotonMesh const&... photons) {
    MeshGrid F_nu = createGrid(nu_obs.size(), t_obs.size(), 0);
    size_t t_num = t_obs.size();
    for (size_t l = 0; l < nu_obs.size(); ++l) {
        calcSpecificFlux(F_nu.data() + l * t_num, t_obs, nu_obs[l], photons...);
    }
    return F_nu;
}

template <typename... PhotonMesh>
Array Observer::flux(Array const& t_obs, Array const& band_pass_freq, PhotonMesh const&... photons) {
    Array nu_obs = boundaryToCenterLog(band_pass_freq);
    MeshGrid F_nu = specificFlux(t_obs, nu_obs, photons...);
    Array flux = zeros(t_obs.size());
    for (size_t i = 0; i < F_nu.size(); ++i) {
        double dnu = band_pass_freq[i + 1] - band_pass_freq[i];
        for (size_t j = 0; j < flux.size(); ++j) {
            flux[j] += dnu * F_nu[i][j];
        }
    }
    return flux;
}
#endif