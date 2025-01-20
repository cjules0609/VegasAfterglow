#ifndef _OBSERVER_
#define _OBSERVER_

#include <functional>
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

    void observe(MeshGrid const& Gamma, double theta_obs, double lumi_dist, double z);

    template <typename... PhotonMesh>
    Array specific_flux(Array const& t_bins, double nu_obs, PhotonMesh const&... photons);

    template <typename... PhotonMesh>
    MeshGrid specific_flux(Array const& t_bins, Array const& nu_obs, PhotonMesh const&... photons);

    template <typename... PhotonMesh>
    Array flux(Array const& t_bins, Array const& band_pass_freq, PhotonMesh const&... photons);

    template <typename... PhotonMesh>
    MeshGrid spectrum(Array const& t_bins, Array const& band_pass_freq, PhotonMesh const&... photons);

    MeshGrid3d doppler;
    MeshGrid3d t_obs;
    MeshGrid3d flux_grid;

    double theta_obs{0};
    double z{0};
    double lumi_dist{1};

   private:
    template <typename Iter, typename... PhotonMesh>
    void calc_specific_flux(Iter f_nu, Array const& t_bins, double nu_obs, PhotonMesh const&... photons) const;
    void calc_t_obs_grid(MeshGrid const& Gamma);
    void optimize_idx_search(Array const& t_bins);
    int find_idx(Array const& T, double t) const;
    MeshGrid3d lg_t_obs;
    Coord const& coord;
    double t_max{0};
    double t_min{0};
    double t_space{1};
    MeshGrid3d* t_tab{&t_obs};
    bool optimized_search{false};
    size_t effective_phi_size{1};
};

template <typename Iter, typename... PhotonMesh>
void Observer::calc_specific_flux(Iter f_nu, Array const& t_bins, double nu_obs, const PhotonMesh&... photons) const {
    size_t const thread_num = th_pool.get_thread_count();
    MeshGrid f_grid = create_grid(thread_num, t_bins.size() - 1, 0);

    size_t const theta_size = t_obs[0].size();
    size_t const r_size = t_obs[0][0].size();

    th_pool.detach_blocks(0, effective_phi_size, [&](size_t start, size_t end) {
        auto const thread_id = BS::this_thread::get_index().value();
        for (size_t i = start; i < end; i++) {
            double dphi = 2 * con::pi;
            if (effective_phi_size != 1) {
                dphi = coord.phi_b[i + 1] - coord.phi_b[i];
            }
            for (size_t j = 0; j < theta_size; j++) {
                for (size_t k = 0; k < r_size; k++) {
                    int const bin_id = find_idx(t_bins, (*t_tab)[i][j][k]);
                    if (bin_id == -1) {  // Skip unobserved and non-emission regions
                        continue;
                    }

                    double const D = doppler[i][j][k];
                    double const nu_com = (1 + this->z) * nu_obs / D;

                    // Compute differential energy per unit solid angle in comoving frame
                    double const dE_nu_com = dphi * (photons[j][k].E_nu(nu_com) + ...);

                    f_grid[thread_id][bin_id] += D * D * dE_nu_com;
                }
            }
        }
    });
    th_pool.wait();

    for (size_t i = 0; i < thread_num; i++) {
        for (size_t j = 0; j < t_bins.size() - 1; j++) {
            f_nu[j] += f_grid[i][j];
        }
    }

    double coef = (1 + this->z) / (lumi_dist * lumi_dist * 4 * con::pi);

    for (size_t i = 0; i < t_bins.size() - 1; ++i) {
        double dt_obs = t_bins[i + 1] - t_bins[i];
        f_nu[i] *= coef / dt_obs;
    }
}

template <typename... PhotonMesh>
Array Observer::specific_flux(Array const& t_bins, double nu_obs, PhotonMesh const&... photons) {
    Array F_nu = zeros(t_bins.size() - 1);
    optimize_idx_search(t_bins);
    calc_specific_flux(F_nu, t_bins, nu_obs, photons...);
    return F_nu;
}

template <typename... PhotonMesh>
MeshGrid Observer::specific_flux(Array const& t_bins, Array const& nu_obs, PhotonMesh const&... photons) {
    MeshGrid F_nu = create_grid(nu_obs.size(), t_bins.size() - 1, 0);
    size_t bin_num = (t_bins.size() - 1);
    optimize_idx_search(t_bins);

    for (size_t l = 0; l < nu_obs.size(); ++l) {
        calc_specific_flux(F_nu.data() + l * bin_num, t_bins, nu_obs[l], photons...);
    }
    return F_nu;
}

template <typename... PhotonMesh>
Array Observer::flux(Array const& t_bins, Array const& band_pass_freq, PhotonMesh const&... photons) {
    Array nu_obs = boundary2centerlog(band_pass_freq);
    MeshGrid F_nu = specific_flux(t_bins, nu_obs, photons...);
    Array flux = zeros(t_bins.size() - 1);
    for (size_t i = 0; i < F_nu.size(); ++i) {
        double dnu = band_pass_freq[i + 1] - band_pass_freq[i];
        for (size_t j = 0; j < flux.size(); ++j) {
            flux[j] += dnu * F_nu[i][j];
        }
    }
    return flux;
}

#endif