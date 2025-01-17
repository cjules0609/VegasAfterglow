#ifndef _OBSERVER_
#define _OBSERVER_

#include <iostream>
#include <stdexcept>
#include <thread>
#include <vector>

#include "macros.h"
#include "mesh.h"

class Observer {
   public:
    Observer(Coord const& coord);

    Observer() = delete;

    void observe(MeshGrid const& Gamma, double theta_obs, double lumi_dist, double z);

    template <typename... PhotonMesh>
    MeshGrid3d specific_flux_grid(double nu_obs, PhotonMesh const&... rad_ptc);

    template <typename... PhotonMesh>
    Array specific_flux(Array const& t_bins, double nu_obs, PhotonMesh const&... photons) const;

    template <typename... PhotonMesh>
    MeshGrid specific_flux(Array const& t_bins, Array const& nu_obs, PhotonMesh const&... photons) const;

    template <typename... PhotonMesh>
    Array flux(Array const& t_bins, Array const& band_pass_freq, PhotonMesh const&... photons) const;

    template <typename... PhotonMesh>
    MeshGrid spectrum(Array const& t_bins, Array const& band_pass_freq, PhotonMesh const&... photons) const;

    MeshGrid3d doppler;
    MeshGrid3d t_obs;
    double theta_obs{0};
    double z{0};
    double lumi_dist{1};

   private:
    template <typename... PhotonMesh>
    void calc_specific_flux(Array& f_nu, Array const& t_bins, double nu_obs, PhotonMesh const&... photons) const;
    void calc_t_obs_grid(MeshGrid const& Gamma);
    // void calc_sorted_EAT_surface();
    // std::vector<std::pair<size_t, double>> eat_s;
    Coord const& coord;
    size_t effective_phi_size{1};
};

/*
template <typename... PhotonMesh>
void Observer::calc_specific_flux(Array& f_nu, Array const& t_bins, double nu_obs, const PhotonMesh&... photons) const {
    if (eat_s.empty()) {
        throw std::runtime_error("EAT surface is not defined. Please call observe() method first.");
    }
    // Variables for tracking and indexing
    size_t now_idx = 0;

    for (auto eat : eat_s) {
        size_t idx = eat.first;
        double t_obs = this->t_obs.data()[idx];
        if (t_obs < t_bins[0] || t_obs == 0 || std::isinf(t_obs)) {  // Skip unobserved and non-emission regions
            continue;
        }
        size_t i = idx / (coord.theta.size() * coord.r.size());
        size_t j = (idx % (coord.theta.size() * coord.r.size())) / coord.r.size();
        size_t k = idx % coord.r.size();
        // Compute the differential angle in phi
        double dphi = coord.phi_b[i + 1] - coord.phi_b[i];
        if (effective_phi_size == 1) {
            dphi = 2 * con::pi;
        }

        // Compute doppler shift
        double D = doppler.data()[idx];
        double nu_com = (1 + this->z) * nu_obs / D;

        // Compute differential energy per unit solid angle in comoving frame
        double dE_nu_dOmega_com = dphi / (2 * con::pi) * (photons[j][k].E_nu(nu_com) + ...) / (4 * con::pi);

        // Convert to observer frame
        double dE_nu_dOmega_obs = D * D * dE_nu_dOmega_com;

        // Accumulate flux in the appropriate time bin
        if (t_obs > t_bins[now_idx + 1] && now_idx + 1 < t_bins.size()) {
            ++now_idx;
        }

        if (now_idx < f_nu.size()) {
            f_nu[now_idx] += dE_nu_dOmega_obs;
        } else {
            break;
        }
    }

    // Normalize the flux for each bin
    for (size_t i = 0; i < f_nu.size(); ++i) {
        double dt_obs = t_bins[i + 1] - t_bins[i];
        f_nu[i] *= (1 + this->z) / (lumi_dist * lumi_dist) / dt_obs;
    }
}*/

int find_idx(Array const& T, double t);

template <typename... PhotonMesh>
void Observer::calc_specific_flux(Array& f_nu, Array const& t_bins, double nu_obs, const PhotonMesh&... photons) const {
    // Variables for tracking and indexing
    for (size_t i = 0; i < effective_phi_size; i++) {
        for (size_t j = 0; j < t_obs[0].size(); j++) {
            for (size_t k = 0; k < t_obs[0][0].size(); k++) {
                double t_obs = this->t_obs[i][j][k];
                int bin_id = find_idx(t_bins, t_obs);
                if (bin_id == -1 || !std::isfinite(t_obs)) {  // Skip unobserved and non-emission regions
                    continue;
                }

                double dphi = coord.phi_b[i + 1] - coord.phi_b[i];
                if (effective_phi_size == 1) {
                    dphi = 2 * con::pi;
                }
                double D = doppler[i][j][k];
                double nu_com = (1 + this->z) * nu_obs / D;

                // Compute differential energy per unit solid angle in comoving frame
                double dE_nu_dOmega_com = dphi / (2 * con::pi) * (photons[j][k].E_nu(nu_com) + ...) / (4 * con::pi);

                // Convert to observer frame
                double dE_nu_dOmega_obs = D * D * dE_nu_dOmega_com;

                f_nu[bin_id] += dE_nu_dOmega_obs;
            }
        }
    }

    for (size_t i = 0; i < f_nu.size(); ++i) {
        double dt_obs = t_bins[i + 1] - t_bins[i];
        f_nu[i] *= (1 + this->z) / (lumi_dist * lumi_dist) / dt_obs;
    }
}

template <typename... PhotonMesh>
MeshGrid3d Observer::specific_flux_grid(double nu_obs, PhotonMesh const&... photons) {
    MeshGrid3d F_nu_obs = create_3d_grid(effective_phi_size, coord.theta.size(), coord.r.size());

    for (size_t i = 0; i < F_nu_obs.size(); ++i) {
        double dphi = coord.phi_b[i + 1] - coord.phi_b[i];
        if (F_nu_obs.size() == 1) {
            dphi = 2 * con::pi;
        }
        for (size_t j = 0; j < F_nu_obs[0].size(); ++j) {
            for (size_t k = 0; k < F_nu_obs[0][0].size(); ++k) {
                double D = this->doppler[i][j][k];
                double nu_com = (1 + z) * nu_obs / D;
                F_nu_obs[i][j][k] = (1 + z) / (lumi_dist * lumi_dist) * D * D * D * dphi / (2 * con::pi) *
                                    (photons[j][k].L_nu(nu_com) + ...) / (4 * con::pi);
            }
        }
    }
    return F_nu_obs;
}

template <typename... PhotonMesh>
Array Observer::specific_flux(Array const& t_bins, double nu_obs, PhotonMesh const&... photons) const {
    Array F_nu = zeros(t_bins.size() - 1);
    calc_specific_flux(F_nu, t_bins, nu_obs, photons...);
    return F_nu;
}

template <typename... PhotonMesh>
MeshGrid Observer::specific_flux(Array const& t_bins, Array const& nu_obs, PhotonMesh const&... photons) const {
    MeshGrid F_nu = create_grid(nu_obs.size(), t_bins.size() - 1, 0);
    Array f_nu = zeros(t_bins.size() - 1);
    for (size_t l = 0; l < nu_obs.size(); ++l) {
        //  threads.emplace_back(&Observer::calc_specific_flux<PhotonMesh...>, this, std::ref(F_nu[l]),
        //  std::cref(t_bins),nu_obs[l], std::cref(photons)...);
        std::fill(f_nu.begin(), f_nu.end(), 0);
        calc_specific_flux(f_nu, t_bins, nu_obs[l], photons...);
        std::copy(f_nu.begin(), f_nu.end(), F_nu[l].begin());
    }
    // for (auto& thread : threads)
    //     thread.join();
    return F_nu;
}

template <typename... PhotonMesh>
Array Observer::flux(Array const& t_bins, Array const& band_pass_freq, PhotonMesh const&... photons) const {
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