#ifndef _OBSERVER_
#define _OBSERVER_

#include <iostream>
#include <stdexcept>
#include <thread>
#include <vector>

#include "mesh.h"
#include "shock.h"

struct EATinfo {
    double t_obs;
    size_t i;  // index of phi
    size_t j;  // index of theta
    size_t k;  // index of r
    bool operator<(EATinfo const& other) const { return t_obs < other.t_obs; }
};

using EATsurface = std::vector<EATinfo>;

class Observer {
   public:
    void observe(Coord const& coord, Shock const& shock, double theta_obs, double z);

    template <typename... RadPhotonMesh>
    MeshGrid3d gen_F_nu_grid(double nu_obs, RadPhotonMesh const&... rad_ptc);

    template <typename... RadPhotonMesh>
    MeshGrid gen_F_nu(Array const& t_bins, Array const& nu_obs, RadPhotonMesh const&... photons) const;

    template <typename... RadPhotonMesh>
    Array gen_flux(Array const& t_bins, Array const& band_pass_freq, RadPhotonMesh const&... photons) const;

    template <typename... RadPhotonMesh>
    MeshGrid spectrum(double nu_min, double nu_max, double t, RadPhotonMesh const&... photons) const;

    double theta_obs{0};
    double z{0};
    double D_L{1};
    MeshGrid3d doppler;
    MeshGrid3d t_obs;

   private:
    template <typename... RadPhotonMesh>
    void gen_F_nu_(Array& F_nu, Array const& t_bins, double nu_obs, RadPhotonMesh const&... photons) const;
    void gen_phi_grid(Coord const& coord, double theta_obs);
    void calc_doppler_grid(Coord const& coord, MeshGrid const& Gamma);
    void calc_t_obs_grid(Coord const& coord, MeshGrid const& Gamma);
    void calc_sorted_EAT_surface(Coord const& coord, MeshGrid3d const& t_obs);
    void calc_D_L(double z);
    double first_non_zero_time() const;
    EATsurface eat_s;
    Array phi_b;
    Array phi;
};

template <typename... RadPhotonMesh>
void Observer::gen_F_nu_(std::vector<double>& f_nu, const std::vector<double>& t_bins, double nu_obs,
                         const RadPhotonMesh&... photons) const {
    if (eat_s.empty()) {
        throw std::runtime_error("EAT surface is not defined. Please call observe() method first.");
    }

    // Variables for tracking and indexing
    size_t next_idx = 0;

    for (const auto& rad_unit : eat_s) {
        double t_obs = rad_unit.t_obs;
        if (t_obs == 0) {  // Skip non-emission regions
            continue;
        }

        // Compute the differential angle in phi
        double dphi = this->phi_b[rad_unit.i + 1] - this->phi_b[rad_unit.i];

        // Compute doppler shift
        double D = doppler[rad_unit.i][rad_unit.j][rad_unit.k];
        double nu_com = (1 + this->z) * nu_obs / D;

        // Compute differential energy per unit solid angle in comoving frame
        double dE_nu_dOmega_com =
            dphi / (2 * con::pi) * (photons[rad_unit.j][rad_unit.k].E_nu(nu_com) + ...) / (4 * con::pi);

        // Convert to observer frame
        double dE_nu_dOmega_obs = D * D * dE_nu_dOmega_com;

        // Accumulate flux in the appropriate time bin
        if (t_obs > t_bins[next_idx + 1] && next_idx + 1 < t_bins.size()) {
            ++next_idx;
        }

        if (next_idx < f_nu.size()) {
            f_nu[next_idx] += dE_nu_dOmega_obs;
        } else {
            break;
        }
    }

    // Normalize the flux for each bin
    for (size_t i = 0; i < f_nu.size(); ++i) {
        double dt_obs = t_bins[i + 1] - t_bins[i];
        f_nu[i] *= (1 + this->z) / (D_L * D_L) / dt_obs;
    }
}

template <typename... RadPhotonMesh>
MeshGrid Observer::gen_F_nu(Array const& t_bins, Array const& nu_obs, RadPhotonMesh const&... photons) const {
    if (eat_s.empty()) {
        throw std::runtime_error("EAT surface is not defined. Please call observe() method first.");
    }

    MeshGrid F_nu = create_grid(nu_obs.size(), t_bins.size() - 1, 0);
    std::vector<std::thread> threads;
    for (size_t l = 0; l < nu_obs.size(); ++l) {
        threads.emplace_back(&Observer::gen_F_nu_<RadPhotonMesh...>, this, std::ref(F_nu[l]), std::cref(t_bins),
                             nu_obs[l], std::cref(photons)...);
    }
    for (auto& thread : threads) {
        thread.join();
    }
    return F_nu;
}

template <typename... RadPhotonMesh>
Array Observer::gen_flux(Array const& t_bins, Array const& band_pass_freq, RadPhotonMesh const&... photons) const {
    Array nu_obs = boundary2centerlog(band_pass_freq);
    MeshGrid F_nu = gen_F_nu(t_bins, nu_obs, photons...);
    Array flux = zeros(t_bins.size() - 1);
    for (size_t i = 0; i < F_nu.size(); ++i) {
        double dnu = band_pass_freq[i + 1] - band_pass_freq[i];
        for (size_t j = 0; j < flux.size(); ++j) {
            flux[j] += dnu * F_nu[i][j];
        }
    }
    return flux;
}

template <typename... RadPhotonMesh>
MeshGrid3d Observer::gen_F_nu_grid(double nu_obs, RadPhotonMesh const&... photons) {
    if (eat_s.empty()) {
        throw std::runtime_error("EAT surface is not defined. Please call observe() method first.");
    }
    MeshGrid3d F_nu_obs = create_3d_grid_like(doppler, 0);

    for (size_t i = 0; i < F_nu_obs.size(); ++i) {
        for (size_t j = 0; j < F_nu_obs[0].size(); ++j) {
            for (size_t k = 0; k < F_nu_obs[0][0].size(); ++k) {
                double D = this->doppler[i][j][k];
                double nu_com = (1 + this->z) * nu_obs / D;
                double dphi = this->phi_b[i + 1] - this->phi_b[i];

                F_nu_obs[i][j][k] = (1 + this->z) / (this->D_L * this->D_L) * D * D * D * dphi / (2 * con::pi) *
                                    (photons[j][k].L_nu(nu_com) + ...) / (4 * con::pi);
            }
        }
    }
    return F_nu_obs;
}

#endif