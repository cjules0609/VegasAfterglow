#ifndef _OBSERVER_
#define _OBSERVER_

#include <iostream>
#include <stdexcept>
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
    MeshGrid gen_F_nu(Array const& t_bins, Array const& nu_obs, RadPhotonMesh const&... rad_ptc) const;

    template <typename... RadPhotonMesh>
    Array gen_Flux(Array const& t_bins, double band_filter_low, double band_filter_hi, size_t freq_resolution,
                   RadPhotonMesh const&... rad_ptc) const;

    template <typename... RadPhotonMesh>
    MeshGrid spectrum(double nu_min, double nu_max, double t, RadPhotonMesh const&... rad_ptc) const;

    double theta_obs{0};
    double z{0};
    double D_L{1};
    MeshGrid3d doppler;
    MeshGrid3d t_obs;

   private:
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
MeshGrid Observer::gen_F_nu(Array const& t_bins, Array const& nu_obs, RadPhotonMesh const&... rad_ptc) const {
    if (eat_s.empty()) {
        throw std::runtime_error("EAT surface is not defined. Please call observe() method first.");
    }

    MeshGrid F_nu = create_grid(nu_obs.size(), t_bins.size() - 1, 0);

    for (size_t l = 0; l < nu_obs.size(); ++l) {
        for (size_t m = 0, n = 0; m < eat_s.size(); ++m) {
            double t_ = eat_s[m].t_obs;
            size_t i_ = eat_s[m].i;
            size_t j_ = eat_s[m].j;
            size_t k_ = eat_s[m].k;

            if (t_ == 0) {  // non-emission region
                continue;
            }

            double dphi = this->phi_b[i_ + 1] - this->phi_b[i_];

            double doppler_ = doppler[i_][j_][k_];

            double nu_prime = (1 + this->z) * nu_obs[l] / doppler_;

            double dE_nu_dOmega_com = dphi / (2 * con::pi) * (rad_ptc[j_][k_].E_nu(nu_prime) + ...) / (4 * con::pi);

            double dE_nu_dOmega_obs = doppler_ * doppler_ * dE_nu_dOmega_com;

            if (t_ > t_bins[n + 1]) {
                n++;
            }
            if (n < F_nu[0].size()) {
                F_nu[l][n] += dE_nu_dOmega_obs;
                /* if (t_ / con::day > 50 && t_ / con::day < 500) {
                        double dt = (rad_ptc[j_][k_].dt_com + ...);
                        std::cout << t_ / con::sec << ' ' << dt / doppler_ << ' ' << t_bins[n + 1] - t_bins[n] << '\n';
                    }
                    */
            } else {
                break;
            }
        }
        for (size_t i = 0; i < F_nu[0].size(); ++i) {
            double dt = t_bins[i + 1] - t_bins[i];
            F_nu[l][i] *= (1 + this->z) / (D_L * D_L) / dt;
        }
    }
    // exit(0);
    return F_nu;
}

template <typename... RadPhotonMesh>
Array Observer::gen_Flux(Array const& t_bins, double band_filter_low, double band_filter_hi, size_t freq_resol,
                         RadPhotonMesh const&... rad_ptc) const {
    Array nu_obs_b = logspace(band_filter_low, band_filter_hi, freq_resol + 1);
    Array nu_obs = boundary2centerlog(nu_obs_b);
    MeshGrid F_nu = gen_F_nu(t_bins, nu_obs, rad_ptc...);
    Array flux = zeros(t_bins.size() - 1);
    for (size_t i = 0; i < F_nu.size(); ++i) {
        double dnu = nu_obs_b[i + 1] - nu_obs_b[i];
        for (size_t j = 0; j < flux.size(); ++j) {
            flux[j] += dnu * F_nu[i][j];
        }
    }
    return flux;
}

template <typename... RadPhotonMesh>
MeshGrid3d Observer::gen_F_nu_grid(double nu_obs, RadPhotonMesh const&... rad_ptc) {
    if (eat_s.empty()) {
        throw std::runtime_error("EAT surface is not defined. Please call observe() method first.");
    }
    MeshGrid3d F_nu_obs = create_3d_grid_like(doppler, 0);

    for (size_t i = 0; i < F_nu_obs.size(); ++i) {
        for (size_t j = 0; j < F_nu_obs[0].size(); ++j) {
            for (size_t k = 0; k < F_nu_obs[0][0].size(); ++k) {
                double doppler_ = this->doppler[i][j][k];
                double nu_prime = nu_obs / doppler_;
                double dphi = this->phi_b[i + 1] - this->phi_b[i];

                F_nu_obs[i][j][k] = (1 + this->z) / (this->D_L * this->D_L) * doppler_ * doppler_ * doppler_ * dphi /
                                    (2 * con::pi) * (rad_ptc[j][k].L_nu(nu_prime) + ...) / (4 * con::pi);
            }
        }
    }
    return F_nu_obs;
}

#endif