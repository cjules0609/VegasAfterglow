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
    MeshGrid gen_light_curve(size_t time_resolution, Array const& nu_obs, RadPhotonMesh const&... rad_ptc) const;

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
    EATsurface eat_s;
    Array phi_b;
    Array phi;
};

template <typename... RadPhotonMesh>
MeshGrid Observer::gen_light_curve(size_t time_resolution, Array const& nu_obs, RadPhotonMesh const&... rad_ptc) const {
    if (eat_s.empty()) {
        throw std::runtime_error("EAT surface is not defined. Please call observe() method first.");
    }
    Array t_bin = logspace(eat_s.front().t_obs, eat_s.back().t_obs, time_resolution + 1);
    Array t_c = boundary2center(t_bin);
    MeshGrid F_nu = create_grid(nu_obs.size() + 1, time_resolution, 0);

    for (size_t l = 0; l < nu_obs.size(); ++l) {
        for (size_t m = 0, n = 0; m < eat_s.size(); ++m) {
            double t_ = eat_s[m].t_obs;
            size_t i_ = eat_s[m].i;
            size_t j_ = eat_s[m].j;
            size_t k_ = eat_s[m].k;

            double dphi = this->phi_b[i_ + 1] - this->phi_b[i_];

            double doppler_ = doppler[i_][j_][k_];

            double nu_prime = nu_obs[l] / doppler_;

            double dL_nu_dOmega_com = dphi / (2 * con::pi) * (rad_ptc[j_][k_].L_nu(nu_prime) + ...) / (4 * con::pi);

            double dL_nu_dOmega_obs = doppler_ * doppler_ * doppler_ * dL_nu_dOmega_com;

            if (t_ > t_bin[n + 1]) {
                n++;
            }
            if (n < F_nu[0].size()) {
                F_nu[l + 1][n] += dL_nu_dOmega_obs;
            } else {
                break;
            }
        }
        for (size_t i = 0; i < F_nu[0].size(); ++i) {
            if (l == 0) {
                F_nu[0][i] = t_c[i];
            }
            F_nu[l + 1][i] *= (1 + this->z) / (D_L * D_L);
        }
    }
    return F_nu;
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