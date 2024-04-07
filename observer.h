#ifndef _OBSERVER_
#define _OBSERVER_

#include <iostream>
#include <stdexcept>
#include <vector>

#include "mesh.h"
#include "shock.h"
struct Indexes {
    size_t i, j, k;
    bool operator<(Indexes const& other) const { return i < other.i; }
};

using EATsurface = std::vector<std::pair<double, Indexes>>;

class Observer {
   public:
    void observe(Coord const& coord, Shock const& shock, double theta_obs);

    template <typename... RadParticleMesh>
    MeshGrid3d gen_I_nu_grid(double nu_obs, RadParticleMesh const&... rad_ptc);

    template <typename... RadParticleMesh>
    MeshGrid gen_light_curve(size_t time_resolution, Array const& nu_obs, RadParticleMesh const&... rad_ptc) const;

    template <typename... RadParticleMesh>
    MeshGrid spectrum(double nu_min, double nu_max, double t, RadParticleMesh const&... rad_ptc) const;

    double theta_obs{0};
    double D_L{1};
    MeshGrid3d doppler;
    MeshGrid3d t_obs;

   private:
    void calc_doppler_grid(Coord const& coord, MeshGrid const& Gamma);
    void calc_t_obs_grid(Coord const& coord, MeshGrid const& Gamma);
    void calc_sorted_EAT_surface(MeshGrid3d const& t_obs);
    void calc_emission_volume(Coord const& coord, MeshGrid const& Gamma, MeshGrid const& D_com);
    EATsurface eat_s;
    MeshGrid3d emission_V;
};

template <typename... RadParticleMesh>
MeshGrid Observer::gen_light_curve(size_t time_resolution, Array const& nu_obs,
                                   RadParticleMesh const&... rad_ptc) const {
    if (eat_s.empty()) {
        throw std::runtime_error("EAT surface is not defined. Please call observe() method first.");
    }
    Array t_bin = logspace(eat_s.front().first, eat_s.back().first, time_resolution + 1);
    Array t_c = boundary2center(t_bin);
    MeshGrid F_nu = create_grid(nu_obs.size() + 1, time_resolution, 0);

    for (size_t l = 0; l < nu_obs.size(); ++l) {
        for (size_t m = 0, n = 0; m < eat_s.size(); ++m) {
            double t_ = eat_s[m].first;
            size_t i_ = eat_s[m].second.i;
            size_t j_ = eat_s[m].second.j;
            size_t k_ = eat_s[m].second.k;

            double doppler_ = doppler[i_][j_][k_];

            double nu_prime = nu_obs[l] / doppler_;

            // emissivity in the observer frame
            double j_nu_tot = doppler_ * doppler_ * (rad_ptc[j_][k_].j_nu(nu_prime) + ...);

            double beaming = doppler_ * doppler_;

            double F_nu_ = this->emission_V[i_][j_][k_] * j_nu_tot * beaming;

            if (t_ > t_bin[n + 1]) {
                n++;
            }
            if (n < F_nu[0].size()) {
                F_nu[l + 1][n] += F_nu_;
            } else {
                break;
            }
        }
        for (size_t i = 0; i < F_nu[0].size(); ++i) {
            if (l == 0) {
                F_nu[0][i] = t_c[i];
            }
            F_nu[l + 1][i] /= ((t_bin[i + 1] - t_bin[i]) * (D_L * D_L));
        }
    }
    return F_nu;
}

/*template <typename... RadParticleMesh>
MeshGrid Observer::spectrum(double nu_min, double nu_max, double t, RadParticleMesh const&... rad_com) const {
    if (eat_s.empty()) {
        throw std::runtime_error("EAT surface is not defined. Please call observe() method first.");
    }
    size_t data_points = 100;
    Array nu_bin = logspace(nu_min, nu_max, data_points + 1);
    Array nu_c = boundary2center(nu_bin);
    MeshGrid L_nu = create_grid(2, data_points, 0);

    double t_min = eat_s.front().first;
    double t_max = eat_s.back().first;
}*/

template <typename... RadParticleMesh>
MeshGrid3d Observer::gen_I_nu_grid(double nu_obs, RadParticleMesh const&... rad_ptc) {
    if (eat_s.empty()) {
        throw std::runtime_error("EAT surface is not defined. Please call observe() method first.");
    }
    MeshGrid3d I_nu_obs = create_3d_grid_like(doppler, 0);

    for (size_t i = 0; i < I_nu_obs.size(); ++i) {
        for (size_t j = 0; j < I_nu_obs[0].size(); ++j) {
            for (size_t k = 0; k < I_nu_obs[0][0].size(); ++k) {
                double doppler = this->doppler[i][j][k];
                double nu_prime = nu_obs / doppler;
                double j_nu_tot = (rad_ptc[j][k].j_nu(nu_prime) + ...);

                I_nu_obs[i][j][k] = doppler * doppler * doppler * j_nu_tot;
            }
        }
    }
    return I_nu_obs;
}

#endif