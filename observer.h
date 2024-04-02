#ifndef _OBSERVER_
#define _OBSERVER_

#include <vector>

#include "mesh.h"

struct Indexes {
    size_t i, j, k;
    bool operator<(Indexes const& other) const { return i < other.i; }
};

using EATsurface = std::vector<std::pair<double, Indexes>>;

class Observer {
   public:
    void observe(Coord const& coord, MeshGrid const& Gamma, double theta_obs);

    template <typename... RadMesh>
    MeshGrid3d Observer::I_nu_history(double nu_obs, RadMesh const&... rad_com);

    template <typename... RadMesh>
    MeshGrid light_curve(size_t data_points, double nu_obs, RadMesh const&... rad_com) const;

    double theta_obs;
    MeshGrid3d Doppler;
    MeshGrid3d t_obs;

   private:
    void calc_doppler_grid(Coord const& coord, MeshGrid const& Gamma);
    void calc_t_obs_grid(Coord const& coord, MeshGrid const& Gamma);
    void get_sorted_EAT_surface(MeshGrid3d const& t_obs);
    void calc_emission_surface(Coord const& coord);
    EATsurface eat_s;
    MeshGrid3d emission_S;
};

template <typename... RadMesh>
MeshGrid Observer::light_curve(size_t data_points, double nu_obs, RadMesh const&... rad_com) const {
    if (eat_s.empty()) {
        throw std::runtime_error("EAT surface is not defined. Please call observe() method first.");
    }
    Array t_bin = logspace(eat_s.front().first, eat_s.back().first, data_points + 1);
    Array t_c = boundary2center(t_bin);
    MeshGrid L_nu = createGrid(2, data_points, 0);

    for (size_t m = 0, n = 0; m < eat_s.size(); ++m) {
        double t_ = eat_s[m].first;
        size_t i_ = eat_s[m].second.i;
        size_t j_ = eat_s[m].second.j;
        size_t k_ = eat_s[m].second.k;

        double doppler_ = Doppler[i_][j_][k_];

        double nu_prime = nu_obs / doppler_;
        double I_nu_tot = (rad_com[j_][k_].I_nu(nu_prime) + ...);
        double dL_nu_ = this->emission_S[i_][j_][k_] * doppler_ * doppler_ * doppler_ * I_nu_tot;

        if (t_ > t_bin[n + 1]) {
            n++;
        }
        if (n < L_nu[0].size()) {
            L_nu[1][n] += dL_nu_;
        } else {
            break;
        }
    }
    for (size_t i = 0; i < L_nu[0].size(); ++i) {
        L_nu[0][i] = t_c[i];
        L_nu[1][i] /= (t_bin[i + 1] - t_bin[i]);
    }
    return L_nu;
}

template <typename... RadMesh>
MeshGrid3d Observer::I_nu_history(double nu_obs, RadMesh const&... rad_com) {
    MeshGrid3d I_nu_obs = createGrid3d(coord.r, coord.theta, coord.phi);

    for (size_t i = 0; i < coord.phi.size(); ++i) {
        for (size_t j = 0; j < coord.theta.size(); ++j) {
            for (size_t k = 0; k < coord.r.size(); ++k) {
                double doppler = this->Doppler[i][j][k];
                double nu_prime = nu_obs / doppler;
                double I_nu_tot = (rad_com[j][k].I_nu(nu_prime) + ...);

                I_nu_obs[i][j][k] = doppler * doppler * doppler * I_nu_tot;
            }
        }
    }
    return I_nu_obs;
}

#endif