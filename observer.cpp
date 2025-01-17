#include "observer.h"

#include <boost/numeric/odeint.hpp>
#include <cmath>

#include "macros.h"
#include "physics.h"
#include "utilities.h"
Observer::Observer(Coord const& coord)
    : coord(coord),
      doppler(create_3d_grid(coord.phi.size(), coord.theta.size(), coord.r.size())),
      t_obs(create_3d_grid(coord.phi.size(), coord.theta.size(), coord.r.size())) {}

void Observer::observe(MeshGrid const& Gamma, double theta_obs, double lumi_dist, double z) {
    this->theta_obs = theta_obs;
    this->z = z;
    this->lumi_dist = lumi_dist;
    if (theta_obs == 0) {
        effective_phi_size = 1;
    } else {
        effective_phi_size = coord.phi.size();
    }
    calc_t_obs_grid(Gamma);
    // calc_sorted_EAT_surface();
}

void Observer::calc_t_obs_grid(MeshGrid const& Gamma) {
    double cos_obs = cos(theta_obs);
    double sin_obs = sin(theta_obs);
    for (size_t i = 0; i < effective_phi_size; ++i) {
        double phi_ = coord.phi[i];
        double cosphi = cos(phi_);
        for (size_t j = 0; j < coord.theta.size(); ++j) {
            double theta_ = coord.theta[j];
            double cos_ = sin(theta_) * cosphi * sin_obs + cos(theta_) * cos_obs;
            for (size_t k = 0; k < coord.r.size(); ++k) {
                double gamma_ = Gamma[j][k];
                double beta = gamma_to_beta(gamma_);
                if (k == 0) {
                    t_obs[i][j][k] = coord.r_b[0] * (1 / beta - cos_) / con::c;
                } else {
                    double dr = coord.r_b[k] - coord.r_b[k - 1];
                    this->t_obs[i][j][k] = dr * (1 / beta - cos_) / con::c + t_obs[i][j][k - 1];
                }
                this->doppler[i][j][k] = 1 / (gamma_ * (1 - beta * cos_));
            }
        }
    }
}

int find_idx(Array const& T, double t) {
    int n = T.size();
    if (t <= T[0] || t >= T[n - 1]) {
        return -1;  // t is out of range
    }

    int low = 0, high = n - 1;
    while (low <= high) {
        int mid = low + (high - low) / 2;

        if (T[mid] < t && t < T[mid + 1]) {
            return mid;  // Found the index
        } else if (T[mid] < t) {
            low = mid + 1;  // Search in the right half
        } else {
            high = mid - 1;  // Search in the left half
        }
    }

    return -1;  // No valid index found
}
/*
void Observer::calc_sorted_EAT_surface() {
    size_t theta_size = t_obs[0].size();
    size_t r_size = t_obs[0][0].size();
    this->eat_s.resize(effective_phi_size * theta_size * r_size);
    for (size_t i = 0; i < eat_s.size(); ++i) {
        eat_s[i] = std::make_pair(i, t_obs.data()[i]);
    }
    std::sort(eat_s.begin(), eat_s.end(), [&](auto& a, auto& b) { return a.second < b.second; });
}*/
