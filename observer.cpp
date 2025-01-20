#include "observer.h"

#include <boost/numeric/odeint.hpp>
#include <cmath>

#include "macros.h"
#include "physics.h"
#include "utilities.h"

Observer::Observer(Coord const& coord)
    : coord(coord),
      doppler(create_3d_grid(coord.phi.size(), coord.theta.size(), coord.r.size())),
      t_obs(create_3d_grid(coord.phi.size(), coord.theta.size(), coord.r.size())),
      flux_grid(create_3d_grid(coord.phi.size(), coord.theta.size(), coord.r.size())),
      lg_t_obs(create_3d_grid(coord.phi.size(), coord.theta.size(), coord.r.size())) {}

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
}

inline double fast_log2(double x) {
    uint64_t bx = std::bit_cast<uint64_t>(x);
    uint64_t ex = bx >> 52;
    int64_t t = static_cast<int64_t>(ex) - 1023;
    bx = 0x3FF0000000000000ULL | (bx & 0x000FFFFFFFFFFFFFULL);
    x = std::bit_cast<double>(bx);

    return -2.153626303138227 + (3.0478808242332507 + (-1.0518747250923068 + 0.15824921903511038 * x) * x) * x + t;
}

void Observer::calc_t_obs_grid(MeshGrid const& Gamma) {
    double cos_obs = cos(theta_obs);
    double sin_obs = sin(theta_obs);
    th_pool.detach_blocks(0, effective_phi_size, [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
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
                    this->lg_t_obs[i][j][k] = log2(t_obs[i][j][k]);
                    this->doppler[i][j][k] = 1 / (gamma_ * (1 - beta * cos_));
                }
            }
        }
    });
    th_pool.wait();
}

int Observer::find_idx(Array const& T, double t) const {
    if (t <= t_min || t >= t_max || !std::isfinite(t)) {
        return -1;  // t is out of range
    }

    if (optimized_search) {
        return static_cast<int>(std::floor((t - t_min) / t_space));
    } else {  // binary search
        int n = T.size();
        int low = 0, high = n - 1;
        while (low <= high) {
            int mid = low + (high - low) / 2;

            if (T[mid] < t && t < T[mid + 1]) {
                return mid;
            } else if (T[mid] < t) {
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
        return -1;
    }
}

void Observer::optimize_idx_search(Array const& t_bins) {
    if (is_log_scale(t_bins)) {  // optimization for log scale
        t_min = log2(t_bins[0]);
        t_max = log2(t_bins[t_bins.size() - 1]);
        t_space = log2(t_bins[1]) - t_min;
        t_tab = &lg_t_obs;
        optimized_search = true;
    } else if (is_linear_scale(t_bins)) {  // optimization for linear scale
        t_min = t_bins[0];
        t_max = t_bins[t_bins.size() - 1];
        t_space = t_bins[1] - t_min;
        t_tab = &t_obs;
        optimized_search = true;
    } else {
        t_min = t_bins[0];
        t_max = t_bins[t_bins.size() - 1];
        t_space = 0;
        t_tab = &t_obs;
        optimized_search = false;
    }
}
