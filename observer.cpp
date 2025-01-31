#include "observer.h"

#include <boost/numeric/odeint.hpp>
#include <cmath>

#include "macros.h"
#include "physics.h"
#include "utilities.h"

Observer::Observer(Coord const& coord)
    : coord(coord),
      doppler(create3DGrid(coord.phi.size(), coord.theta.size(), coord.r.size())),
      t_obs_grid(create3DGrid(coord.phi.size(), coord.theta.size(), coord.r.size())),
      log_t_obs(create3DGrid(coord.phi.size(), coord.theta.size(), coord.r.size())),
      log_doppler(create3DGrid(coord.phi.size(), coord.theta.size(), coord.r.size())),
      log_r(zeros(coord.r.size())),
      dphi(zeros(coord.phi.size())) {
    for (size_t i = 0; i < coord.r.size(); ++i) {
        log_r[i] = std::log(coord.r[i]);
    }
}

void Observer::calcObsTimeGrid(MeshGrid const& Gamma, MeshGrid const& t_eng) {
    double cos_obs = std::cos(theta_obs);
    double sin_obs = std::sin(theta_obs);
    for (size_t i = 0; i < effective_phi_size; ++i) {
        double cos_phi = std::cos(coord.phi[i]);
        for (size_t j = 0; j < coord.theta.size(); ++j) {
            double cos_theta = std::cos(coord.theta[j]);
            double sin_theta = std::sqrt(1 - cos_theta * cos_theta);
            double cos_v = sin_theta * cos_phi * sin_obs + cos_theta * cos_obs;

            for (size_t k = 0; k < coord.r.size(); ++k) {
                double gamma_ = Gamma[j][k];
                double beta = gammaTobeta(gamma_);
                doppler[i][j][k] = 1 / (gamma_ * (1 - beta * cos_v));
                if (gamma_ == 1) {
                    t_obs_grid[i][j][k] = std::numeric_limits<double>::infinity();
                } else {
                    t_obs_grid[i][j][k] = (t_eng[j][k] + (1 - cos_v) * coord.r[k] / con::c) / (1 + z);
                }
                log_t_obs[i][j][k] = fastLog(t_obs_grid[i][j][k]);
                log_doppler[i][j][k] = fastLog(doppler[i][j][k]);
            }
        }
    }
}

double Observer::dopplerInterp(double log_t, size_t i, size_t j, size_t k) const {
    double d_lo = log_doppler[i][j][k];
    double d_hi = log_doppler[i][j][k + 1];
    double t_lo = log_t_obs[i][j][k];
    double t_hi = log_t_obs[i][j][k + 1];
    return fastExp(d_lo + (d_hi - d_lo) * (log_t - t_lo) / (t_hi - t_lo));
}

double Observer::radiusInterp(double log_t, size_t i, size_t j, size_t k) const {
    double r_lo = log_r[k];
    double r_hi = log_r[k + 1];
    double t_lo = log_t_obs[i][j][k];
    double t_hi = log_t_obs[i][j][k + 1];
    return fastExp(r_lo + (r_hi - r_lo) * (log_t - t_lo) / (t_hi - t_lo));
}

// call photon.I_nu(nu) is expensive, and not all grids are needed, so we don't use pre-calculated lookup tables
double Observer::intensityInterp(double log_t, size_t i, size_t j, size_t k, double log_I_lo, double log_I_hi) const {
    double t_lo = log_t_obs[i][j][k];
    double t_hi = log_t_obs[i][j][k + 1];

    return fastExp(log_I_lo + (log_I_hi - log_I_lo) * (log_t - t_lo) / (t_hi - t_lo));
}
