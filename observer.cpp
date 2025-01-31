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
            double cos_v = std::sin(coord.theta[j]) * cos_phi * sin_obs + std::cos(coord.theta[j]) * cos_obs;

            for (size_t k = 0; k < coord.r.size(); ++k) {
                double gamma_ = Gamma[j][k];
                double beta = gammaTobeta(gamma_);
                doppler[i][j][k] = 1 / (gamma_ * (1 - beta * cos_v));
                if (gamma_ == 1) {
                    t_obs_grid[i][j][k] = std::numeric_limits<double>::infinity();
                } else {
                    t_obs_grid[i][j][k] = (t_eng[j][k] + (1 - cos_v) * coord.r[k] / con::c) * (1 + z);
                }
            }
        }
    }
}
