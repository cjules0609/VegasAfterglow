#include "observer.h"

#include <boost/numeric/odeint.hpp>
#include <cmath>

#include "macros.h"
#include "physics.h"
#include "utilities.h"

void LogScaleInterp::reset() { idx_hi = 0; }

bool LogScaleInterp::isfinite() const { return std::isfinite(log_t_lo) && std::isfinite(log_t_hi); }

double LogScaleInterp::interpRadius(double log_t) const {
    return fastExp(log_r_lo + (log_r_hi - log_r_lo) * (log_t - log_t_lo) / (log_t_hi - log_t_lo));
}

double LogScaleInterp::interpIntensity(double log_t) const {
    return fastExp(log_I_lo + (log_I_hi - log_I_lo) * (log_t - log_t_lo) / (log_t_hi - log_t_lo));
}

double LogScaleInterp::interpDoppler(double log_t) const {
    return fastExp(log_d_lo + (log_d_hi - log_d_lo) * (log_t - log_t_lo) / (log_t_hi - log_t_lo));
}

void Observer::calcSolidAngle() {
    for (size_t i = 0; i < eff_phi_size; ++i) {
        double dphi = coord.phi_b[i + 1] - coord.phi_b[i];
        if (eff_phi_size == 1) {
            dphi = 2 * con::pi;
        }
        for (size_t j = 0; j < coord.theta.size(); ++j) {
            dOmega[i][j] = std::fabs(std::cos(coord.theta_b[j + 1]) - std::cos(coord.theta_b[j])) * dphi;
        }
    }
}

Observer::Observer(Coord const& coord)
    : coord(coord),
      doppler(create3DGrid(coord.phi.size(), coord.theta.size(), coord.r.size())),
      t_obs_grid(create3DGrid(coord.phi.size(), coord.theta.size(), coord.r.size())),
      dOmega(createGrid(coord.phi.size(), coord.theta.size())),
      log_r(zeros(coord.r.size())) {
    for (size_t i = 0; i < coord.r.size(); ++i) {
        log_r[i] = std::log(coord.r[i]);
    }

    calcSolidAngle();
}

void Observer::calcObsTimeGrid(MeshGrid3d const& Gamma, MeshGrid3d const& t_eng) {
    double cos_obs = std::cos(theta_obs);
    double sin_obs = std::sin(theta_obs);
    for (size_t i = 0; i < eff_phi_size; ++i) {
        double cos_phi = std::cos(coord.phi[i]);
        for (size_t j = 0; j < coord.theta.size(); ++j) {
            double cos_v = std::sin(coord.theta[j]) * cos_phi * sin_obs + std::cos(coord.theta[j]) * cos_obs;

            for (size_t k = 0; k < coord.r.size(); ++k) {
                double gamma_ = Gamma[i * interp.jet_3d][j][k];
                double t_eng_ = t_eng[i * interp.jet_3d][j][k];
                double beta = gammaTobeta(gamma_);
                doppler[i][j][k] = 1 / (gamma_ * (1 - beta * cos_v));
                if (gamma_ == 1) {
                    t_obs_grid[i][j][k] = std::numeric_limits<double>::infinity();
                } else {
                    t_obs_grid[i][j][k] = (t_eng_ + (1 - cos_v) * coord.r[k] / con::c) * (1 + z);
                }
            }
        }
    }
}
