#include "observer.h"

#include <boost/numeric/odeint.hpp>
#include <cmath>

#include "macros.h"
#include "physics.h"
#include "utilities.h"
void Observer::observe(Coord const& coord, Shock const& shock, double theta_obs, double lumi_dist, double z) {
    this->theta_obs = theta_obs;
    this->z = z;
    this->lumi_dist = lumi_dist;
    gen_phi_grid(coord, theta_obs);
    calc_doppler_grid(coord, shock.Gamma);
    calc_t_obs_grid(coord, shock.Gamma);
    calc_sorted_EAT_surface(coord, t_obs);
}

void Observer::gen_phi_grid(Coord const& coord, double theta_obs) {
    if (theta_obs == 0) {
        phi_b = {0, 2 * con::pi};
    } else {
        phi_b = coord.phi_b;
    }
    phi = boundary2center(phi_b);
}

double Observer::first_non_zero_time() const {
    for (auto const& eat : eat_s) {
        if (eat.t_obs != 0) {
            return eat.t_obs;
        }
    }
    return 0;
}

void Observer::calc_t_obs_grid(Coord const& coord, MeshGrid const& Gamma) {
    t_obs = create_3d_grid(this->phi.size(), coord.theta.size(), coord.r.size());
    using namespace boost::numeric::odeint;
    double atol = 0;
    double rtol = 1e-6;
    auto stepper = make_dense_output(atol, rtol, runge_kutta_dopri5<double>());

    double delta_r = (coord.r_b[1] - coord.r_b[0]) / 100;

    for (size_t i = 0; i < this->phi.size(); ++i) {
        double phi_ = this->phi[i];
        for (size_t j = 0; j < coord.theta.size(); ++j) {
            double theta_ = coord.theta[j];
            double cos_ = sin(theta_) * cos(phi_) * sin(theta_obs) + cos(theta_) * cos(theta_obs);

            auto eqn = [&](double const& t, double& dtdr, double r) {
                double Gamma_ = interp(r, coord.r, Gamma[j]);
                double beta = gamma_to_beta(Gamma_);
                dtdr = (1 / beta - cos_) / con::c;
            };

            double Gamma0 = Gamma[j][0];
            if (Gamma0 == 1) {
                std::fill(t_obs[i][j].begin(), t_obs[i][j].end(), 0);
                continue;
            }

            double beta0 = gamma_to_beta(Gamma0);
            t_obs[i][j][0] = coord.r[0] * (1 - beta0 * cos_) / con::c / beta0;

            stepper.initialize(t_obs[i][j][0], coord.r[0], delta_r);
            for (size_t k = 0; stepper.current_time() <= coord.r.back();) {
                stepper.do_step(eqn);
                for (; stepper.current_time() > coord.r[k + 1] && k + 1 < coord.r.size();) {
                    k++;
                    stepper.calc_state(coord.r[k], t_obs[i][j][k]);
                }
            }
        }
    }
}

void Observer::calc_doppler_grid(Coord const& coord, MeshGrid const& Gamma) {
    this->doppler = create_3d_grid(this->phi.size(), coord.theta.size(), coord.r.size());
    for (size_t i = 0; i < this->phi.size(); ++i) {
        double phi_ = this->phi[i];
        for (size_t j = 0; j < coord.theta.size(); ++j) {
            double theta_ = coord.theta[j];
            double cos_ = sin(theta_) * cos(phi_) * sin(theta_obs) + cos(theta_) * cos(theta_obs);
            for (size_t k = 0; k < Gamma[j].size(); ++k) {
                double gamma_ = Gamma[j][k];
                double beta = gamma_to_beta(gamma_);
                this->doppler[i][j][k] = 1 / (gamma_ * (1 - beta * cos_));
            }
        }
    }
}

void Observer::calc_sorted_EAT_surface(Coord const& coord, MeshGrid3d const& t_obs) {
    this->eat_s.resize(t_obs.size() * t_obs[0].size() * t_obs[0][0].size());
    size_t idx = 0;
    for (size_t i = 0; i < t_obs.size(); ++i) {
        for (size_t j = 0; j < t_obs[0].size(); ++j) {
            for (size_t k = 0; k < t_obs[0][0].size(); ++k) {
                this->eat_s[idx++] = EATinfo{t_obs[i][j][k], i, j, k};
            }
        }
    }
    std::sort(eat_s.begin(), eat_s.end());
}
