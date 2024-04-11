#include "observer.h"

#include <boost/numeric/odeint.hpp>
#include <cmath>

#include "macros.h"
#include "utilities.h"

void Observer::observe(Coord& coord, Shock const& shock, double theta_obs, double z) {
    if (theta_obs == 0 && coord.phi.size() >= 1) {
        std::cout << "on-axis observer, reducing phi grid size to 1\n";
        coord.phi_b = Array{0, 2 * con::pi};
        coord.phi = Array{con::pi};
    }

    this->theta_obs = theta_obs;
    calc_doppler_grid(coord, shock.Gamma);
    calc_t_obs_grid(coord, shock.Gamma);
    calc_sorted_EAT_surface(t_obs);
    calc_emission_volume(coord, shock.Gamma, shock.D_com);
    calc_D_L(z);
}

void Observer::calc_D_L(double z) {
    using namespace boost::numeric::odeint;
    double atol = 0;
    double rtol = 1e-6;
    auto stepper = bulirsch_stoer_dense_out<double>{atol, rtol};
    // auto stepper = make_dense_output(1.0e-6, 1.0e-6, runge_kutta_cash_karp54<std::vector<double>>());

    auto eqn = [&](double const& y, double& dydz, double z0) {
        dydz = 1 / sqrt(con::Omega_m * (1 + z0) * (1 + z0) * (1 + z0) + con::Omega_L);
    };
    stepper.initialize(0, 0, 1e-6);

    for (; stepper.current_time() <= z;) {
        stepper.do_step(eqn);
    }
    stepper.calc_state(z, this->D_L);
    this->D_L = (1 + z) * this->D_L * con::c / con::H0;
}

void Observer::calc_t_obs_grid(Coord const& coord, MeshGrid const& Gamma) {
    t_obs = create_3d_grid(coord.phi.size(), coord.theta.size(), coord.r.size());
    using namespace boost::numeric::odeint;
    double atol = 0;
    double rtol = 1e-6;
    auto stepper = bulirsch_stoer_dense_out<double>{atol, rtol};

    double dr0 = (coord.r_b[1] - coord.r_b[0]) / 1000;
    for (size_t i = 0; i < coord.phi.size(); ++i) {
        double phi_ = coord.phi[i];
        for (size_t j = 0; j < coord.theta.size(); ++j) {
            double theta_ = coord.theta[j];
            double cos_ = sin(theta_) * cos(phi_) * sin(theta_obs) + cos(theta_) * cos(theta_obs);

            auto eqn = [&](double const& t, double& dtdr, double r) {
                double Gamma_ = interp(r, coord.r, Gamma[j]);
                double beta = sqrt(1 - 1 / Gamma_ / Gamma_);
                dtdr = (1 - beta * cos_) / (beta * con::c);
            };

            double Gamma0 = Gamma[j][0];
            double beta0 = sqrt(1 - 1 / Gamma0 / Gamma0);
            t_obs[i][j][0] = coord.r[0] * (1 - beta0 * cos_) / con::c / beta0;

            stepper.initialize(t_obs[i][j][0], coord.r[0], dr0);
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
    this->doppler = create_3d_grid(coord.phi.size(), coord.theta.size(), coord.r.size());
    for (size_t i = 0; i < coord.phi.size(); ++i) {
        double phi_ = coord.phi[i];
        for (size_t j = 0; j < coord.theta.size(); ++j) {
            double theta_ = coord.theta[j];
            double cos_ = sin(theta_) * cos(phi_) * sin(theta_obs) + cos(theta_) * cos(theta_obs);
            for (size_t k = 0; k < Gamma[j].size(); ++k) {
                double gamma_ = Gamma[j][k];
                double beta = sqrt(1 - 1 / gamma_ / gamma_);
                this->doppler[i][j][k] = 1 / (gamma_ * (1 - beta * cos_));
            }
        }
    }
}

void Observer::calc_emission_volume(Coord const& coord, MeshGrid const& Gamma, MeshGrid const& D_com) {
    this->emission_V = create_3d_grid(coord.phi.size(), coord.theta.size(), coord.r.size());
    for (size_t i = 0; i < coord.phi.size(); ++i) {
        for (size_t j = 0; j < coord.theta.size(); ++j) {
            for (size_t k = 0; k < coord.r.size(); ++k) {
                double dcos = std::cos(coord.theta_b[j + 1]) - std::cos(coord.theta_b[j]);
                double dphi = coord.phi_b[i + 1] - coord.phi_b[i];
                double dOmega = std::fabs(dphi * dcos);
                double r_ = coord.r[k];
                double Gamma_ = Gamma[j][k];
                double shock_width = D_com[j][k] * Gamma_;
                this->emission_V[i][j][k] = r_ * r_ * dOmega * shock_width;
            }
        }
    }
}

void Observer::calc_sorted_EAT_surface(MeshGrid3d const& t_obs) {
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
