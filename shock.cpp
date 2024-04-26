
#include "shock.h"

#include <boost/numeric/odeint.hpp>

#include "macros.h"

Shock::Shock(Coord const& coord, double eps_e, double eps_B, double xi, double zeta)
    : t_com(create_grid(coord.theta.size(), coord.r.size(), 0)),
      Gamma(create_grid(coord.theta.size(), coord.r.size(), 1)),
      B(create_grid(coord.theta.size(), coord.r.size(), 0)),
      width(create_grid(coord.theta.size(), coord.r.size(), 0)),
      n_p(create_grid(coord.theta.size(), coord.r.size(), 0)),
      eps_e{eps_e},
      eps_B{eps_B},
      xi{xi},
      zeta{zeta} {}

void co_moving_B(Shock& shock, Coord const& coord) {
    for (size_t j = 0; j < coord.theta.size(); ++j) {
        for (size_t k = 0; k < coord.r.size(); ++k) {
            double eps_B = shock.eps_B;
            double Gamma = shock.Gamma[j][k];
            double n_2 = shock.n_p[j][k];
            shock.B[j][k] = sqrt(8 * con::pi * eps_B * n_2 * con::mp * (Gamma - 1)) * con::c;
        }
    }
}

// lab frame FS width
void co_moving_FS_shock_width(Shock& shock, Coord const& coord) {
    for (size_t j = 0; j < coord.theta.size(); ++j) {
        for (size_t k = 0; k < coord.r.size(); ++k) {
            double Gamma = shock.Gamma[j][k];
            shock.width[j][k] = coord.r[k] / shock.Gamma[j][k];
        }
    }
}

void FS_n2(Shock& shock, Coord const& coord, Medium const& medium) {
    for (size_t j = 0; j < coord.theta.size(); ++j) {
        for (size_t k = 0; k < coord.r.size(); ++k) {
            shock.n_p[j][k] = 4 * shock.Gamma[j][k] * medium.rho(coord.r[k]) / con::mp;
        }
    }
}

BlastWaveEqn::BlastWaveEqn(Medium const& medium, Jet const& blast, double theta_lo, double theta_hi)
    : medium(medium),
      jet(blast),
      theta_lo(theta_lo),
      theta_hi(theta_hi),
      theta(0.5 * (theta_lo + theta_hi)),
      dOmega(4 * con::pi * std::fabs(std::cos(theta_hi) - std::cos(theta_lo))){};  // bipolar outflow

void BlastWaveEqn::operator()(Array const& y, Array& dydr, double r) {
    double Gamma = y[0];
    double u = y[1];
    double t_eng = y[2];  // engine time
    // double t_com = y[3];  // comoving time
    // double D_RS = y[4];   // comoving frame reverse shock width
    double D_FS = y[5];  // comoving frame forward shock width

    dydr[0] = dGammadr(r, Gamma, u, t_eng);
    dydr[1] = dUdr(r, Gamma, u, t_eng);

    dydr[2] = dtdr_eng(Gamma);
    dydr[3] = dtdr_com(Gamma);
    dydr[4] = dDdr_RS(r, Gamma, D_FS, t_eng);
    dydr[5] = dDdr_FS(Gamma);
};

double BlastWaveEqn::dGammadr(double r, double Gamma, double u, double t_eng) {
    double gamma_eos = (4 * Gamma + 1) / (3 * Gamma);
    double dm = medium.mass(r) * dOmega / (4 * con::pi);
    double dM0 = jet.dEdOmega(theta, t_eng) * dOmega / (jet.Gamma0(theta) * con::c2);
    double Gamma2 = Gamma * Gamma;
    double a1 = dOmega * r * r * medium.rho(r) / Gamma * (Gamma2 - 1) * (gamma_eos * Gamma - gamma_eos + 1);
    double a2 = -(gamma_eos - 1) / Gamma * (gamma_eos * Gamma2 - gamma_eos + 1) * 3 * u / r;
    double b1 = (dM0 + dm) * con::c2;
    double b2 = (gamma_eos * gamma_eos * (Gamma2 - 1) + 3 * gamma_eos - 2) * u / Gamma2;
    return -(a1 + a2) / (b1 + b2);
};

double BlastWaveEqn::dUdr(double r, double Gamma, double u, double t_eng) {
    double gamma_eos = (4 * Gamma + 1) / (3 * Gamma);
    double E = dOmega * r * r * medium.rho(r) * con::c2;
    return (1 - medium.eps_e) * (Gamma - 1) * E - (gamma_eos - 1) * (3 / r - dGammadr(r, Gamma, u, t_eng) / Gamma) * u;
};

double BlastWaveEqn::dtdr_eng(double Gamma) {
    double Gb = std::sqrt(Gamma * Gamma - 1);
    return (Gamma - Gb) / (Gb * con::c);
};

double BlastWaveEqn::dtdr_com(double Gamma) { return 1 / (std::sqrt(Gamma * Gamma - 1) * con::c); };  // co-moving time

double BlastWaveEqn::dDdr_FS(double Gamma) {
    double gamma_eos = (4 * Gamma + 1) / (3 * Gamma);
    double cs = con::c * std::sqrt(gamma_eos - 1);
    return 0;  // cs * dtdr_com(Gamma);
};            

double BlastWaveEqn::dDdr_RS(double r, double Gamma, double D_com, double t_eng) {
    /*double Gamma0 = jet.Gamma0(theta);
    double dM0 = jet.dEdOmega(theta, t_eng) / (jet.Gamma0(theta) * con::c2);
    double n4 = dM0 / (r * r * D_com / Gamma);
    double n1 = medium.rho(r) / con::mp;
    double gamma_eos = (4 * Gamma + 1) / (3 * Gamma);
    double gamma34 = (Gamma / Gamma0 + Gamma0 / Gamma) / 2;
    double n3 = n4 * (gamma_eos * gamma34 + 1) / (gamma_eos - 1);
    return 1 / (Gamma0 * Gamma * std::sqrt(n4 / n1) * (1 - Gamma0 * n4 / Gamma / n3));*/
    return 0;
};  

void solve_single_shell(Array const& r, Array& Gamma, Array& t_com, Array& D_FS, Array& D_RS, double u0,
                        BlastWaveEqn const& eqn) {
    using namespace boost::numeric::odeint;
    double atol = 0;     // integrator absolute tolerance
    double rtol = 1e-9;  // integrator relative tolerance
    auto stepper = bulirsch_stoer_dense_out<std::vector<double>>{atol, rtol};
    Array state{0, 0, 0, 0, 0, 0};

    double dr = (r[1] - r[0]) / 1000;
    double r0 = r[0];
    double Gamma0 = Gamma[0];
    double beta0 = sqrt(Gamma0 * Gamma0 - 1) / Gamma0;
    double t_com0 = t_com[0];

    // engine time used to calculate the energy injection
    double t_eng0 = r0 * (1 - beta0) / beta0 / con::c;

    // initialize the integrator
    stepper.initialize(Array{Gamma0, u0, t_eng0, t_com0, 0, 0}, r0, dr);
    // integrate the shell over r
    for (int i = 0; stepper.current_time() <= r.back();) {
        stepper.do_step(eqn);

        for (; stepper.current_time() > r[i + 1] && i + 1 < r.size();) {
            i++;
            stepper.calc_state(r[i], state);
            Gamma[i] = state[0];
            // state[1]: u blast wave internal energy
            // state[2]: engine time
            t_com[i] = state[3];
            D_RS[i] = state[4];
            D_FS[i] = state[5];
        }
    }
}

std::pair<Shock, Shock> gen_shocks(Coord const& coord, Jet const& jet, Medium const& medium) {
    Shock shock_f(coord, medium.eps_e, medium.eps_B, medium.xi, medium.zeta);  // forward shock
    Shock shock_r(coord, medium.eps_e, medium.eps_B, medium.xi, medium.zeta);  // reverse shock

    for (size_t i = 0; i < coord.theta.size(); ++i) {
        auto eqn = BlastWaveEqn(medium, jet, coord.theta_b[i], coord.theta_b[i + 1]);
        double Gamma0 = jet.Gamma0(coord.theta[i]);

        shock_f.Gamma[i][0] = Gamma0;
        shock_f.t_com[i][0] = coord.r[0] / sqrt(Gamma0 * Gamma0 - 1) / con::c;

        // shell solid angle
        double dcos = std::fabs(std::cos(coord.theta_b[i + 1]) - std::cos(coord.theta_b[i]));
        double dphi = 2 * con::pi;
        double dOmega = 2 * dphi * dcos;  // bipolar outflow

        // initial internal energy
        double u0 = (Gamma0 - 1) * medium.mass(coord.r[0]) * dOmega / (4 * con::pi) * con::c2;
        solve_single_shell(coord.r, shock_f.Gamma[i], shock_f.t_com[i], shock_f.width[i], shock_r.width[i], u0, eqn);
    }
    FS_n2(shock_f, coord, medium);
    co_moving_B(shock_f, coord);
    co_moving_FS_shock_width(shock_f, coord);

    return std::make_pair(shock_r, shock_f);
}
