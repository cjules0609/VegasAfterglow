
#include "shock.h"

#include <boost/numeric/odeint.hpp>

#include "macros.h"
#include "physics.h"
Shock::Shock(Coord const& coord, double eps_e, double eps_B, double xi, double zeta)
    : t_com(create_grid(coord.theta.size(), coord.r.size(), 0)),
      t_com_b(create_grid(coord.theta.size(), coord.r_b.size(), 0)),
      Gamma(create_grid(coord.theta.size(), coord.r.size(), 1)),
      e_th(create_grid(coord.theta.size(), coord.r.size(), 1)),
      B(create_grid(coord.theta.size(), coord.r.size(), 0)),
      width(create_grid(coord.theta.size(), coord.r.size(), 0)),
      n_p(create_grid(coord.theta.size(), coord.r.size(), 0)),
      eps_e{eps_e},
      eps_B{eps_B},
      xi{xi},
      zeta{zeta} {}

double co_moving_B(double eps_B, double e_thermal) { return sqrt(8 * con::pi * eps_B * e_thermal); }

double n_down_str(double Gamma_relative, double n_up_str, double cs_up_str) {
    double u4 = sqrt(Gamma_relative * Gamma_relative - 1);
    double us4 = 4 * u4 * sqrt((1 + u4 * u4) / (8 * u4 * u4 + 9));  // check this
    double M = us4 / cs_up_str;
    if (M > 10) {
        return 4 * Gamma_relative * n_up_str;
    } else if (M >= 1) {
        double ad_idx = adiabatic_index(Gamma_relative);
        return n_up_str * (ad_idx + 1) * M * M / ((ad_idx - 1) * M * M + 2);
    } else {
        return 0;  // shock cannot be formed
    }
}

inline double e_thermal_down_str(double Gamma_relative, double n_down_str) {
    return n_down_str * (Gamma_relative - 1) * con::mp * con::c2;
}

double sound_speed(double pressure, double ad_idx, double rho_rest) {
    return sqrt(ad_idx * pressure / (rho_rest * con::c2 + ad_idx / (ad_idx - 1) * pressure)) * con::c;
}

BlastWaveEqn::BlastWaveEqn(Medium const& medium, Jet const& blast, double theta_lo, double theta_hi)
    : medium(medium),
      jet(blast),
      theta_lo(theta_lo),
      theta_hi(theta_hi),
      theta(0.5 * (theta_lo + theta_hi)),
      dOmega(2 * con::pi * std::fabs(std::cos(theta_hi) - std::cos(theta_lo))){};

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
    dydr[5] = dDdr_FS(r, Gamma);
};

double BlastWaveEqn::dGammadr(double r, double Gamma, double u, double t_eng) {
    double ad_idx = adiabatic_index(Gamma);
    double dm = medium.mass(r) * dOmega / (4 * con::pi);
    double dM0 = jet.dEdOmega(theta, t_eng) * dOmega / (jet.Gamma0_profile(theta) * con::c2);
    double Gamma2 = Gamma * Gamma;
    double a1 = dOmega * r * r * medium.rho(r) / Gamma * (Gamma2 - 1) * (ad_idx * Gamma - ad_idx + 1);
    double a2 = -(ad_idx - 1) / Gamma * (ad_idx * Gamma2 - ad_idx + 1) * 3 * u / r;
    double b1 = (dM0 + dm) * con::c2;
    double b2 = (ad_idx * ad_idx * (Gamma2 - 1) + 3 * ad_idx - 2) * u / Gamma2;
    return -(a1 + a2) / (b1 + b2);
};

double BlastWaveEqn::dUdr(double r, double Gamma, double u, double t_eng) {
    double ad_idx = adiabatic_index(Gamma);
    double E = dOmega * r * r * medium.rho(r) * con::c2;
    return (1 - medium.eps_e) * (Gamma - 1) * E - (ad_idx - 1) * (3 / r - dGammadr(r, Gamma, u, t_eng) / Gamma) * u;
};

double BlastWaveEqn::dtdr_eng(double Gamma) {
    double Gb = std::sqrt(Gamma * Gamma - 1);
    return (Gamma - Gb) / (Gb * con::c);
};

double BlastWaveEqn::dtdr_com(double Gamma) { return 1 / (std::sqrt(Gamma * Gamma - 1) * con::c); };  // co-moving time

double BlastWaveEqn::dDdr_FS(double r, double Gamma) {
    double ad_idx = adiabatic_index(Gamma);
    double n1 = medium.rho(r) / con::mp;
    double n2 = n_down_str(Gamma, n1, medium.cs);
    double p2 = (ad_idx - 1) * e_thermal_down_str(Gamma, n2);
    double cs = sound_speed(p2, ad_idx, n2 * con::mp);
    return cs * dtdr_com(Gamma);
};

double BlastWaveEqn::dDdr_RS(double r, double Gamma, double D_com, double t_eng) {
    double Gamma0 = jet.Gamma0_profile(theta);
    double n4 = jet.dEdOmega(theta, t_eng) / (Gamma0 * con::mp * con::c2 * r * r * D_com);
    double n1 = medium.rho(r) / con::mp;
    // double gamma34 = sqrt(Gamma0 / (2 * sqrt(n4 / n1)));
    double Gamma34 = (Gamma0 / Gamma + Gamma / Gamma0) / 2;
    double ad_idx3 = adiabatic_index(Gamma34);
    double n3 = n4 * (ad_idx3 * Gamma34 + 1) / (ad_idx3 - 1);

    return Gamma / (Gamma0 * std::sqrt(n4 / n1) * (1 - Gamma0 * n4 / Gamma / n3));
};

void solve_single_shell(size_t j, Array const& r_b, Array const& r, Shock& shock_r, Shock& shock_f,
                        BlastWaveEqn const& eqn) {
    using namespace boost::numeric::odeint;
    double atol = 0;     // integrator absolute tolerance
    double rtol = 1e-6;  // integrator relative tolerance
    // auto stepper = bulirsch_stoer_dense_out<std::vector<double>>{atol, rtol};
    auto stepper = make_dense_output(atol, rtol, runge_kutta_dopri5<std::vector<double>>());

    Array state{0, 0, 0, 0, 0, 0};

    double r0 = r[0];
    double Gamma0 = shock_f.Gamma[j][0];
    double beta0 = gamma_to_beta(Gamma0);
    double t_eng0 = r0 * (1 - beta0) / beta0 / con::c;
    double u0 = (Gamma0 - 1) * eqn.medium.mass(r0) * eqn.dOmega / (4 * con::pi) * con::c2;
    double t_com0 = shock_f.t_com[j][0];
    double FS_width0 = shock_f.width[j][0];
    double RS_width0 = shock_r.width[j][0];

    // initialize the integrator
    double dr = (r[1] - r[0]) / 1000;
    stepper.initialize(Array{Gamma0, u0, t_eng0, t_com0, RS_width0, FS_width0}, r0, dr);
    // integrate the shell over r

    double N3 = 0;
    double E_th = 0;
    bool RS_crossed = false;
    for (int k = 0, k1 = 0; stepper.current_time() <= r.back();) {
        stepper.do_step(eqn);

        for (; stepper.current_time() > r[k + 1] && k + 1 < r.size();) {
            k++;
            stepper.calc_state(r[k], state);
            double Gamma = state[0];
            // state[1]: u blast wave internal energy
            double t_eng = state[2];
            double t_com = state[3];
            double D_RS = state[4];
            double D_FS = state[5];

            double n1 = eqn.medium.rho(r[k]) / con::mp;
            double Gamma12 = Gamma;
            shock_f.Gamma[j][k] = Gamma;

            shock_f.t_com[j][k] = t_com;
            // shock_f.width[theta_idx][i] = r[i] / Gamma / 12;  // D_FS;
            shock_f.width[j][k] = D_FS;
            shock_f.n_p[j][k] = n_down_str(Gamma12, n1, eqn.medium.cs);
            shock_f.e_th[j][k] = e_thermal_down_str(Gamma12, shock_f.n_p[j][k]);
            shock_f.B[j][k] = co_moving_B(shock_f.eps_B, shock_f.e_th[j][k]);

            shock_r.Gamma[j][k] = Gamma;
            shock_r.t_com[j][k] = t_com;

            // reverse shock crossing
            if (D_RS < D_FS && !RS_crossed) {
                double Gamma34 = (Gamma0 / Gamma + Gamma / Gamma0) / 2;
                shock_r.width[j][k] = D_RS;
                double n4 = eqn.jet.dEdOmega(eqn.theta, t_eng) / (Gamma0 * con::mp * con::c2 * r[k] * r[k] * D_FS);
                shock_r.n_p[j][k] = n4 * 4 * Gamma34;  // n_down(gamma34, n4, cs_4); can be improved
                shock_r.e_th[j][k] = e_thermal_down_str(Gamma34, shock_r.n_p[j][k]);

            } else {  // reverse shock crossed
                if (!RS_crossed) {
                    N3 = shock_r.n_p[j][k] * D_RS * r[k] * r[k];
                    E_th = shock_r.e_th[j][k] * D_RS * r[k] * r[k];
                    RS_crossed = true;
                }
                shock_r.width[j][k] = D_FS;
                shock_r.n_p[j][k] = N3 / (D_FS * r[k] * r[k]);
                shock_r.e_th[j][k] = E_th / (D_FS * r[k] * r[k]);
            }
            shock_r.B[j][k] = co_moving_B(shock_r.eps_B, shock_r.e_th[j][k]);
        }

        for (; stepper.current_time() > r_b[k1 + 1] && k1 + 1 < r_b.size();) {
            k1++;
            stepper.calc_state(r_b[k1], state);
            shock_f.t_com_b[j][k1] = shock_r.t_com_b[j][k1] = state[3];
        }
    }
}

std::pair<Shock, Shock> gen_shocks(Coord const& coord, Jet const& jet, Medium const& medium) {
    Shock shock_f(coord, medium.eps_e, medium.eps_B, medium.xi, medium.zeta);  // forward shock
    Shock shock_r(coord, medium.eps_e, medium.eps_B, medium.xi, medium.zeta);  // reverse shock

    for (size_t i = 0; i < coord.theta.size(); ++i) {
        auto eqn = BlastWaveEqn(medium, jet, coord.theta_b[i], coord.theta_b[i + 1]);
        double Gamma0 = jet.Gamma0_profile(coord.theta[i]);

        double n1 = medium.rho(coord.r[0]) / con::mp;

        shock_f.Gamma[i][0] = Gamma0;

        shock_f.t_com[i][0] = coord.r[0] / sqrt(Gamma0 * Gamma0 - 1) / con::c;
        shock_f.t_com_b[i][0] = coord.r_b[0] / sqrt(Gamma0 * Gamma0 - 1) / con::c;
        shock_f.width[i][0] = con::c * jet.duration * Gamma0;
        shock_f.n_p[i][0] = n_down_str(Gamma0, n1, medium.cs);
        shock_f.e_th[i][0] = e_thermal_down_str(Gamma0, shock_f.n_p[i][0]);
        shock_f.B[i][0] = co_moving_B(shock_f.eps_B, shock_f.e_th[i][0]);

        shock_r.Gamma[i][0] = shock_f.Gamma[i][0];
        shock_r.t_com[i][0] = shock_f.t_com[i][0];
        shock_r.t_com_b[i][0] = shock_f.t_com_b[i][0];
        shock_r.width[i][0] = 0;
        shock_r.n_p[i][0] = 0;
        shock_r.e_th[i][0] = 0;
        shock_r.B[i][0] = 0;

        solve_single_shell(i, coord.r_b, coord.r, shock_r, shock_f, eqn);
    }
    return std::make_pair(shock_r, shock_f);
}
