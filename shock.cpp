
#include "shock.h"

#include <boost/numeric/odeint.hpp>

#include "afterglow.h"
#include "macros.h"
#include "mesh.h"
#include "physics.h"
#include "utilities.h"

inline double coMovingWeibelB(double eps_B, double e_thermal) { return sqrt(8 * con::pi * eps_B * e_thermal); }

inline double e_ThermalDownStr(double gamma_rel, double n_down_str) {
    return n_down_str * (gamma_rel - 1) * con::mp * con::c2;
}

double u_DownStr(double gamma_rel, double sigma) {
    double ad_idx = adiabaticIndex(gamma_rel);
    if (sigma == 0) {
        return std::sqrt((gamma_rel - 1) * (ad_idx - 1) * (ad_idx - 1) / (ad_idx * (2 - ad_idx) * (gamma_rel - 1) + 2));
    } else {
        double A = ad_idx * (2 - ad_idx) * (gamma_rel - 1) + 2;
        double B = -(gamma_rel + 1) *
                       ((2 - ad_idx) * (ad_idx * gamma_rel * gamma_rel + 1) + ad_idx * (ad_idx - 1) * gamma_rel) *
                       sigma -
                   (gamma_rel - 1) * (ad_idx * (2 - ad_idx) * (gamma_rel * gamma_rel - 2) + 2 * gamma_rel + 3);
        double C = (gamma_rel + 1) * (ad_idx * (1 - ad_idx / 4) * (gamma_rel * gamma_rel - 1) + 1) * sigma * sigma +
                   (gamma_rel * gamma_rel - 1) * (2 * gamma_rel - (2 - ad_idx) * (ad_idx * gamma_rel - 1)) * sigma +
                   (gamma_rel + 1) * (gamma_rel - 1) * (gamma_rel - 1) * (ad_idx - 1) * (ad_idx - 1);
        double D =
            -(gamma_rel - 1) * (gamma_rel + 1) * (gamma_rel + 1) * (2 - ad_idx) * (2 - ad_idx) * sigma * sigma / 4;
        /*double x0 = (-B - sqrt(B * B - 3 * A * C)) / 3 / A;
        double x1 = (-B + sqrt(B * B - 3 * A * C)) / 3 / A;
        return sqrt(
            root_bisection([=](double x) -> double { return A * x * x * x + B * x * x + C * x + D; }, x0, x1, 1e-6));*/
        double b = B / A;
        double c = C / A;
        double d = D / A;
        double P = c - b * b / 3;
        double Q = 2 * b * b * b / 27 - b * c / 3 + d;
        double u = std::sqrt(-P / 3);
        double uds = 2 * u * std::cos(std::acos((3 * Q / (2 * P * u)) - 2 * con::pi) / 3) - b / 3;
        return std::sqrt(uds);
    }
}

double u_UpStr(double u_down, double gamma_rel) {
    return std::sqrt((1 + u_down * u_down) * (gamma_rel * gamma_rel - 1)) + u_down * gamma_rel;
}

double u_UpStr2u_DownStr(double gamma_rel, double sigma) {
    double u_down_s_ = u_DownStr(gamma_rel, sigma);
    double u_up_s_ = u_UpStr(u_down_s_, gamma_rel);
    double ratio_u = u_up_s_ / u_down_s_;
    if (u_down_s_ == 0) {
        ratio_u = (7 * gamma_rel + 1) / (gamma_rel + 1);  // (g_hat+1)/(g_hat-1)
    }
    return ratio_u;
}

double n_DownStr(double n_up_str, double gamma_rel, double sigma) {
    return n_up_str * u_UpStr2u_DownStr(gamma_rel, sigma);
}

double soundSpeed(double pressure, double ad_idx, double rho_rest) {
    return std::sqrt(ad_idx * pressure / (rho_rest * con::c2 + ad_idx / (ad_idx - 1) * pressure)) * con::c;
}

double dtdr_Engine(double Gamma) {
    double Gb = std::sqrt(Gamma * Gamma - 1);
    return std::fabs(Gamma - Gb) / (Gb * con::c);
}

double dtdr_CoMoving(double Gamma) { return 1 / (std::sqrt(Gamma * Gamma - 1) * con::c); };  // co-moving time

double dDdr_Jet(double Gamma) {  // does not applies to initial non-relativistic jet
    double cs = con::c / std::sqrt(3);
    return cs * dtdr_CoMoving(Gamma) / Gamma;
}
//----------------------------------------reverse shock---------------------------------------
double fa(double gamma34, double u3s_, double sigma) {
    return 1 - sigma * (gamma34 + 1) /
                   (u3s_ * u3s_ * gamma34 + u3s_ * std::sqrt((1 + u3s_ * u3s_) * (gamma34 * gamma34 - 1))) / 2;
}

double fc(double p2, double pB3) {
    double p3 = p2 - pB3;
    return p2 / p3;
}

double calc_n4(double dEdOmega, double Gamma0, double r, double D_jet_lab, double sigma) {
    return dEdOmega / (Gamma0 * con::mp * con::c2 * r * r * Gamma0 * D_jet_lab) / (1 + sigma);
}

double calc_pB4(double n4, double sigma) { return sigma * n4 * con::mp * con::c2 / 2; }
//--------------------------------------------------------------------------------------------
Shock::Shock(size_t phi_size, size_t theta_size, size_t r_size, double eps_e, double eps_B)
    : t_com(create3DGrid(phi_size, theta_size, r_size, 0)),
      t_eng(create3DGrid(phi_size, theta_size, r_size, 0)),
      Gamma_rel(create3DGrid(phi_size, theta_size, r_size, 1)),
      B(create3DGrid(phi_size, theta_size, r_size, 0)),
      column_num_den(create3DGrid(phi_size, theta_size, r_size, 0)),
      eps_e(eps_e),
      eps_B(eps_B),
      phi_size(phi_size),
      theta_size(theta_size),
      r_size(r_size) {}

ForwardShockEqn::ForwardShockEqn(Medium const& medium, Ejecta const& jet, Ejecta const& inject, double phi,
                                 double theta, double eps_e)
    : medium(medium),
      jet(jet),
      inject(inject),
      phi(phi),
      theta(theta),
      eps_e(eps_e),
      jet_sigma0(jet.sigma0(phi, theta, 0)),
      inject_sigma0(inject.sigma0(phi, theta, 0)),
      spreading_factor(1),
      jet_Gamma0(jet.Gamma0(phi, theta, 0)),
      inject_Gamma0(inject.Gamma0(phi, theta, 0)),
      gamma4(jet_Gamma0) {};
// dM0dOmega(jet.dE0dOmega(theta) / (jet.Gamma0(theta) * con::c2)) {};

void ForwardShockEqn::operator()(State const& y, State& dydr, double r) {
    double Gamma = y[0];
    double u = y[1];
    double t_eng = y[2];  // engine time
    // double t_com = y[3];  // co-moving time
    // double D_jet = y[4];  // co-moving jet shell width

    dydr[0] = dGammadr(r, Gamma, u, t_eng);
    dydr[1] = dUdr(r, Gamma, u, t_eng);
    dydr[2] = dtdr_Engine(Gamma);
    dydr[3] = dtdr_CoMoving(Gamma);
    dydr[4] = dDdr_Jet(this->gamma4);
}

double ForwardShockEqn::dGammadr(double r, double Gamma, double u, double t_eng) {
    double ad_idx = adiabaticIndex(Gamma);
    double dM = jet.dEdOmega(phi, theta, t_eng) / (jet_Gamma0 * (1 + jet_sigma0) * con::c2);
    double dm = medium.mass(r) / (4 * con::pi);
    double dm_inj = inject.dEdOmega(phi, theta, t_eng) / (inject_Gamma0 * (1 + inject_sigma0) * con::c2);

    double Gamma2 = Gamma * Gamma;

    double a1 = -Gamma * (Gamma2 - 1) * (ad_idx * Gamma - ad_idx + 1) * r * r * medium.rho(r) * con::c2;
    double a2 = Gamma * (ad_idx - 1) * (ad_idx * Gamma2 - ad_idx + 1) * 3 * u / r;
    double L_jet = jet.dLdOmega(phi, theta, t_eng);  // most of case =0;
    double L_inj = inject.dLdOmega(phi, theta, t_eng);
    double a3 = Gamma2 * dtdr_Engine(Gamma) *
                (L_jet * (1 - Gamma / (jet_Gamma0 * (1 + jet_sigma0))) +
                 L_inj * (1 - Gamma / (inject_Gamma0 * (1 + inject_sigma0))));

    double b1 = Gamma2 * (dM + dm + dm_inj) * con::c2;
    double b2 = (ad_idx * ad_idx * (Gamma2 - 1) + 3 * ad_idx - 2) * u;

    return (a1 + a2 + a3) / (b1 + b2);
}

double ForwardShockEqn::dUdr(double r, double Gamma, double u, double t_eng) {
    double ad_idx = adiabaticIndex(Gamma);
    double E = r * r * medium.rho(r) * con::c2;
    return (1 - eps_e) * (Gamma - 1) * E -
           (ad_idx - 1) * (3 / r - dGammadr(r, Gamma, u, t_eng) / Gamma) * u * spreading_factor;
}

bool checkReverseShockGeneration(FRShockEqn const& eqn, double r, double gamma, double t_eng, double D_jet) {
    double n4 = calc_n4(eqn.jet.dEdOmega(eqn.phi, eqn.theta, t_eng), gamma, r, D_jet, eqn.sigma);
    double n1 = eqn.medium.rho(r) / con::mp;
    return eqn.sigma < 8. / 3 * gamma * gamma * n1 / n4;
}

double calc_gamma3(double r, double n1, double n4, double gamma4, double sigma) {
    double C = n4 / n1 * (1 + sigma);
    double gamma3 =
        gamma4 * std::sqrt(((C - 2) - 2 * std::sqrt(C * (gamma4 * gamma4 - 1) + 1)) / (C - 4 * gamma4 * gamma4));
    return gamma3;
}

double FRShockEqn::dN3drPerOmega(double r, double n1, double n4, double gamma3) {
    double gamma34 = (gamma4 / gamma3 + gamma3 / gamma4) / 2;
    double ratio_u = u_UpStr2u_DownStr(gamma34, this->sigma);

    /*double ad_idx2 = adiabatic_index(Gamma);
    double ad_idx3 = adiabatic_index(Gamma34);
    double u3s_ = u_down_str(Gamma34, this->sigma);
    double u4s_ = u_up_str(u3s_, Gamma34);
    double n2 = n_down_str(n1, Gamma, this->sigma);
    double e2 = e_thermal_down_str(Gamma, n2);
    double p2 = (ad_idx2 - 1) * e2;
    double pB4 = calc_pB4(n4, this->sigma);
    double pB3 = pB4 * ratio_u * ratio_u;
    double f_a = fa(Gamma34, u3s_, this->sigma);
    double f_b = ratio_u / ((ad_idx3 * Gamma34 + 1) / (ad_idx3 - 1));
    double f_c = fc(p2, pB3);
    double F = f_a * f_b * f_c;*/
    double n3 = n4 * ratio_u;
    double dxdr = 1. / (gamma4 * std::sqrt((1 + this->sigma) * n4 / n1) * std::fabs(1 - gamma4 * n4 / gamma3 / n3));
    return n3 * r * r * gamma3 * dxdr;
}

FRShockEqn::FRShockEqn(Medium const& medium, Ejecta const& jet, Ejecta const& inject, double phi, double theta)
    : medium(medium),
      jet(jet),
      inject(inject),
      phi(phi),
      theta(theta),
      sigma(jet.sigma0(phi, theta, 0)),
      gamma4(jet.Gamma0(phi, theta, 0)) {}

void FRShockEqn::operator()(State const& y, State& dydr, double r) {
    // y[0] leave blank
    // double D_rs = y[1];
    double t_eng = y[2];
    // double t_com = y[3];
    double D_jet_lab = y[4];

    double n4 = calc_n4(jet.dEdOmega(phi, theta, t_eng), gamma4, r, D_jet_lab, this->sigma);
    double n1 = medium.rho(r) / con::mp;

    double gamma3 = calc_gamma3(r, n1, n4, this->gamma4, this->sigma);
    dydr[0] = 0;
    dydr[1] = dN3drPerOmega(r, n1, n4, gamma3);
    dydr[2] = dtdr_Engine(gamma3);
    dydr[3] = dtdr_CoMoving(gamma3);
    dydr[4] = dDdr_Jet(this->gamma4);
}

void updateShockState(Shock& shock, size_t i, size_t j, size_t k, double r, double Gamma_rel, double t_com,
                      double t_eng, double dMdOmega_up, double n_up_str, double sigma) {
    if (Gamma_rel > 1) {
        double ratio_u = u_UpStr2u_DownStr(Gamma_rel, sigma);
        double pB_up = calc_pB4(n_up_str, sigma);
        double pB_down = pB_up * ratio_u * ratio_u;
        double n_down_str = n_up_str * ratio_u;
        double co_moving_width = dMdOmega_up / (r * r * n_down_str * con::mp);
        double e_th = e_ThermalDownStr(Gamma_rel, n_down_str);
        shock.Gamma_rel[i][j][k] = Gamma_rel;
        shock.t_com[i][j][k] = t_com;
        shock.t_eng[i][j][k] = t_eng;
        shock.column_num_den[i][j][k] = n_down_str * co_moving_width;
        shock.B[i][j][k] = coMovingWeibelB(shock.eps_B, e_th) + std::sqrt(pB_down * 8 * con::pi);
    } else {
        shock.Gamma_rel[i][j][k] = 1;
        shock.t_com[i][j][k] = 0;
        shock.t_eng[i][j][k] = 0;
        shock.column_num_den[i][j][k] = 0;
        shock.B[i][j][k] = 0;
    }
}

struct CrossState {
    double gamma_rel;
    double r;
    double column_num_den;
    double B;
};

void Blandford_McKee(size_t i, size_t j, size_t k, Shock& shock, CrossState const& state_c, double r, double t_com,
                     double t_eng) {
    double const g = 2.;
    shock.t_com[i][j][k] = t_com;
    shock.t_eng[i][j][k] = t_eng;

    // shock.Gamma[j][k] = (state_c.gamma - 1) * std::pow(r / state_c.r, -g) + 1;
    //  shock.n_p[j][k] = state_c.n3 * std::pow(r / state_c.r, -6 * (3 + g) / 7);
    // shock.e_th[j][k] = state_c.e3 * std::pow(r / state_c.r, -8 * (3 + g) / 7);
    // shock.width_eff[j][k] = state_c.D_eff * std::pow(r / state_c.r, (6. * (3. + g) - 14.) / 7.);

    shock.Gamma_rel[i][j][k] = (state_c.gamma_rel - 1) * std::pow(r / state_c.r, -g) + 1;
    shock.column_num_den[i][j][k] = state_c.column_num_den * std::pow(r / state_c.r, -2);
    shock.B[i][j][k] = state_c.B * std::pow(r / state_c.r, (6. * (3. + g) - 14.) / 7.);
}

void setForwardInit(ForwardShockEqn& eqn, ForwardShockEqn::State& state, double r0) {
    double gamma2 = eqn.jet.Gamma0(eqn.phi, eqn.theta, 0);
    double u0 = (gamma2 - 1) * eqn.medium.mass(r0) / (4 * con::pi) * con::c2;
    double beta0 = gammaTobeta(gamma2);
    double t_eng0 = r0 * (1 - beta0) / beta0 / con::c;
    double t_com0 = r0 / std::sqrt(gamma2 * gamma2 - 1) / con::c;
    double D_jet0 = con::c * eqn.jet.duration;
    state = {gamma2, u0, t_eng0, t_com0, D_jet0};
}

void setReverseInit(FRShockEqn& eqn, FRShockEqn::State& state, double r0) {
    double gamma4 = eqn.jet.Gamma0(eqn.phi, eqn.theta, 0);
    double u0 = (gamma4 - 1) * eqn.medium.mass(r0) / (4 * con::pi) * con::c2;
    double beta0 = gammaTobeta(gamma4);
    double t_eng0 = r0 * (1 - beta0) / beta0 / con::c;
    double t_com0 = r0 / std::sqrt(gamma4 * gamma4 - 1) / con::c;
    double D_jet0 = con::c * eqn.jet.duration;
    double dN3dOmega = 0;
    state = {0., dN3dOmega, t_eng0, t_com0, D_jet0};
}

void updateForwardShock(size_t i, size_t j, int k, double r_k, ForwardShockEqn& eqn, const FRShockEqn::State& state,
                        Shock& f_shock) {
    double n1 = eqn.medium.rho(r_k) / con::mp;
    double Gamma = state[0];
    double t_eng = state[2];
    double t_com = state[3];
    double dM1dOmega = eqn.medium.mass(r_k) / (4 * con::pi);

    updateShockState(f_shock, i, j, k, r_k, Gamma, t_com, t_eng, dM1dOmega, n1, eqn.jet_sigma0);
}

void solveForwardShell(size_t i, size_t j, const Array& r_b, const Array& r, Shock& f_shock, ForwardShockEqn& eqn) {
    using namespace boost::numeric::odeint;

    double atol = 0, rtol = 1e-6, r0 = r_b[0];
    double dr = (r[1] - r[0]) / 100;
    ForwardShockEqn::State state;
    setForwardInit(eqn, state, r0);
    if (state[0] <= con::Gamma_cut) {  // initial low Lorentz factor
        return;
    }

    auto stepper = bulirsch_stoer_dense_out<ForwardShockEqn::State>{atol, rtol};
    stepper.initialize(state, r0, dr);

    double t_com_last = state[3];

    double r_back = r[r.size() - 1];

    for (int k = 0; stepper.current_time() <= r_back;) {
        stepper.do_step(eqn);

        while (k < r.size() && stepper.current_time() > r[k]) {
            stepper.calc_state(r[k], state);
            updateForwardShock(i, j, k, r[k], eqn, state, f_shock);
            ++k;
        }
    }
}

void solveFRShell(size_t i, size_t j, Array const& r_b, Array const& r, Shock& f_shock, Shock& r_shock,
                  ForwardShockEqn& eqn_f, FRShockEqn& eqn_r) {
    using namespace boost::numeric::odeint;
    double atol = 0, rtol = 1e-6, r0 = r_b[0];
    double dr = (r[1] - r[0]) / 100;
    std::array<double, 5> state;

    auto stepper = bulirsch_stoer_dense_out<ForwardShockEqn::State>{atol, rtol};
    // auto stepper = make_dense_output(atol, rtol, runge_kutta_dopri5<std::vector<double>>());

    setForwardInit(eqn_f, state, r0);
    if (state[0] <= 1 + 1e-6) {  // initial low Lorentz factor
        return;
    }
    double gamma3 = eqn_f.gamma4;
    double t_eng = state[2];
    double t_com = state[3];
    double D_jet = state[4];

    bool RS_crossing = checkReverseShockGeneration(eqn_r, r0, gamma3, t_eng, D_jet);
    bool crossed = false;

    if (RS_crossing) {
        setReverseInit(eqn_r, state, r0);
    }

    double t_com_last = state[3];
    double n1 = 0, n4 = 0;
    double dN3dOmega = 0;
    double dN4dOmega =
        eqn_r.jet.dEdOmega(eqn_r.phi, eqn_r.theta, t_eng) / (eqn_r.gamma4 * con::mp * con::c2 * (1 + eqn_r.sigma));

    CrossState state_c;
    stepper.initialize(state, r0, dr);
    // integrate the shell over r
    double r_back = r[r.size() - 1];
    for (int k = 0; stepper.current_time() <= r_back;) {
        RS_crossing ? stepper.do_step(eqn_r) : stepper.do_step(eqn_f);

        for (; stepper.current_time() > r[k] && k < r.size(); k++) {
            n1 = eqn_f.medium.rho(r[k]) / con::mp;
            stepper.calc_state(r[k], state);

            t_eng = state[2];
            t_com = state[3];
            D_jet = state[4];
            if (RS_crossing) {
                n4 = calc_n4(eqn_r.jet.dEdOmega(eqn_r.phi, eqn_r.theta, t_eng), eqn_r.gamma4, r[k], D_jet, eqn_r.sigma);
                gamma3 = calc_gamma3(r[k], n1, n4, eqn_r.gamma4, eqn_r.sigma);
            } else {
                gamma3 = state[0];
            }

            double dM2dOmega = eqn_f.medium.mass(r[k]) / (4 * con::pi);
            updateShockState(f_shock, i, j, k, r[k], gamma3, t_com, t_eng, dM2dOmega, n1, eqn_f.jet_sigma0);

            if (!crossed && RS_crossing) {
                dN3dOmega = state[1];
                double gamma34 = (eqn_r.gamma4 / gamma3 + gamma3 / eqn_r.gamma4) / 2;
                double dM3dOmega = dN3dOmega * con::mp * con::c2;
                updateShockState(r_shock, i, j, k, r[k], gamma34, t_com, t_eng, dM3dOmega, n4, eqn_r.sigma);

                crossed = dN3dOmega >= dN4dOmega;
                if (crossed) {
                    RS_crossing = false;
                    state_c = {r_shock.Gamma_rel[i][j][k], r[k], r_shock.column_num_den[i][j][k], r_shock.B[i][j][k]};
                    double u0 = (gamma3 - 1) * eqn_r.medium.mass(r[k]) / (4 * con::pi) * con::c2;
                    state = {gamma3, u0, t_eng, t_com, D_jet};
                    eqn_f.gamma4 = gamma3;
                    stepper.initialize(state, r[k], dr);
                }
            } else if (!crossed && !RS_crossing) {
                RS_crossing = checkReverseShockGeneration(eqn_r, r0, gamma3, t_eng, D_jet);
                if (RS_crossing) {
                    state = {0., 0., t_eng, t_com, D_jet};
                    stepper.initialize(state, r[k], dr);
                }
            } else {
                Blandford_McKee(i, j, k, r_shock, state_c, r[k], t_com, t_eng);
            }
        }
    }
}

/*
void solve_single_shell(size_t j, Array const& r_b, Array const& r, Shock& f_shock, Shock& r_shock,
                        ForwardShockEqn& eqn) {
    // integrate the shell over r
    for (int k = 0, k1 = 1; stepper.current_time() <= r_b.back();) {
        stepper.do_step(eqn);

        if (eqn.jet.spreading) {  // check if this works for sigma!=0. should be sound speed in region 3;
            double dr = stepper.current_time_step();
            double r_last = stepper.current_time() - dr;
            double Gamma_axis = loglog_interp(r_last, r, f_shock.Gamma[0]);
            if (k == 0) {
                Gamma_axis = eqn.jet.Gamma0_profile(0);
            }

            double ad_idx = adiabatic_index(Gamma_axis);
            double n1 = eqn.medium.rho(r_last) / con::mp;
            double n2 = n_down_str(n1, Gamma_axis, 0);
            double p2 = (ad_idx - 1) * e_thermal_down_str(Gamma_axis, n2);
            double jet_cs = sound_speed(p2, ad_idx, n2 * con::mp);
            if (theta_c < con::pi / 2) {
                theta_c += jet_cs / con::c * dr / (r_last * Gamma_axis * Gamma_axis * theta_c);
                // theta_c += jet_cs / con::c * dr / (r_last * Gamma_axis);
            }
            double t_eng = state[2];
            double dEdOmega_init = eqn.jet.dEdOmega(eqn.theta, t_eng);
            double dEdOmega_spread = eqn.jet.dEdOmega_spread(eqn.theta, theta_c, t_eng);

            if (dEdOmega_init == 0) {
                eqn.spreading_factor = 1;
            } else {
                eqn.spreading_factor = dEdOmega_spread / dEdOmega_init;
            }
        }
    }
}*/

Shock genForwardShock(Coord const& coord, Ejecta const& jet, Medium const& medium, double eps_e, double eps_B,
                      Ejecta const& inject) {
    auto [phi_size, theta_size, r_size] = coord.shape();
    Shock f_shock(1, theta_size, r_size, eps_e, eps_B);

    for (size_t j = 0; j < theta_size; ++j) {
        auto eqn = ForwardShockEqn(medium, jet, inject, 0, coord.theta[j], eps_e);
        solveForwardShell(0, j, coord.r_b, coord.r, f_shock, eqn);
    }

    return f_shock;
}

std::pair<Shock, Shock> genFRShocks(Coord const& coord, Ejecta const& jet, Medium const& medium, double eps_e,
                                    double eps_B, Ejecta const& inject) {
    auto [phi_size, theta_size, r_size] = coord.shape();

    Shock f_shock(1, theta_size, r_size, eps_e, eps_B);
    Shock r_shock(1, theta_size, r_size, eps_e, eps_B);

    for (size_t j = 0; j < theta_size; ++j) {
        auto eqn_f = ForwardShockEqn(medium, jet, inject, 0, coord.theta[j], eps_e);
        auto eqn_r = FRShockEqn(medium, jet, inject, 0, coord.theta[j]);
        solveFRShell(0, j, coord.r_b, coord.r, f_shock, r_shock, eqn_f, eqn_r);
    }

    return std::make_pair(std::move(f_shock), std::move(r_shock));
}

Shock genForwardShock3D(Coord const& coord, Ejecta const& jet, Medium const& medium, double eps_e, double eps_B,
                        Ejecta const& inject) {
    auto [phi_size, theta_size, r_size] = coord.shape();

    Shock f_shock(phi_size, theta_size, r_size, eps_e, eps_B);
    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            auto eqn = ForwardShockEqn(medium, jet, inject, coord.phi[i], coord.theta[j], eps_e);
            solveForwardShell(i, j, coord.r_b, coord.r, f_shock, eqn);
        }
    }
    return f_shock;
}

std::pair<Shock, Shock> genFRShocks3D(Coord const& coord, Ejecta const& jet, Medium const& medium, double eps_e,
                                      double eps_B, Ejecta const& inject) {
    auto [phi_size, theta_size, r_size] = coord.shape();

    Shock f_shock(phi_size, theta_size, r_size, eps_e, eps_B);
    Shock r_shock(phi_size, theta_size, r_size, eps_e, eps_B);
    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            auto eqn_f = ForwardShockEqn(medium, jet, inject, coord.phi[i], coord.theta[j], eps_e);
            auto eqn_r = FRShockEqn(medium, jet, inject, coord.phi[i], coord.theta[j]);
            solveFRShell(i, j, coord.r_b, coord.r, f_shock, r_shock, eqn_f, eqn_r);
        }
    }
    return std::make_pair(std::move(f_shock), std::move(r_shock));
}
