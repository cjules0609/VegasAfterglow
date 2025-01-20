
#include "synchrotron.h"

#include <cmath>

#include "afterglow.h"
#include "inverse-compton.h"
#include "macros.h"
#include "physics.h"
#include "utilities.h"
Y_IC::Y_IC(double nu_m, double nu_c, double B, double Y_T) {
    gamma_hat_m = con::me * con::c2 / (con::h * nu_m);
    gamma_hat_c = con::me * con::c2 / (con::h * nu_c);
    this->Y_T = Y_T;
    nu_hat_m = syn_nu(gamma_hat_m, B);
    nu_hat_c = syn_nu(gamma_hat_c, B);

    if (nu_hat_m <= nu_hat_c) {
        regime = 1;  // fast IC cooling
    } else {
        regime = 2;  // slow IC cooling
    }
}

Y_IC::Y_IC(double Y_T) {
    this->Y_T = Y_T;
    regime = 3;
}

Y_IC::Y_IC() {
    nu_hat_m = 0;
    nu_hat_c = 0;
    gamma_hat_m = 0;
    gamma_hat_c = 0;
    Y_T = 0;
    regime = 0;
}

inline double fast_exp(double x) {
    /*constexpr double a = (1ll << 52) / 0.6931471805599453;
    constexpr double b = (1ll << 52) * (1023 - 0.04367744890362246);
    x = a * x + b;

    // Remove these lines if bounds checking is not needed
    // constexpr double c = (1ll << 52);
    // constexpr double d = (1ll << 52) * 2047;
    // if (x < c || x > d)
    //    x = (x < c) ? 0.0 : d;

    // With C++20 one can use std::bit_cast instead
    uint64_t n = static_cast<uint64_t>(x);
    memcpy(&x, &n, 8);

    return x;*/
    return exp(x);
}

inline double fast_pow(double a, double b) {
    uint64_t bits = std::bit_cast<uint64_t>(a);

    uint64_t exponent = static_cast<uint64_t>(b * ((static_cast<int64_t>((bits >> 52) & 0x7FF) - 1023)) + 1023);

    bits = (exponent << 52);

    return std::bit_cast<double>(bits);
}
/*
inline double fast_pow(double a, double b) {
    union {
        double d;
        int x[2];
    } u = {a};
    u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
    u.x[0] = 0;
    return u.d;
    // return pow(a, b);
}*/

inline double pow52(double a) { return sqrt(a * a * a * a * a); }

inline double pow43(double a) { return cbrt(a * a * a * a); }

inline double pow23(double a) { return cbrt(a * a); }

double Y_IC::as_gamma(double gamma, double p) const {
    switch (regime) {
        case 1:
            if (gamma <= gamma_hat_m) {
                return Y_T;
            } else if (gamma <= gamma_hat_c) {
                return Y_T / sqrt(gamma / gamma_hat_m);
            } else {
                return Y_T * pow43(gamma_hat_c / gamma) / sqrt(gamma_hat_c / gamma_hat_m);
            }
            break;
        case 2:
            if (gamma <= gamma_hat_c) {
                return Y_T;
            } else if (gamma <= gamma_hat_m) {
                return Y_T * fast_pow(gamma / gamma_hat_c, (p - 3) / 2);
            } else {
                return Y_T * pow43(gamma_hat_m / gamma) * fast_pow(gamma_hat_m / gamma_hat_c, (p - 3) / 2);
            }
            break;
        case 3:
            return Y_T;
            break;
        default:
            return 0;
            break;
    }
}

double Y_IC::as_nu(double nu, double p) const {
    switch (regime) {
        case 1:
            if (nu <= nu_hat_m) {
                return Y_T;
            } else if (nu <= nu_hat_c) {
                return Y_T * fast_pow(nu / nu_hat_m, -0.25);
            } else {
                return Y_T * pow23(nu_hat_c / nu) * fast_pow(nu / nu_hat_m, -0.25);
            }
            break;
        case 2:
            if (nu <= nu_hat_c) {
                return Y_T;
            } else if (nu <= nu_hat_m) {
                return Y_T * fast_pow(nu / nu_hat_c, (p - 3) / 4);
            } else {
                return Y_T * pow23(nu_hat_m / nu) * fast_pow(nu_hat_m / nu_hat_c, (p - 3) / 4);
            }
            break;
        case 3:
            return Y_T;
            break;
        default:
            return 0;
            break;
    }
}

double Y_IC::Y_T_tilt(std::vector<Y_IC> const& Ys) {
    double Y_tilt = 0;
    for (auto& Y : Ys) {
        Y_tilt += Y.Y_T;
    }
    return Y_tilt;
}

double Y_IC::Y_tilt_at_gamma(std::vector<Y_IC> const& Ys, double gamma, double p) {
    double Y_tilt = 0;
    for (auto& Y : Ys) {
        Y_tilt += Y.as_gamma(gamma, p);
    }
    return Y_tilt;
}

double Y_IC::Y_tilt_at_nu(std::vector<Y_IC> const& Ys, double nu, double p) {
    double Y_tilt = 0;
    for (auto& Y : Ys) {
        Y_tilt += Y.as_nu(nu, p);
    }
    return Y_tilt;
}

SynPhotonsMesh create_syn_photons_grid(size_t theta_size, size_t r_size) {
    SynPhotonsMesh grid(boost::extents[theta_size][r_size]);
    return grid;
}

SynElectronsMesh create_syn_electrons_grid(size_t theta_size, size_t r_size) {
    SynElectronsMesh grid(boost::extents[theta_size][r_size]);
    return grid;
}

double SynElectrons::N(double gamma) const {
    if (gamma < gamma_c) {
        return N_tot * gamma_spectrum_(gamma);
    } else {
        return N_tot * gamma_spectrum_(gamma) * (1 + Y_c) / (1 + Y_IC::Y_tilt_at_gamma(Ys, gamma, p));
    }
}

double SynElectrons::n(double gamma) const {
    if (gamma < gamma_c) {
        return n_tot * gamma_spectrum_(gamma);
    } else {
        return n_tot * gamma_spectrum_(gamma) * (1 + Y_c) / (1 + Y_IC::Y_tilt_at_gamma(Ys, gamma, p));
    }
}

inline bool order(double a, double b, double c) { return a <= b && b <= c; };

size_t get_regime(double a, double c, double m) {
    if (order(a, m, c)) {
        return 1;
    } else if (order(m, a, c)) {
        return 2;
    } else if (order(a, c, m)) {
        return 3;
    } else if (order(c, a, m)) {
        return 4;
    } else if (order(m, c, a)) {
        return 5;
    } else if (order(c, m, a)) {
        return 6;
    }
    return 0;
}

double SynElectrons::gamma_spectrum_(double gamma) const {
    switch (regime) {
        case 1:  // same as case 2
        case 2:
            if (gamma <= gamma_m) {
                return 0;
            } else if (gamma <= gamma_c) {
                return (p - 1) * fast_pow(gamma / gamma_m, -p) / gamma_m;
            } else {
                return (p - 1) * fast_pow(gamma / gamma_m, -p) * gamma_c / (gamma * gamma_m) *
                       fast_exp(-gamma / gamma_M);
            }
            break;
        case 3:
            if (gamma <= gamma_c) {
                return 0;
            } else if (gamma <= gamma_m) {
                return gamma_c / (gamma * gamma);
            } else {
                return gamma_c / (gamma * gamma_m) * fast_pow(gamma / gamma_m, -p) * fast_exp(-gamma / gamma_M);
            }
            break;
        case 4:  // Gao, Lei, Wu and Zhang 2013 Eq 18
            if (gamma <= gamma_a) {
                return 3 * gamma * gamma / (gamma_a * gamma_a * gamma_a);
            } else if (gamma <= gamma_m) {
                return gamma_c / (gamma * gamma);
            } else {
                return gamma_c / (gamma * gamma_m) * fast_pow(gamma / gamma_m, -p) * fast_exp(-gamma / gamma_M);
            }
            break;
        case 5:  // Gao, Lei, Wu and Zhang 2013 Eq 19
            if (gamma <= gamma_a) {
                return 3 * gamma * gamma / (gamma_a * gamma_a * gamma_a);
            } else {
                return (p - 1) * gamma_c / (gamma * gamma_m) * fast_pow(gamma / gamma_m, -p) *
                       fast_exp(-gamma / gamma_M);
            }
            break;
        case 6:  // Gao, Lei, Wu and Zhang 2013 Eq 20
            if (gamma <= gamma_a) {
                return 3 * gamma * gamma / (gamma_a * gamma_a * gamma_a);
            } else {
                return fast_pow(gamma_m, p - 1) * gamma_c * fast_pow(gamma, -(p + 1)) * fast_exp(-gamma / gamma_M);
            }
            break;
        default:
            return 0;
    }
}

double SynPhotons::L_nu(double nu) const {
    if (nu < nu_c) {
        return L_nu_peak * spectrum_(nu);
    } else {
        return L_nu_peak * spectrum_(nu) * (1 + Y_c) / (1 + Y_IC::Y_tilt_at_nu(Ys, nu, p));
    }
}

double SynPhotons::E_nu(double nu) const {
    if (nu < nu_c) {
        return E_nu_peak * spectrum_(nu);
    } else {
        return E_nu_peak * spectrum_(nu) * (1 + Y_c) / (1 + Y_IC::Y_tilt_at_nu(Ys, nu, p));
    }
}

void SynPhotons::update_constant() {
    // Update constants based on spectral parameters
    a_m_1_3 = cbrt(nu_a / nu_m);                       // a_m_1_3 represents (nu_a / nu_m)^(1/3)
    c_m_1_2 = sqrt(nu_c / nu_m);                       // c_m_1_2 represents (nu_c / nu_m)^(1/2)
    m_a_pa4_2 = fast_pow(nu_m / nu_a, (p + 4) / 2);    // m_a_pa4_2 represents (nu_m / nu_a)^((p+4)/2)
    a_m_mpa1_2 = fast_pow(nu_a / nu_m, (-p + 1) / 2);  // a_m_mpa1_2 represents (nu_a / nu_m)^((-p+1)/2)
    a_c_1_3 = cbrt(nu_a / nu_c);                       // a_c_1_3 represents (nu_a / nu_c)^(1/3)
    a_m_1_2 = sqrt(nu_a / nu_m);                       // a_m_1_2 represents (nu_a / nu_m)^(1/2)
    R4 = sqrt(nu_c / nu_a) / 3;                        // R4 is a scaling factor based on (nu_c / nu_a)^(1/2) / 3
    R6 = sqrt(nu_c / nu_a) * fast_pow(nu_m / nu_a, (p - 1) / 2) / 3;  // R6 is used in regime 6 for scaling
    R5 = (p - 1) * R6;                                                // R5 scales R6 by a factor of (p - 1)
}

double SynPhotons::spectrum_(double nu) const {
    switch (regime) {
        case 1:
            if (nu <= nu_a) {
                return a_m_1_3 * (nu / nu_a) * (nu / nu_a);
            } else if (nu <= nu_m) {
                return cbrt(nu / nu_m);
            } else if (nu <= nu_c) {
                return fast_pow(nu / nu_m, -(p - 1) / 2);
            } else {
                return c_m_1_2 * fast_pow(nu / nu_m, -p / 2) * fast_exp(-nu / nu_M);
            }
            break;
        case 2:
            if (nu <= nu_m) {
                return m_a_pa4_2 * (nu / nu_m) * (nu / nu_m);
            } else if (nu <= nu_a) {
                return a_m_mpa1_2 * pow52(nu / nu_a);  //  fast_pow(nu / nu_a, 5. / 2);
            } else if (nu <= nu_c) {
                return fast_pow(nu / nu_m, -(p - 1) / 2);
            } else {
                return c_m_1_2 * fast_pow(nu / nu_m, -p / 2) * fast_exp(-nu / nu_M);
            }
            break;
        case 3:
            if (nu <= nu_a) {
                return a_c_1_3 * (nu / nu_a) * (nu / nu_a);
            } else if (nu <= nu_c) {
                return cbrt(nu / nu_c);
            } else if (nu <= nu_m) {
                return sqrt(nu_c / nu);
            } else {
                return c_m_1_2 * fast_pow(nu / nu_m, -p / 2) * fast_exp(-nu / nu_M);
            }
            break;
        case 4:
            if (nu <= nu_a) {
                return (nu / nu_a) * (nu / nu_a);
            } else if (nu <= nu_m) {
                return R4 * sqrt(nu_a / nu);
            } else {
                return R4 * a_m_1_2 * fast_pow(nu / nu_m, -p / 2) * fast_exp(-nu / nu_M);
            }
            break;
        case 5:
            if (nu <= nu_a) {
                return (nu / nu_a) * (nu / nu_a);
            } else {
                return R5 * fast_pow(nu / nu_a, -p / 2) * fast_exp(-nu / nu_M);
            }
            break;
        case 6:
            if (nu <= nu_a) {
                return (nu / nu_a) * (nu / nu_a);
            } else {
                return R6 * fast_pow(nu / nu_a, -p / 2) * fast_exp(-nu / nu_M);
            }
            break;

        default:
            break;
    }
}

// single electron power in the co-moving frame
double syn_p_nu_peak(double B, double p) { return (p - 1) / 2 * sqrt(3) * con::e3 * B / (con::me * con::c2); }

double syn_nu(double gamma, double B) {
    double nu = 3 * con::e * B / (4 * con::pi * con::me * con::c) * gamma * gamma;
    return nu;
}

double syn_gamma(double nu, double B) {
    double gamma = sqrt(nu * 4 * con::pi * con::me * con::c / (3 * con::e * B));
    return gamma;
}

double syn_gamma_M(double B, double zeta, std::vector<Y_IC> const& Ys, double p) {
    double Y0 = Y_IC::Y_T_tilt(Ys);
    double gamma_M = sqrt(6 * con::pi * con::e / (con::sigmaT * B * zeta * (1 + Y0)));
    double Y1 = Y_IC::Y_tilt_at_gamma(Ys, gamma_M, p);

    for (; fabs((Y1 - Y0) / Y0) > 1e-5;) {
        gamma_M = sqrt(6 * con::pi * con::e / (con::sigmaT * B * zeta * (1 + Y1)));
        Y0 = Y1;
        Y1 = Y_IC::Y_tilt_at_gamma(Ys, gamma_M, p);
    }

    return gamma_M;
}

double syn_gamma_m(double e_th, double gamma_M, double eps_e, double n_e, double p) {
    double gamma_bar_minus_1 = eps_e * e_th / (n_e * con::me * con::c2);
    double gamma_m_minus_1 = 1;
    if (p > 2) {
        gamma_m_minus_1 = (p - 2) / (p - 1) * gamma_bar_minus_1;
    } else if (p < 2) {
        // need to check in non-relativistic limit
        gamma_m_minus_1 = pow((2 - p) / (p - 1) * gamma_bar_minus_1 * pow(gamma_M, p - 2), 1 / (p - 1));
    } else {
        gamma_m_minus_1 = root_bisection(
            [=](double x) -> double {
                return (x * log(gamma_M) - (x + 1) * log(x) - gamma_bar_minus_1 - log(gamma_M));
            },
            0, gamma_M);
    }
    return gamma_m_minus_1 + 1;
}

double syn_gamma_c(double t_com, double B, std::vector<Y_IC> const& Ys, double p) {
    // t_com = (6*pi*gamma*me*c^2) /(gamma^2*beta^2*sigma_T*c*B^2*(1 + Y_tilt))
    // double gamma_c = 6 * con::pi * con::me * con::c / (con::sigmaT * B * B * (1 + Y_tilt) * t_com);

    double Y0 = Y_IC::Y_T_tilt(Ys);
    double gamma_bar = 6 * con::pi * con::me * con::c / (con::sigmaT * B * B * (1 + Y0) * t_com);
    double gamma_c = (gamma_bar + sqrt(gamma_bar * gamma_bar + 4)) / 2;
    double Y1 = Y_IC::Y_tilt_at_gamma(Ys, gamma_c, p);

    for (; fabs((Y1 - Y0) / Y0) > 1e-5;) {
        gamma_bar = 6 * con::pi * con::me * con::c / (con::sigmaT * B * B * (1 + Y1) * t_com);
        gamma_c = (gamma_bar + sqrt(gamma_bar * gamma_bar + 4)) / 2;
        Y0 = Y1;
        Y1 = Y_IC::Y_tilt_at_gamma(Ys, gamma_c, p);
    }

    return gamma_c;
}

double syn_gamma_a(double Gamma, double B, double I_syn_peak, double gamma_m, double gamma_c, double gamma_M) {
    double gamma_peak = std::min(gamma_m, gamma_c);
    double nu_peak = syn_nu(gamma_peak, B);
    double ad_idx = adiabatic_index(Gamma);

    double kT = (gamma_peak - 1) * con::me * con::c2 * (ad_idx - 1);
    // 2kT(nv_a/c)^2 = I_peak*(nu_a/nu_peak)^(1/3)
    double nu_a = pow(I_syn_peak * con::c2 / cbrt(nu_peak) / kT / 2, 3. / 5);

    // the nu_peak is not the real peak, peak at nu_a; kT = (gamma_a-1) * me *c^2*(ad_idx-1), I_syn = I_peak;
    if (fabs(gamma_peak - 1) < 1e-6 || nu_a > nu_peak) {
        /*nu_a = pow(I_syn_peak / con::me / 2 / (ad_idx - 1) / sqrt(4 * con::pi / 3 * con::me * con::c / con::e /
         * B),2.0 / 5);*/ //this works only for gamma >> 1
        double nu_M = syn_nu(gamma_M, B);
        double A = sqrt(4 * con::pi / 3 * con::me * con::c / con::e / B);
        double B = I_syn_peak / (2 * con::me * (ad_idx - 1));
        nu_a = root_bisection([=](double x) -> double { return A * x * x * x * x * x - x * x * x * x - B; },
                              sqrt(nu_peak), sqrt(nu_M));
        nu_a *= nu_a;
    }
    double gamma_a = syn_gamma(nu_a, B);
    if (gamma_a < 1) {
        gamma_a = 1;
    }
    return gamma_a;
}

double syn_nu_E_peak(double nu_a, double nu_m, double nu_c) {
    if (order(nu_a, nu_m, nu_c)) {
        return nu_c;
    } else if (order(nu_m, nu_a, nu_c)) {
        return nu_c;
    } else if (order(nu_a, nu_c, nu_m)) {
        return nu_m;
    } else if (order(nu_c, nu_a, nu_m)) {
        if (nu_a * nu_a < nu_m * nu_c) {
            return nu_m;
        } else {
            return nu_a;
        }
    } else if (order(nu_m, nu_c, nu_a)) {
        return nu_a;
    } else if (order(nu_c, nu_m, nu_a)) {
        return nu_a;
    }
}

double syn_nu_E_peak(SynPhotons const& ph) { return syn_nu_E_peak(ph.nu_a, ph.nu_m, ph.nu_c); }

double syn_gamma_N_peak(double gamma_a, double gamma_m, double gamma_c) {
    double gamma_peak = std::min(gamma_m, gamma_c);
    if (gamma_a > gamma_c) {
        return gamma_a;
    } else {
        return gamma_peak;
    }
}

double syn_gamma_N_peak(SynElectrons const& e) { return syn_gamma_N_peak(e.gamma_a, e.gamma_m, e.gamma_c); }

double calc_surface_element(double r, double theta_start, double theta_end) {
    double dcos = std::fabs(cos(theta_end) - cos(theta_start));
    return r * r * dcos;
}

// update electron specific gamma based on new Y parameter
void update_electrons_4_Y(SynElectronsMesh& e, Shock const& shock) {
    th_pool.detach_blocks(0, e.size(), [&](size_t start, size_t end) {
        for (size_t j = start; j < end; ++j) {
            for (size_t k = 0; k < e[0].size(); ++k) {
                double Gamma = shock.Gamma[j][k];
                double t_com = shock.t_com[j][k];
                double B = shock.B[j][k];
                double p = e[j][k].p;
                auto& Ys = e[j][k].Ys;

                auto& electron = e[j][k];

                electron.gamma_M = syn_gamma_M(B, shock.zeta, Ys, p);
                electron.gamma_c = syn_gamma_c(t_com, B, Ys, p);
                electron.gamma_a =
                    syn_gamma_a(Gamma, B, electron.I_nu_peak, electron.gamma_m, electron.gamma_c, electron.gamma_M);
                electron.regime = get_regime(electron.gamma_a, electron.gamma_c, electron.gamma_m);
                electron.gamma_N_peak = syn_gamma_N_peak(electron.gamma_a, electron.gamma_m, electron.gamma_c);
                electron.Y_c = Y_IC::Y_tilt_at_gamma(Ys, electron.gamma_c, p);
            }
        }
    });
    th_pool.wait();
}

SynElectronsMesh gen_syn_electrons(Coord const& coord, Shock const& shock) {
    SynElectronsMesh e = create_syn_electrons_grid(coord.theta.size(), coord.r.size());

    th_pool.detach_blocks(0, coord.theta.size(), [&](size_t start, size_t end) {
        for (size_t j = start; j < end; ++j) {
            for (size_t k = 0; k < coord.r.size(); ++k) {
                constexpr double gamma_syn_limit = 5;
                double dS = calc_surface_element(coord.r[k], coord.theta_b[j], coord.theta_b[j + 1]);
                double Gamma = shock.Gamma[j][k];
                double e_th = shock.e_th[j][k];
                double t_com = shock.t_com[j][k];
                double B = shock.B[j][k];
                double D = shock.width_eff[j][k];
                double n_e = shock.n_p[j][k] * shock.xi;

                auto& electron = e[j][k];

                electron.gamma_M = syn_gamma_M(B, shock.zeta, e[j][k].Ys, shock.p);
                electron.gamma_m = syn_gamma_m(e_th, electron.gamma_M, shock.eps_e, n_e, shock.p);

                // fraction of synchrotron electron; the rest are cyclotron
                double f = std::min(std::pow((gamma_syn_limit - 1) / (electron.gamma_m - 1), 1 - shock.p), 1.0);
                if (f < 1) {
                    electron.gamma_m = gamma_syn_limit;
                }

                electron.n_tot = n_e;
                electron.N_tot = n_e * dS * D * f;

                electron.I_nu_peak = syn_p_nu_peak(B, shock.p) * electron.N_tot / dS / (4 * con::pi);
                electron.gamma_c = syn_gamma_c(t_com, B, e[j][k].Ys, shock.p);
                electron.gamma_a =
                    syn_gamma_a(Gamma, B, electron.I_nu_peak, electron.gamma_m, electron.gamma_c, electron.gamma_M);
                electron.regime = get_regime(electron.gamma_a, electron.gamma_c, electron.gamma_m);
                electron.p = shock.p;
                electron.gamma_N_peak = syn_gamma_N_peak(electron.gamma_a, electron.gamma_m, electron.gamma_c);
            }
        }
    });
    th_pool.wait();
    return e;
}

SynPhotonsMesh gen_syn_photons(SynElectronsMesh const& e, Coord const& coord, Shock const& shock) {
    SynPhotonsMesh ph = create_syn_photons_grid(coord.theta.size(), coord.r.size());

    th_pool.detach_blocks(0, coord.theta.size(), [&](size_t start, size_t end) {
        for (size_t j = start; j < end; ++j) {
            for (size_t k = 0; k < coord.r.size(); ++k) {
                double B = shock.B[j][k];

                ph[j][k].L_nu_peak = syn_p_nu_peak(B, e[j][k].p) * e[j][k].N_tot;
                ph[j][k].E_nu_peak = ph[j][k].L_nu_peak * shock.dt_com[j][k];
                ph[j][k].nu_M = syn_nu(e[j][k].gamma_M, B);
                ph[j][k].nu_m = syn_nu(e[j][k].gamma_m, B);
                ph[j][k].nu_c = syn_nu(e[j][k].gamma_c, B);
                ph[j][k].nu_a = syn_nu(e[j][k].gamma_a, B);
                ph[j][k].regime = e[j][k].regime;
                ph[j][k].nu_E_peak = syn_nu_E_peak(ph[j][k]);
                ph[j][k].p = e[j][k].p;
                ph[j][k].Ys = e[j][k].Ys;
                ph[j][k].Y_c = e[j][k].Y_c;
                ph[j][k].update_constant();
            }
        }
    });
    th_pool.wait();
    return ph;
}

SynPhotonsMesh gen_syn_photons(Coord const& coord, Shock const& shock) {
    SynElectronsMesh e = create_syn_electrons_grid(coord.theta.size(), coord.r.size());
    return gen_syn_photons(e, coord, shock);
}
