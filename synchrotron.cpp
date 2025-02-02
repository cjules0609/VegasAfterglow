
#include "synchrotron.h"

#include <cmath>

#include "afterglow.h"
#include "inverse-compton.h"
#include "macros.h"
#include "physics.h"
#include "utilities.h"

InverseComptonY::InverseComptonY(double nu_m, double nu_c, double B, double Y_T) {
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

InverseComptonY::InverseComptonY(double Y_T) {
    this->Y_T = Y_T;
    regime = 3;
}

InverseComptonY::InverseComptonY() {
    nu_hat_m = 0;
    nu_hat_c = 0;
    gamma_hat_m = 0;
    gamma_hat_c = 0;
    Y_T = 0;
    regime = 0;
}

double InverseComptonY::as_gamma(double gamma, double p) const {
    switch (regime) {
        case 3:
            return Y_T;
            break;
        case 1:
            if (gamma <= gamma_hat_m) {
                return Y_T;
            } else if (gamma <= gamma_hat_c) {
                return Y_T / std::sqrt(gamma / gamma_hat_m);
            } else
                return Y_T * pow43(gamma_hat_c / gamma) / std::sqrt(gamma_hat_c / gamma_hat_m);

            break;
        case 2:
            if (gamma <= gamma_hat_c) {
                return Y_T;
            } else if (gamma <= gamma_hat_m) {
                return Y_T * fastPow(gamma / gamma_hat_c, (p - 3) / 2);
            } else
                return Y_T * pow43(gamma_hat_m / gamma) * fastPow(gamma_hat_m / gamma_hat_c, (p - 3) / 2);

            break;
        default:
            return 0;
            break;
    }
}

double InverseComptonY::as_nu(double nu, double p) const {
    switch (regime) {
        case 3:
            return Y_T;
            break;
        case 1:
            if (nu <= nu_hat_m) {
                return Y_T;
            } else if (nu <= nu_hat_c) {
                return Y_T * std::sqrt(std::sqrt(nu_hat_m / nu));
            } else
                return Y_T * pow23(nu_hat_c / nu) * std::sqrt(std::sqrt(nu_hat_m / nu));

            break;
        case 2:
            if (nu <= nu_hat_c) {
                return Y_T;
            } else if (nu <= nu_hat_m) {
                return Y_T * fastPow(nu / nu_hat_c, (p - 3) / 4);
            } else
                return Y_T * pow23(nu_hat_m / nu) * fastPow(nu_hat_m / nu_hat_c, (p - 3) / 4);

            break;
        default:
            return 0;
            break;
    }
}

double InverseComptonY::Y_Thompson(std::vector<InverseComptonY> const& Ys) {
    double Y_tilt = 0;
    for (auto& Y : Ys) {
        Y_tilt += Y.Y_T;
    }
    return Y_tilt;
}

double InverseComptonY::Y_tilt_gamma(std::vector<InverseComptonY> const& Ys, double gamma, double p) {
    double Y_tilt = 0;
    for (auto& Y : Ys) {
        Y_tilt += Y.as_gamma(gamma, p);
    }
    return Y_tilt;
}

double InverseComptonY::Y_tilt_nu(std::vector<InverseComptonY> const& Ys, double nu, double p) {
    double Y_tilt = 0;
    for (auto& Y : Ys) {
        Y_tilt += Y.as_nu(nu, p);
    }
    return Y_tilt;
}

SynPhotonGrid createSynPhotonGrid(size_t phi_size, size_t theta_size, size_t r_size) {
    SynPhotonGrid grid(boost::extents[phi_size][theta_size][r_size]);
    return grid;
}

SynElectronGrid createSynElectronGrid(size_t phi_size, size_t theta_size, size_t r_size) {
    SynElectronGrid grid(boost::extents[phi_size][theta_size][r_size]);
    return grid;
}

double SynElectrons::columnNumDen(double gamma) const {
    if (gamma < gamma_c) {
        return column_num_den * gammaSpectrum(gamma);
    } else {
        return column_num_den * gammaSpectrum(gamma) * (1 + Y_c) / (1 + InverseComptonY::Y_tilt_gamma(Ys, gamma, p));
    }
}

inline bool order(double a, double b, double c) { return a <= b && b <= c; };

size_t getRegime(double a, double c, double m) {
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
    } else
        return 0;
}

double SynElectrons::gammaSpectrum(double gamma) const {
    switch (regime) {
        case 1:  // same as case 2
        case 2:
            if (gamma <= gamma_m) {
                return 0;
            } else if (gamma <= gamma_c) {
                return (p - 1) * fastPow(gamma / gamma_m, -p) / gamma_m;
            } else
                return (p - 1) * fastPow(gamma / gamma_m, -p) * gamma_c / (gamma * gamma_m) * fastExp(-gamma / gamma_M);

            break;
        case 3:
            if (gamma <= gamma_c) {
                return 0;
            } else if (gamma <= gamma_m) {
                return gamma_c / (gamma * gamma);
            } else
                return gamma_c / (gamma * gamma_m) * fastPow(gamma / gamma_m, -p) * fastExp(-gamma / gamma_M);

            break;
        case 4:  // Gao, Lei, Wu and Zhang 2013 Eq 18
            if (gamma <= gamma_a) {
                return 3 * gamma * gamma / (gamma_a * gamma_a * gamma_a);
            } else if (gamma <= gamma_m) {
                return gamma_c / (gamma * gamma);
            } else
                return gamma_c / (gamma * gamma_m) * fastPow(gamma / gamma_m, -p) * fastExp(-gamma / gamma_M);

            break;
        case 5:  // Gao, Lei, Wu and Zhang 2013 Eq 19
            if (gamma <= gamma_a) {
                return 3 * gamma * gamma / (gamma_a * gamma_a * gamma_a);
            } else
                return (p - 1) * gamma_c / (gamma * gamma_m) * fastPow(gamma / gamma_m, -p) * fastExp(-gamma / gamma_M);

            break;
        case 6:  // Gao, Lei, Wu and Zhang 2013 Eq 20
            if (gamma <= gamma_a) {
                return 3 * gamma * gamma / (gamma_a * gamma_a * gamma_a);
            } else
                return fastPow(gamma_m, p - 1) * gamma_c * fastPow(gamma, -(p + 1)) * fastExp(-gamma / gamma_M);

            break;
        default:
            return 0;
    }
}

double SynPhotons::I_nu(double nu) const {
    if (nu < nu_c) {
        return e->I_nu_peak * spectrum(nu);
    } else {
        return e->I_nu_peak * spectrum(nu) * (1 + e->Y_c) / (1 + InverseComptonY::Y_tilt_nu(e->Ys, nu, e->p));
    }
}

void SynPhotons::updateConstant() {
    // Update constants based on spectral parameters
    double p = e->p;
    a_m_1_3 = std::cbrt(nu_a / nu_m);                 // a_m_1_3 represents (nu_a / nu_m)^(1/3)
    c_m_1_2 = std::sqrt(nu_c / nu_m);                 // c_m_1_2 represents (nu_c / nu_m)^(1/2)
    m_a_pa4_2 = fastPow(nu_m / nu_a, (p + 4) / 2);    // m_a_pa4_2 represents (nu_m / nu_a)^((p+4)/2)
    a_m_mpa1_2 = fastPow(nu_a / nu_m, (-p + 1) / 2);  // a_m_mpa1_2 represents (nu_a / nu_m)^((-p+1)/2)
    a_c_1_3 = std::cbrt(nu_a / nu_c);                 // a_c_1_3 represents (nu_a / nu_c)^(1/3)
    a_m_1_2 = std::sqrt(nu_a / nu_m);                 // a_m_1_2 represents (nu_a / nu_m)^(1/2)
    R4 = std::sqrt(nu_c / nu_a) / 3;                  // R4 is a scaling factor based on (nu_c / nu_a)^(1/2) / 3
    R6 = std::sqrt(nu_c / nu_a) * fastPow(nu_m / nu_a, (p - 1) / 2) / 3;  // R6 is used in regime 6 for scaling
    // R5 = (p - 1) * R6;                                                // R5 scales R6 by a factor of (p - 1)
}

double SynPhotons::spectrum(double nu) const {
    double p = e->p;
    switch (e->regime) {
        case 1:
            if (nu <= nu_a) {
                return a_m_1_3 * (nu / nu_a) * (nu / nu_a);
            }
            if (nu <= nu_m) {
                return std::cbrt(nu / nu_m);
            }
            if (nu <= nu_c) {
                return fastPow(nu / nu_m, -(p - 1) / 2);
            }
            return c_m_1_2 * fastPow(nu / nu_m, -p / 2) * fastExp(-nu / nu_M);

            break;
        case 2:
            if (nu <= nu_m) {
                return m_a_pa4_2 * (nu / nu_m) * (nu / nu_m);
            }
            if (nu <= nu_a) {
                return a_m_mpa1_2 * pow52(nu / nu_a);  //  fast_pow(nu / nu_a, 5. / 2);
            }
            if (nu <= nu_c) {
                return fastPow(nu / nu_m, -(p - 1) / 2);
            }
            return c_m_1_2 * fastPow(nu / nu_m, -p / 2) * fastExp(-nu / nu_M);

            break;
        case 3:
            if (nu <= nu_a) {
                return a_c_1_3 * (nu / nu_a) * (nu / nu_a);
            }
            if (nu <= nu_c) {
                return std::cbrt(nu / nu_c);
            }
            if (nu <= nu_m) {
                return std::sqrt(nu_c / nu);
            }
            return c_m_1_2 * fastPow(nu / nu_m, -p / 2) * fastExp(-nu / nu_M);

            break;
        case 4:
            if (nu <= nu_a) {
                return (nu / nu_a) * (nu / nu_a);
            }
            if (nu <= nu_m) {
                return R4 * std::sqrt(nu_a / nu);
            }
            return R4 * a_m_1_2 * fastPow(nu / nu_m, -p / 2) * fastExp(-nu / nu_M);

            break;
        case 5:
            if (nu <= nu_a) {
                return (nu / nu_a) * (nu / nu_a);
            }
            return (p - 1) * R6 * fastPow(nu / nu_a, -p / 2) * fastExp(-nu / nu_M);

            break;
        case 6:
            if (nu <= nu_a) {
                return (nu / nu_a) * (nu / nu_a);
            }
            return R6 * fastPow(nu / nu_a, -p / 2) * fastExp(-nu / nu_M);

            break;

        default:
            break;
    }
}

// single electron power in the co-moving frame
double syn_p_nu_peak(double B, double p) { return (p - 1) / 2 * std::sqrt(3) * con::e3 * B / (con::me * con::c2); }

double syn_nu(double gamma, double B) {
    double nu = 3 * con::e * B / (4 * con::pi * con::me * con::c) * gamma * gamma;
    return nu;
}

double syn_gamma(double nu, double B) {
    double gamma = std::sqrt(nu * 4 * con::pi * con::me * con::c / (3 * con::e * B));
    return gamma;
}

double syn_gamma_M(double B, std::vector<InverseComptonY> const& Ys, double p) {
    if (B == 0) {
        return std::numeric_limits<double>::infinity();
    }
    double Y0 = InverseComptonY::Y_Thompson(Ys);
    double gamma_M = std::sqrt(6 * con::pi * con::e / (con::sigmaT * B * (1 + Y0)));
    double Y1 = InverseComptonY::Y_tilt_gamma(Ys, gamma_M, p);

    for (; std::fabs((Y1 - Y0) / Y0) > 1e-5;) {
        gamma_M = std::sqrt(6 * con::pi * con::e / (con::sigmaT * B * (1 + Y1)));
        Y0 = Y1;
        Y1 = InverseComptonY::Y_tilt_gamma(Ys, gamma_M, p);
    }

    return gamma_M;
}

double syn_gamma_m(double Gamma_rel, double gamma_M, double eps_e, double p, double xi) {
    double gamma_bar_minus_1 = eps_e * (Gamma_rel - 1) * con::mp / (xi * con::me);
    double gamma_m_minus_1 = 1;
    if (p > 2) {
        gamma_m_minus_1 = (p - 2) / (p - 1) * gamma_bar_minus_1;
    } else if (p < 2) {
        // need to check in non-relativistic limit
        gamma_m_minus_1 = std::pow((2 - p) / (p - 1) * gamma_bar_minus_1 * std::pow(gamma_M, p - 2), 1 / (p - 1));
    } else {
        gamma_m_minus_1 = rootBisection(
            [=](double x) -> double {
                return (x * std::log(gamma_M) - (x + 1) * std::log(x) - gamma_bar_minus_1 - std::log(gamma_M));
            },
            0, gamma_M);
    }
    return gamma_m_minus_1 + 1;
}

double syn_gamma_c(double t_com, double B, std::vector<InverseComptonY> const& Ys, double p) {
    // t_com = (6*pi*gamma*me*c^2) /(gamma^2*beta^2*sigma_T*c*B^2*(1 + Y_tilt))
    // double gamma_c = 6 * con::pi * con::me * con::c / (con::sigmaT * B * B * (1 + Y_tilt) * t_com);

    double Y0 = InverseComptonY::Y_Thompson(Ys);
    double gamma_bar = 6 * con::pi * con::me * con::c / (con::sigmaT * B * B * (1 + Y0) * t_com);
    double gamma_c = (gamma_bar + std::sqrt(gamma_bar * gamma_bar + 4)) / 2;
    double Y1 = InverseComptonY::Y_tilt_gamma(Ys, gamma_c, p);

    for (; std::fabs((Y1 - Y0) / Y0) > 1e-5;) {
        gamma_bar = 6 * con::pi * con::me * con::c / (con::sigmaT * B * B * (1 + Y1) * t_com);
        gamma_c = (gamma_bar + std::sqrt(gamma_bar * gamma_bar + 4)) / 2;
        Y0 = Y1;
        Y1 = InverseComptonY::Y_tilt_gamma(Ys, gamma_c, p);
    }

    return gamma_c;
}

double syn_gamma_a(double Gamma_rel, double B, double I_syn_peak, double gamma_m, double gamma_c, double gamma_M) {
    double gamma_peak = std::min(gamma_m, gamma_c);
    double nu_peak = syn_nu(gamma_peak, B);
    double ad_idx = adiabaticIndex(Gamma_rel);

    double kT = (gamma_peak - 1) * con::me * con::c2 * (ad_idx - 1);
    // 2kT(nv_a/c)^2 = I_peak*(nu_a/nu_peak)^(1/3)
    double nu_a = std::pow(I_syn_peak * con::c2 / std::cbrt(nu_peak) / kT / 2, 3. / 5);

    // the nu_peak is not the real peak, peak at nu_a; kT = (gamma_a-1) * me *c^2*(ad_idx-1), I_syn = I_peak;
    if (fabs(gamma_peak - 1) < 1e-6 || nu_a > nu_peak) {
        /*nu_a = pow(I_syn_peak / con::me / 2 / (ad_idx - 1) / sqrt(4 * con::pi / 3 * con::me * con::c / con::e /
         * B),2.0 / 5);*/ //this works only for gamma >> 1
        double nu_M = syn_nu(gamma_M, B);
        double A = std::sqrt(4 * con::pi / 3 * con::me * con::c / con::e / B);
        double B = I_syn_peak / (2 * con::me * (ad_idx - 1));
        nu_a = rootBisection([=](double x) -> double { return A * x * x * x * x * x - x * x * x * x - B; },
                             std::sqrt(nu_peak), std::sqrt(nu_M));
        nu_a *= nu_a;
    }
    double gamma_a = syn_gamma(nu_a, B);
    if (gamma_a < 1) {
        gamma_a = 1;
    }
    return gamma_a;
}
/*
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
}*/

// double syn_nu_E_peak(SynPhotons const& ph) { return syn_nu_E_peak(ph.nu_a, ph.nu_m, ph.nu_c); }

double syn_gamma_N_peak(double gamma_a, double gamma_m, double gamma_c) {
    double gamma_peak = std::min(gamma_m, gamma_c);
    if (gamma_a > gamma_c) {
        return gamma_a;
    } else {
        return gamma_peak;
    }
}

double syn_gamma_N_peak(SynElectrons const& e) { return syn_gamma_N_peak(e.gamma_a, e.gamma_m, e.gamma_c); }

// update electron specific gamma based on new Y parameter
void updateElectrons4Y(SynElectronGrid& e, Shock const& shock) {
    auto [phi_size, theta_size, r_size] = shock.shape();

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < r_size; ++k) {
                double Gamma_rel = shock.Gamma_rel[i][j][k];
                double t_com = shock.t_com[i][j][k];
                double B = shock.B[i][j][k];
                double p = e[i][j][k].p;
                auto& Ys = e[i][j][k].Ys;

                auto& electron = e[i][j][k];

                electron.gamma_M = syn_gamma_M(B, Ys, p);
                electron.gamma_c = syn_gamma_c(t_com, B, Ys, p);
                electron.gamma_a =
                    syn_gamma_a(Gamma_rel, B, electron.I_nu_peak, electron.gamma_m, electron.gamma_c, electron.gamma_M);
                electron.regime = getRegime(electron.gamma_a, electron.gamma_c, electron.gamma_m);
                electron.gamma_N_peak = syn_gamma_N_peak(electron.gamma_a, electron.gamma_m, electron.gamma_c);
                electron.Y_c = InverseComptonY::Y_tilt_gamma(Ys, electron.gamma_c, p);
            }
        }
    }
}

SynElectronGrid genSynElectrons(Shock const& shock, double p, double xi) {
    auto [phi_size, theta_size, r_size] = shock.shape();

    SynElectronGrid electrons = createSynElectronGrid(phi_size, theta_size, r_size);

    constexpr double gamma_syn_limit = 3;

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < r_size; ++k) {
                double Gamma_rel = shock.Gamma_rel[i][j][k];
                double t_com = shock.t_com[i][j][k];
                double B = shock.B[i][j][k];
                double Sigma = shock.column_num_den[i][j][k];

                auto& e = electrons[i][j][k];

                e.gamma_M = syn_gamma_M(B, electrons[i][j][k].Ys, p);
                e.gamma_m = syn_gamma_m(Gamma_rel, e.gamma_M, shock.eps_e, p, xi);
                // fraction of synchrotron electron; the rest are cyclotron
                double f = 1.;
                if (1 < e.gamma_m && e.gamma_m < gamma_syn_limit) {
                    f = std::min(std::pow((gamma_syn_limit - 1) / (e.gamma_m - 1), 1 - p), 1.0);
                    e.gamma_m = gamma_syn_limit;
                }
                e.column_num_den = Sigma * f;
                e.I_nu_peak = syn_p_nu_peak(B, p) * e.column_num_den / (4 * con::pi);
                e.gamma_c = syn_gamma_c(t_com, B, electrons[i][j][k].Ys, p);
                e.gamma_a = syn_gamma_a(Gamma_rel, B, e.I_nu_peak, e.gamma_m, e.gamma_c, e.gamma_M);
                e.regime = getRegime(e.gamma_a, e.gamma_c, e.gamma_m);
                e.gamma_N_peak = syn_gamma_N_peak(e.gamma_a, e.gamma_m, e.gamma_c);
                e.p = p;
            }
        }
    }
    return electrons;
}

SynPhotonGrid genSynPhotons(Shock const& shock, SynElectronGrid const& e) {
    auto [phi_size, theta_size, r_size] = shock.shape();

    SynPhotonGrid ph = createSynPhotonGrid(phi_size, theta_size, r_size);

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < r_size; ++k) {
                double B = shock.B[i][j][k];

                ph[i][j][k].nu_M = syn_nu(e[i][j][k].gamma_M, B);
                ph[i][j][k].nu_m = syn_nu(e[i][j][k].gamma_m, B);
                ph[i][j][k].nu_c = syn_nu(e[i][j][k].gamma_c, B);
                ph[i][j][k].nu_a = syn_nu(e[i][j][k].gamma_a, B);
                ph[i][j][k].e = &(e[i][j][k]);
                ph[i][j][k].updateConstant();
            }
        }
    }
    return ph;
}
/*
SynPhotonGrid genSynPhotons(Shock const& shock, double p, double xi) {
    size_t phi_size = shock.Gamma_rel.size();
    size_t theta_size = shock.Gamma_rel[0].size();
    size_t r_size = shock.Gamma_rel[0][0].size();
    SynElectronGrid e = createSynElectronGrid(phi_size, theta_size, r_size);
    return genSynPhotons(e, shock);
}*/
