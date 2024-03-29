#include <algorithm>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

using MeshGrid = std::vector<std::vector<double>>;
using MeshGrid3d = std::vector<std::vector<std::vector<double>>>;
using Array = std::vector<double>;
using Profile = std::function<double(double)>;
using Profile2d = std::function<double(double, double)>;

static Profile2d noInjection = Profile2d([](double theta, double t_lab) { return 0; });

template <typename Fun>
auto root_bisection(Fun f, decltype(f(0)) low, decltype(f(0)) high) -> decltype(f(0)) {
    using Scalar = decltype(f(0));

    for (; fabs((high - low)) > fabs(high) * 1e-6;) {
        Scalar mid = 0.5 * (high + low);
        if (f(mid) > 0)
            high = mid;
        else
            low = mid;
    }
    return 0.5 * (high + low);
}

namespace consts {
    constexpr double c = 1;  //
    constexpr double mp = 1.67e-24 / 2e33;
    constexpr double me = mp / 1836;
    // constexpr double kB = 1.38e-16;
    constexpr double e = 4.8e-10 / 4.472136e16 / 5.809475e19 * 500;  // M^(1/2)L^(3/2)/T
    constexpr double pi = 3.14159265358979323846;
    constexpr double sigmaT = 6.65e-25 / 1.5e13 / 1.5e13;
}  // namespace consts

MeshGrid createGrid(Array const& r, Array const& theta, double val = 0) {
    return MeshGrid(theta.size(), Array(r.size(), val));
}

MeshGrid3d createGrid3d(Array const& r, Array const& theta, Array const& phi, double val = 0) {
    return MeshGrid3d(phi.size(), MeshGrid(theta.size(), Array(r.size(), val)));
}

Array logspace(double start, double end, size_t num) {
    Array result(num);
    double log_start = std::log10(start);
    double log_end = std::log10(end);
    double step = (log_end - log_start) / (num - 1);
    for (size_t i = 0; i < num; i++) {
        result[i] = std::pow(10, log_start + i * step);
    }
    return result;
}

Array linspace(double start, double end, size_t num) {
    Array result(num);
    double step = (end - start) / (num - 1);
    for (size_t i = 0; i < num; i++) {
        result[i] = start + i * step;
    }
    return result;
}

Array zeros(size_t num) {
    Array result(num);
    std::fill(result.begin(), result.end(), 0);
    return result;
}

Array ones(size_t num) {
    Array result(num);
    std::fill(result.begin(), result.end(), 1);
    return result;
}

class BlastWave {
   public:
    Profile2d dEdOmega;
    Profile Gamma0;
};

class Medium {
   public:
    Medium(Profile rho, Profile mass, double eps_rad, double eps_e, double eps_B, double Y_tilt = 0, double xi = 1)
        : rho(rho), mass(mass), eps_rad(eps_rad), eps_e(eps_e), eps_B(eps_B), Y_tilt{Y_tilt}, xi(xi){};
    double eps_rad{1};
    double eps_e{0.3};
    double eps_B{0.01};
    double Y_tilt{0};
    double zeta{1};
    double xi{1};
    Profile rho;
    Profile mass;
};

class DynamicsEqn {
   public:
    DynamicsEqn(Medium medium, BlastWave blast, double theta_lo, double theta_hi)
        : medium(medium),
          blast(blast),
          theta_lo(theta_lo),
          theta_hi(theta_hi),
          theta(0.5 * (theta_lo + theta_hi)),
          dOmega(4 * consts::pi * std::fabs(std::cos(theta_hi) - std::cos(theta_lo))){};  // bipolar outflow

    void operator()(Array const& y, Array& dydr, double r) {
        double Gamma = y[0];
        double u = y[1];
        double t_obs = y[2];
        double t_lab = y[3];
        dydr[0] = dGammadr(r, Gamma, u, t_lab);
        dydr[1] = dUdr(r, Gamma, u, t_lab);
        dydr[2] = dt_labdr(Gamma);
    };

   private:
    inline double dGammadr(double r, double Gamma, double u, double t_lab) {
        double gamma_eos = (4 * Gamma + 1) / (3 * Gamma);
        double dm = medium.mass(r) * dOmega / (4 * consts::pi);
        double dM0 = blast.dEdOmega(theta, t_lab) * dOmega / (blast.Gamma0(theta) * consts::c * consts::c);
        double Gamma2 = Gamma * Gamma;
        double a1 = dOmega * r * r * medium.rho(r) / Gamma * (Gamma2 - 1) * (gamma_eos * Gamma - gamma_eos + 1);
        double a2 = -(gamma_eos - 1) / Gamma * (gamma_eos * Gamma2 - gamma_eos + 1) * 3 * u / r;
        double b1 = (dM0 + dm) * consts::c * consts::c;
        double b2 = (gamma_eos * gamma_eos * (Gamma2 - 1) + 3 * gamma_eos - 2) * u / Gamma2;
        return -(a1 + a2) / (b1 + b2);
    };

    inline double dUdr(double r, double Gamma, double u, double t_lab) {
        double gamma_eos = (4 * Gamma + 1) / (3 * Gamma);
        double E = dOmega * r * r * medium.rho(r) * consts::c * consts::c;
        return (1 - medium.eps_e * medium.eps_rad) * (Gamma - 1) * E -
               (gamma_eos - 1) * (3 / r - dGammadr(r, Gamma, u, t_lab) / Gamma) * u;
    };

    inline double dt_labdr(double Gamma) {
        double sqr = std::sqrt(Gamma * Gamma - 1);
        return (Gamma - sqr) / (sqr * consts::c);
    };

    Medium medium;
    BlastWave blast;
    double theta_lo;
    double theta_hi;
    double theta;
    double dOmega;
};

template <typename Eqn>
void solveShell(Array const& r, Array& Gamma, Array& u, Eqn const& eqn) {
    using namespace boost::numeric::odeint;
    double atol = 0;
    double rtol = 1e-9;
    auto stepper = bulirsch_stoer_dense_out<std::vector<double>>{atol, rtol};
    Array state{0, 0, 0};
    double dr = (r[1] - r[0]) / 1000;
    int i = 0;
    stepper.initialize(Array{Gamma[0], u[0], 0.0}, r[0], dr);
    for (; stepper.current_time() <= r.back();) {
        auto [t_last, t_current] = stepper.do_step(eqn);
        for (; stepper.current_time() > r[i + 1] && i + 1 < r.size();) {
            i++;
            stepper.calc_state(r[i], state);
            Gamma[i] = state[0];
            u[i] = state[1];
        }
    }
}

Array boundary2center(Array const& boundary) {
    Array center(boundary.size() - 1);
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = 0.5 * (boundary[i] + boundary[i + 1]);
    }
    return center;
}

Array boundary2centerlog(Array const& boundary) {
    Array center(boundary.size() - 1);
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = std::sqrt(boundary[i] * boundary[i + 1]);
    }
    return center;
}

auto solveDynamics(BlastWave const& blast, Medium const& medium, Array const& r, Array const& theta) {
    Array r_c = boundary2centerlog(r);
    Array theta_c = boundary2center(theta);
    MeshGrid Gamma = createGrid(r_c, theta_c, 1);
    MeshGrid U = createGrid(r_c, theta_c, 0);

    for (size_t i = 0; i < theta_c.size(); ++i) {
        double dcos = std::fabs(std::cos(theta[i + 1]) - std::cos(theta[i]));
        double dphi = 2 * consts::pi;
        double dOmega = 2 * dphi * dcos;  // bipolar outflow
        Gamma[i][0] = blast.Gamma0(theta_c[i]);
        U[i][0] = medium.mass(r_c[0]) * dOmega / (4 * consts::pi) * consts::c * consts::c * (Gamma[i][0] - 1);
        auto eqn = DynamicsEqn(medium, blast, theta[i], theta[i + 1]);
        solveShell(r_c, Gamma[i], U[i], eqn);
    }
    return std::make_tuple(Gamma, U);
}

void printArray(Array const& arr) {
    for (auto const& a : arr) {
        std::cout << a << " ";
    }
    std::cout << std::endl;
}

void write2file(MeshGrid3d const& array, std::string const& filename) {
    std::ofstream file(filename);
    for (size_t i = 0; i < array.size(); ++i) {
        for (size_t j = 0; j < array[i].size(); ++j) {
            for (size_t k = 0; k < array[i][j].size(); ++k) {
                file << array[i][j][k] << " ";
            }
            file << std::endl;
        }
    }
}
void write2file(MeshGrid const& grid, std::string const& filename) {
    std::ofstream file(filename);
    for (size_t i = 0; i < grid.size(); ++i) {
        for (size_t j = 0; j < grid[i].size(); ++j) {
            file << grid[i][j] << " ";
        }
        file << std::endl;
    }
}

void write2file(Array const& array, std::string const& filename) {
    std::ofstream file(filename);
    for (size_t i = 0; i < array.size(); ++i) {
        file << array[i] << " ";
    }
}

inline double stepfun(double x) { return x > 0 ? 1 : 0; }

Medium create_ISM(double n_ism, double eps_rad, double eps_e, double eps_B = 0.01) {
    auto rho = [=](double r) { return n_ism * consts::mp; };
    auto m = [=](double r) { return 4 * consts::pi / 3 * r * r * r * n_ism * consts::mp; };
    return Medium{rho, m, eps_rad, eps_e, eps_B};
}

double E_iso2Gamma0(double E_iso, double gamma_max, double E) {
    double u = pow(E / E_iso, 1.0 / 4) * gamma_max;
    double gamma = sqrt(1 + u * u);
    return gamma;
}
BlastWave createIsotropicJet(double E_iso, double Gamma0, Profile2d inject = noInjection) {
    auto dEdOmega = [=](double theta, double t_lab) { return E_iso / (4 * consts::pi) + inject(theta, t_lab); };
    auto Gamma = [=](double theta) { return Gamma0; };
    return BlastWave{dEdOmega, Gamma};
}

BlastWave createTopHatJet(double E_iso, double Gamma0, double theta_c, Profile2d inject = noInjection) {
    auto dEdOmega = [=](double theta, double t_lab) {
        return (theta < theta_c ? (E_iso / (4 * consts::pi)) : 0) + inject(theta, t_lab);
    };
    auto Gamma = [=](double theta) { return theta < theta_c ? Gamma0 : 1 + 1e-6; };
    return BlastWave{dEdOmega, Gamma};
}

BlastWave createPowerLawJet(double E_iso_on_axis, double Gamma0_on_axis, double theta_m, double k,
                            Profile2d inject = noInjection) {
    double e0 = E_iso_on_axis / (4 * consts::pi);
    auto dEdOmega = [=](double theta, double t_lab) {
        return (theta < theta_m ? e0 : e0 * std::pow(theta / theta_m, -k)) + inject(theta, t_lab);
    };
    auto Gamma = [=](double theta) { return E_iso2Gamma0(e0, Gamma0_on_axis, dEdOmega(theta, 0)); };
    return BlastWave{dEdOmega, Gamma};
}

BlastWave createGaussianJet(double E_iso_on_axis, double Gamma0_on_axis, double theta_c,
                            Profile2d inject = noInjection) {
    double e0 = E_iso_on_axis / (4 * consts::pi);
    auto dEdOmega = [=](double theta, double t_lab) {
        return e0 * std::exp(-theta * theta / (2 * theta_c * theta_c)) + inject(theta, t_lab);
    };
    auto Gamma = [=](double theta) { return E_iso2Gamma0(e0, Gamma0_on_axis, dEdOmega(theta, 0)); };
    return BlastWave{dEdOmega, Gamma};
}

Profile2d createIsoPowerLawInjection(double L0, double t0, double t_wait, double q) {
    double t0q = std::pow(t0, q);
    double tw1q = std::pow(t0, 1 - q);
    double Omega = 4 * consts::pi;
    if (std::fabs(q - 1) > 1e-6) {
        return [=](double theta, double t_lab) {
            return stepfun(t_lab - t_wait) * L0 * t0q / (1 - q) * (std::pow(t_lab + t0 - t_wait, 1 - q) - tw1q) / Omega;
        };
    } else {
        return [=](double theta, double t_lab) {
            return stepfun(t_lab - t_wait) * L0 * t0 * std::log((t_lab + t0 - t_wait) / t_wait) / Omega;
        };
    }
}

double interp(double x0, Array const& x, Array const& y) {
    if (x0 < x[0]) {
        return y[0];
    } else if (x0 > x.back()) {
        return y.back();
    } else {
        auto it = std::lower_bound(x.begin(), x.end(), x0);
        size_t i = it - x.begin();
        if (x0 == x[i]) {
            return y[i];
        } else {
            double a = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
            return y[i - 1] + a * (x0 - x[i - 1]);
        }
    }
}

class TobsEqn {
   public:
    void operator()(double const& t, double& dtdr, double r) {
        double Gamma = interp(r, r_, Gamma_);
        double sqr = std::sqrt(Gamma * Gamma - 1);
        dtdr = (Gamma - sqr * cos_) / (sqr * consts::c);
    };

    Array Gamma_;
    Array r_;
    double cos_;
};

MeshGrid3d calc_t_obs(MeshGrid const& Gamma, Array const& r, Array const& theta, Array const& phi, double theta_obs) {
    auto r_c = boundary2centerlog(r);
    auto phi_c = boundary2center(phi);
    auto theta_c = boundary2center(theta);

    MeshGrid3d t_obs = createGrid3d(r_c, theta_c, phi_c);
    using namespace boost::numeric::odeint;
    double atol = 0;
    double rtol = 1e-9;
    auto stepper = bulirsch_stoer_dense_out<double>{atol, rtol};
    TobsEqn eqn;
    eqn.r_ = r_c;
    double dr0 = (r[1] - r[0]) / 1000;

    for (size_t j = 0; j < theta_c.size(); ++j) {
        double theta_ = theta_c[j];
        eqn.Gamma_ = Gamma[j];
        double Gamma0 = Gamma[j][0];
        double beta0 = sqrt(1 - 1 / Gamma0 / Gamma0);
        for (size_t i = 0; i < phi_c.size(); ++i) {
            double phi_ = phi_c[i];
            eqn.cos_ = sin(theta_) * cos(phi_) * sin(theta_obs) + cos(theta_) * cos(theta_obs);
            int k = 0;
            double t_obs0 = r_c[0] * (1 - beta0 * eqn.cos_) / consts::c / beta0;
            stepper.initialize(t_obs0, r_c[0], dr0);
            for (; stepper.current_time() <= r_c.back();) {
                auto [t_last, t_current] = stepper.do_step(eqn);
                for (; stepper.current_time() > r_c[k + 1] && k + 1 < r_c.size();) {
                    k++;
                    stepper.calc_state(r_c[k], t_obs[i][j][k]);
                }
            }
        }
    }
    return t_obs;
}

MeshGrid3d calc_doppler(MeshGrid const& Gamma, Array const& r, Array const& theta, Array const& phi, double theta_obs) {
    auto r_c = boundary2centerlog(r);
    auto phi_c = boundary2center(phi);
    auto theta_c = boundary2center(theta);
    MeshGrid3d doppler = createGrid3d(r_c, theta_c, phi_c);
    for (size_t i = 0; i < phi_c.size(); ++i) {
        double phi_ = phi_c[i];
        for (size_t j = 0; j < theta_c.size(); ++j) {
            double theta_ = theta_c[j];
            double cos_ = sin(theta_) * cos(phi_) * sin(theta_obs) + cos(theta_) * cos(theta_obs);
            for (size_t k = 0; k < Gamma[j].size(); ++k) {
                double gamma_ = Gamma[j][k];
                double beta = sqrt(1 - 1 / gamma_ / gamma_);
                doppler[i][j][k] = 1 / (gamma_ * (1 - beta * cos_));
            }
        }
    }
    return doppler;
}

inline double synchrotron_nu(double gamma, double B) {
    double nu = 3 * consts::e * B / 4 / consts::pi / consts::me / consts::c * gamma * gamma;
    return nu;
}

inline double synchrotron_specific_emissivity(double pel, double nu, double nu_m, double nu_c, double nu_a) {
    if (nu_a < nu_m && nu_m < nu_c) {
        if (nu <= nu_a) {
            return pow(nu_a / nu_m, 1.0 / 3) * pow(nu / nu_a, 2);
        } else if (nu <= nu_m) {
            return pow(nu / nu_m, 1.0 / 3);
        } else if (nu <= nu_c) {
            return pow(nu / nu_m, -(pel - 1) / 2);
        } else {
            return pow(nu_c / nu_m, -(pel - 1) / 2) * pow(nu / nu_c, -pel / 2);
        }
    } else if (nu_m < nu_a && nu_a < nu_c) {
        if (nu <= nu_m) {
            return pow(nu_m / nu_a, (pel + 4) / 2) * pow(nu / nu_m, 2);
        } else if (nu <= nu_a) {
            return pow(nu_a / nu_m, -(pel - 1) / 2) * pow(nu / nu_a, 5.0 / 2);
        } else if (nu <= nu_c) {
            return pow(nu / nu_m, -(pel - 1) / 2);
        } else {
            return pow(nu_c / nu_m, -(pel - 1) / 2) * pow(nu / nu_c, -pel / 2);
        }
    } else if (nu_a < nu_c && nu_c < nu_m) {
        if (nu <= nu_a) {
            return pow(nu_a / nu_c, 1.0 / 3) * pow(nu / nu_a, 2);
        } else if (nu <= nu_c) {
            return pow(nu / nu_c, 1.0 / 3);
        } else if (nu <= nu_m) {
            return pow(nu / nu_c, -1 / 2);
        } else {
            return pow(nu_m / nu_c, -1 / 2) * pow(nu / nu_m, -pel / 2);
        }
    } else if (nu_c < nu_a && nu_a < nu_m) {
        if (nu <= nu_a) {
            return pow(nu / nu_a, 2);
        } else if (nu <= nu_m) {
            return pow(nu / nu_a, -1.0 / 2) * pow(nu_c / nu_a, 1.0 / 2) / 3;
        } else {
            return pow(nu_m / nu_a, -1.0 / 2) * pow(nu / nu_m, -pel / 2) * pow(nu_c / nu_a, 1.0 / 2) / 3;
        }
    } else if (nu_m < nu_c && nu_c < nu_a) {
        if (nu <= nu_a) {
            return pow(nu / nu_a, 2);
        } else {
            return (pel - 1) / 3 * pow(nu / nu_a, -pel / 2) * pow(nu_m / nu_a, (pel - 1) / 2) *
                   pow(nu_c / nu_a, 1.0 / 2);
        }
    } else if (nu_c < nu_m && nu_m < nu_a) {
        if (nu <= nu_a) {
            return pow(nu / nu_a, 2);
        } else {
            return pow(nu / nu_a, -pel / 2) * pow(nu_m / nu_a, (pel - 1) / 2) * pow(nu_c / nu_a, 1.0 / 2) / 3;
        }
    }
}
auto calc_Synchrotron_Radiation(Array const& r, Array const& theta, Array const& phi, MeshGrid const& Gamma,
                                MeshGrid3d const& t_obs, MeshGrid3d const& doppler, Medium const& medium, double pel,
                                double nu_obs) {
    if (pel < 1) {
        std::cout << "pel should be larger than 1" << std::endl;
        exit(1);
    }
    auto r_c = boundary2centerlog(r);
    auto phi_c = boundary2center(phi);
    auto theta_c = boundary2center(theta);
    MeshGrid3d dFnu = createGrid3d(r_c, theta_c, phi_c);
    MeshGrid nv_m = createGrid(r_c, theta_c);
    MeshGrid nv_c = createGrid(r_c, theta_c);
    MeshGrid nv_a = createGrid(r_c, theta_c);
    MeshGrid nv_M = createGrid(r_c, theta_c);
    double dphi = phi[1] - phi[0];
    double dcos = std::fabs(std::cos(theta[1]) - std::cos(theta[0]));
    double dOmega = dphi * dcos;

    for (size_t i = 0; i < phi_c.size(); ++i) {
        double phi_ = phi_c[i];
        for (size_t j = 0; j < theta_c.size(); ++j) {
            double theta_ = theta_c[j];
            for (size_t k = 0; k < r.size(); ++k) {
                double dr = r[k + 1] - r[k];
                double r_ = r_c[k];
                double Gamma_ = Gamma[j][k];
                double doppler_ = doppler[i][j][k];
                double t_obs_ = t_obs[i][j][k];
                double nu_prime = nu_obs / doppler_;

                double Bprime_ =
                    sqrt(8 * consts::pi * medium.eps_B * medium.rho(r_) * (4 * Gamma_ * Gamma_ - 4 * Gamma_)) *
                    consts::c;

                double gamma_M_ =
                    sqrt(6 * consts::pi * consts::e / (consts::sigmaT * Bprime_ * medium.zeta * (1 + medium.Y_tilt)));

                double gamma_m_ = 1;
                if (pel > 2) {
                    gamma_m_ = (pel - 2) / (pel - 1) * medium.eps_e * (Gamma_ - 1) * 1836 / medium.xi;
                } else if (pel > 1) {
                    gamma_m_ = root_bisection(
                        [=](double x) -> double {
                            return (x * log(gamma_M_) - x * log(x) - medium.eps_e * (Gamma_ - 1) * 1836 / medium.xi);
                        },
                        1, gamma_M_);
                }
                /*double gamma_c_ = 6 * consts::pi * consts::me * consts::c /
                                  (consts::sigmaT * doppler_ * t_obs_ * Bprime_ * Bprime_ * (1 + medium.Y_tilt));*/

                double gamma_c_ = 15 * consts::pi * consts::me * consts::c * consts::c * sqrt(Gamma_ * Gamma_ - 1) /
                                      (consts::sigmaT * r_ * Bprime_ * Bprime_ * (1 + medium.Y_tilt)) +
                                  1;

                double tau_abs = 5.0 / 3 * consts::e * medium.rho(r_) / consts::mp * r_ / Bprime_ /
                                 pow(std::min(gamma_c_, gamma_m_), 5);

                double nu_m = synchrotron_nu(gamma_m_, Bprime_);
                double nu_c = synchrotron_nu(gamma_c_, Bprime_);
                double nu_M = synchrotron_nu(gamma_M_, Bprime_);
                // double nu_a = synchrotron_nu(gamma_a_, Bprime_);
                double nu_a = std::min(nu_c, nu_m) * pow(tau_abs, 3.0 / 5);
                if (i == 0) {
                    nv_m[j][k] = nu_m;
                    nv_c[j][k] = nu_c;
                    nv_a[j][k] = nu_a;
                    nv_M[j][k] = nu_M;
                }
                double emissivity_p = (pel - 1) / 2 * sqrt(3) * consts::e * consts::e * consts::e * Bprime_ /
                                      consts::me / consts::c / consts::c * 4 * Gamma_ * medium.rho(r_) / consts::mp *
                                      medium.xi;
                if (nu_prime >= nu_M) {
                    dFnu[i][j][k] = 0;
                } else {
                    dFnu[i][j][k] = dOmega * dr * r_ * r_ * doppler_ * doppler_ * emissivity_p *
                                    synchrotron_specific_emissivity(pel, nu_prime, nu_m, nu_c, nu_a);
                }
            }
        }
    }
    return std::make_tuple(dFnu, nv_m, nv_c, nv_a, nv_M);
}

void f() {
    double E_iso = 1e53 / 9e20 / 2e33;
    double Gamma0 = 300;
    double n_ism = 1 * 1.5e13 * 1.5e13 * 1.5e13;
    double eps_rad = 1;
    double eps_e = 0.1;
    double theta_j = 10.0 / 180 * consts::pi;

    // create model
    auto ism = create_ISM(n_ism, eps_rad, eps_e);
    auto inj = createIsoPowerLawInjection(0 / 9e20 / 2e33, 1000, 1, 0.1);
    auto blast = createTopHatJet(E_iso, Gamma0, theta_j, inj);
    //  auto blast = createPowerLawJet(E_iso, Gamma0, theta_j, 4, inj);
    // auto blast = createGaussianJet(E_iso, Gamma0, theta_j / 6, inj);

    // generate grid
    double M0 = E_iso / (Gamma0 * consts::c * consts::c);
    double R_ES = pow(3 * M0 / (4 * consts::pi * n_ism * consts::mp * Gamma0), 1.0 / 3);
    size_t grid_num = 1000;
    Array r = logspace(R_ES / 100, R_ES * 100, grid_num);
    Array theta = linspace(0, theta_j, 30);
    Array phi = linspace(0, 2 * consts::pi, 37);

    // solve dynamics
    auto [Gamma, U] = solveDynamics(blast, ism, r, theta);

    // calc observables
    double theta_obs = 0.0 / 180 * consts::pi;
    MeshGrid3d doppler = calc_doppler(Gamma, r, theta, phi, theta_obs);
    MeshGrid3d t_obs = calc_t_obs(Gamma, r, theta, phi, theta_obs);

    // calc radiation
    double nu_obs = 1e14 * 500;
    double pel = 2.2;
    auto [dFnu, nv_m, nv_c, nv_a, nv_M] =
        calc_Synchrotron_Radiation(r, theta, phi, Gamma, t_obs, doppler, ism, pel, nu_obs);

    // write to file
    auto r_c = boundary2centerlog(r);
    auto theta_c = boundary2center(theta);
    auto phi_c = boundary2center(phi);
    write2file(r_c, "r.txt");
    write2file(theta_c, "theta.txt");
    write2file(phi_c, "phi.txt");
    write2file(Gamma, "Gamma.txt");
    write2file(doppler, "doppler.txt");
    write2file(t_obs, "t_obs.txt");
    write2file(dFnu, "dFnu.txt");
    write2file(nv_m, "nu_m.txt");
    write2file(nv_c, "nu_c.txt");
    write2file(nv_a, "nu_a.txt");
    write2file(nv_M, "nu_Max.txt");
}

int main() {
    f();
    return 0;
}