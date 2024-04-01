#include <algorithm>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

struct Indexes {
    size_t i, j, k;
    bool operator<(Indexes const& other) const { return i < other.i; }
};

using MeshGrid = std::vector<std::vector<double>>;
using MeshGrid3d = std::vector<std::vector<std::vector<double>>>;
using EATsurface = std::vector<std::pair<double, Indexes>>;
using Array = std::vector<double>;
using Profile = std::function<double(double)>;
using Profile2d = std::function<double(double, double)>;

static Profile2d noInjection = Profile2d([](double theta, double t_lab) { return 0; });

template <typename Fun>
auto root_bisection(Fun f, decltype(f(0)) low, decltype(f(0)) high, decltype(f(0)) eps = 1e-6) -> decltype(f(0)) {
    using Scalar = decltype(f(0));

    for (; fabs((high - low)) > fabs(high) * eps;) {
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
    constexpr double c2 = c * c;
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
class Coord {
   public:
    Coord(Array r_b, Array theta_b, Array phi_b) : r_b(r_b), theta_b(theta_b), phi_b(phi_b) {
        r = boundary2centerlog(r_b);
        theta = boundary2center(theta_b);
        phi = boundary2center(phi_b);
    }
    Coord() = delete;
    Array r_b;
    Array theta_b;
    Array phi_b;
    Array r;
    Array theta;
    Array phi;
};
class Medium {
   public:
    Medium(Profile rho, Profile mass, double eps_rad, double eps_e, double eps_B, double Y_tilt = 0, double xi = 1)
        : rho(rho), mass(mass), eta_rad(eps_rad), eps_e(eps_e), eps_B(eps_B), Y_tilt{Y_tilt}, xi(xi){};
    double eta_rad{1};
    double eps_e{0.1};
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
        double t_eng = y[2];
        // double t_com= y[3];//comoving time
        dydr[0] = dGammadr(r, Gamma, u, t_eng);
        dydr[1] = dUdr(r, Gamma, u, t_eng);
        dydr[2] = dtdr_eng(Gamma);
        dydr[3] = dtdr_com(Gamma);
    };

   private:
    inline double dGammadr(double r, double Gamma, double u, double t_eng) {
        double gamma_eos = (4 * Gamma + 1) / (3 * Gamma);
        double dm = medium.mass(r) * dOmega / (4 * consts::pi);
        double dM0 = blast.dEdOmega(theta, t_eng) * dOmega / (blast.Gamma0(theta) * consts::c2);
        double Gamma2 = Gamma * Gamma;
        double a1 = dOmega * r * r * medium.rho(r) / Gamma * (Gamma2 - 1) * (gamma_eos * Gamma - gamma_eos + 1);
        double a2 = -(gamma_eos - 1) / Gamma * (gamma_eos * Gamma2 - gamma_eos + 1) * 3 * u / r;
        double b1 = (dM0 + dm) * consts::c2;
        double b2 = (gamma_eos * gamma_eos * (Gamma2 - 1) + 3 * gamma_eos - 2) * u / Gamma2;
        return -(a1 + a2) / (b1 + b2);
    };

    inline double dUdr(double r, double Gamma, double u, double t_eng) {
        double gamma_eos = (4 * Gamma + 1) / (3 * Gamma);
        double E = dOmega * r * r * medium.rho(r) * consts::c2;
        return (1 - medium.eps_e * medium.eta_rad) * (Gamma - 1) * E -
               (gamma_eos - 1) * (3 / r - dGammadr(r, Gamma, u, t_eng) / Gamma) * u;
    };

    inline double dtdr_eng(double Gamma) {
        double Gb = std::sqrt(Gamma * Gamma - 1);
        return (Gamma - Gb) / (Gb * consts::c);
    };

    inline double dtdr_com(double Gamma) { return 1 / (sqrt(Gamma * Gamma - 1) * consts::c); };  // co-moving time

    Medium medium;
    BlastWave blast;
    double theta_lo;
    double theta_hi;
    double theta;
    double dOmega;
};

template <typename Eqn>
void solveShell(Array const& r, Array& Gamma, Array& u, Array& t_eng, Array& t_com, Eqn const& eqn) {
    using namespace boost::numeric::odeint;
    double atol = 0;
    double rtol = 1e-9;
    auto stepper = bulirsch_stoer_dense_out<std::vector<double>>{atol, rtol};
    Array state{0, 0, 0, 0};
    double dr = (r[1] - r[0]) / 1000;
    int i = 0;
    stepper.initialize(Array{Gamma[0], u[0], t_eng[0], t_com[0]}, r[0], dr);
    for (; stepper.current_time() <= r.back();) {
        auto [t_last, t_current] = stepper.do_step(eqn);
        for (; stepper.current_time() > r[i + 1] && i + 1 < r.size();) {
            i++;
            stepper.calc_state(r[i], state);
            Gamma[i] = state[0];
            u[i] = state[1];
            t_eng[i] = state[2];
            t_com[i] = state[3];
        }
    }
}

auto solveDynamics(BlastWave const& blast, Medium const& medium, Coord const& coord) {
    MeshGrid Gamma = createGrid(coord.r, coord.theta, 1);
    MeshGrid U = createGrid(coord.r, coord.theta, 0);
    MeshGrid t_eng = createGrid(coord.r, coord.theta, 0);
    MeshGrid t_com = createGrid(coord.r, coord.theta, 0);

    for (size_t i = 0; i < coord.theta.size(); ++i) {
        double dcos = std::fabs(std::cos(coord.theta_b[i + 1]) - std::cos(coord.theta_b[i]));
        double dphi = 2 * consts::pi;
        double dOmega = 2 * dphi * dcos;  // bipolar outflow
        Gamma[i][0] = blast.Gamma0(coord.theta[i]);
        double beta0 = sqrt(Gamma[i][0] * Gamma[i][0] - 1) / Gamma[i][0];
        U[i][0] = medium.mass(coord.r[0]) * dOmega / (4 * consts::pi) * consts::c2 * (Gamma[i][0] - 1);
        t_eng[i][0] = coord.r[0] * (1 - beta0) / beta0 / consts::c;
        t_com[i][0] = coord.r[0] / sqrt(Gamma[i][0] * Gamma[i][0] - 1) / consts::c;
        auto eqn = DynamicsEqn(medium, blast, coord.theta_b[i], coord.theta_b[i + 1]);
        solveShell(coord.r, Gamma[i], U[i], t_eng[i], t_com[i], eqn);
    }
    return std::make_tuple(Gamma, U, t_eng, t_com);
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
inline double shock_width_com(double r, double Gamma) { return r / Gamma / 12; }

MeshGrid3d calc_t_obs(Coord const& coord, MeshGrid const& Gamma, double theta_obs) {
    MeshGrid3d t_obs = createGrid3d(coord.r, coord.theta, coord.phi);
    using namespace boost::numeric::odeint;
    double atol = 0;
    double rtol = 1e-9;
    auto stepper = bulirsch_stoer_dense_out<double>{atol, rtol};
    TobsEqn eqn;
    eqn.r_ = coord.r;
    double dr0 = (coord.r_b[1] - coord.r_b[0]) / 1000;

    for (size_t j = 0; j < coord.theta.size(); ++j) {
        double theta_ = coord.theta[j];
        eqn.Gamma_ = Gamma[j];
        double Gamma0 = Gamma[j][0];
        double beta0 = sqrt(1 - 1 / Gamma0 / Gamma0);
        for (size_t i = 0; i < coord.phi.size(); ++i) {
            double phi_ = coord.phi[i];
            eqn.cos_ = sin(theta_) * cos(phi_) * sin(theta_obs) + cos(theta_) * cos(theta_obs);
            double t_obs0 = coord.r[0] * (1 - beta0 * eqn.cos_) / consts::c / beta0;
            t_obs[i][j][0] = t_obs0;
            int k = 0;
            stepper.initialize(t_obs0, coord.r[0], dr0);
            for (; stepper.current_time() <= coord.r.back();) {
                auto [t_last, t_current] = stepper.do_step(eqn);
                for (; stepper.current_time() > coord.r[k + 1] && k + 1 < coord.r.size();) {
                    k++;
                    stepper.calc_state(coord.r[k], t_obs[i][j][k]);
                }
            }
        }
    }
    return t_obs;
}

MeshGrid3d calc_doppler(Coord const& coord, MeshGrid const& Gamma, double theta_obs) {
    MeshGrid3d doppler = createGrid3d(coord.r, coord.theta, coord.phi);
    for (size_t i = 0; i < coord.phi.size(); ++i) {
        double phi_ = coord.phi[i];
        for (size_t j = 0; j < coord.theta.size(); ++j) {
            double theta_ = coord.theta[j];
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

inline double syn_nu(double gamma, double B) {
    double nu = 3 / (4 * consts::pi) * consts::e * B / (consts::me * consts::c) * gamma * gamma;
    return nu;
}

inline double syn_I_nu(double pel, double nu, double nu_m, double nu_c, double nu_a) {
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
            ;
        }
    } else if (nu_c < nu_m && nu_m < nu_a) {
        if (nu <= nu_a) {
            return pow(nu / nu_a, 2);
        } else {
            return pow(nu / nu_a, -pel / 2) * pow(nu_m / nu_a, (pel - 1) / 2) * pow(nu_c / nu_a, 1.0 / 2) / 3;
        }
    }
}

class SynFreq {
   public:
    MeshGrid nu_m;
    MeshGrid nu_c;
    MeshGrid nu_a;
    MeshGrid nu_M;
};

auto calc_B_field(Coord const& coord, MeshGrid const& Gamma, Medium const& medium) {
    MeshGrid B = createGrid(coord.r, coord.theta);

    for (size_t j = 0; j < coord.theta.size(); ++j) {
        double theta_ = coord.theta[j];
        for (size_t k = 0; k < coord.r.size(); ++k) {
            double dr = coord.r_b[k + 1] - coord.r_b[k];
            double r_ = coord.r[k];
            double Gamma_ = Gamma[j][k];

            B[j][k] =
                sqrt(8 * consts::pi * medium.eps_B * medium.rho(r_) * (4 * Gamma_ * Gamma_ - 4 * Gamma_)) * consts::c;
        }
    }

    return B;
}

auto calc_syn_frequencies(Coord const& coord, MeshGrid const& Gamma, MeshGrid const& t_com, MeshGrid const& B_prime,
                          MeshGrid const& I_syn_peak, Medium const& medium, double pel) {
    if (pel <= 1) {
        std::cout << "pel should be larger than 1" << std::endl;
        exit(1);
    }

    SynFreq syn_freq{createGrid(coord.r, coord.theta), createGrid(coord.r, coord.theta),
                     createGrid(coord.r, coord.theta), createGrid(coord.r, coord.theta)};

    for (size_t j = 0; j < coord.theta.size(); ++j) {
        double theta_ = coord.theta[j];
        for (size_t k = 0; k < coord.r.size(); ++k) {
            // max injection frequency
            double Bprime_ = B_prime[j][k];
            double gamma_M_ =
                sqrt(6 * consts::pi * consts::e / (consts::sigmaT * Bprime_ * medium.zeta * (1 + medium.Y_tilt)));
            double nu_M = syn_nu(gamma_M_, Bprime_);

            // min injection frequency
            double Gamma_ = Gamma[j][k];
            double gamma_m_ = 1;
            if (pel > 2) {
                gamma_m_ = (pel - 2) / (pel - 1) * medium.eps_e * (Gamma_ - 1) * 1836 / medium.xi + 1;
            } else if (pel < 2) {
                gamma_m_ =
                    pow((2 - pel) / (pel - 1) * medium.eps_e * (Gamma_ - 1) * 1836 / medium.xi * pow(gamma_M_, pel - 1),
                        1 / (pel - 1)) +
                    1;
            } else {
                gamma_m_ = root_bisection(
                    [=](double x) -> double {
                        return (x * log(gamma_M_) - (x + 1) * log(x) - medium.eps_e * (Gamma_ - 1) * 1836 / medium.xi -
                                log(gamma_M_));
                    },
                    1, gamma_M_);
            }
            double nu_m = syn_nu(gamma_m_, Bprime_);

            // cooling frequency
            double t_com_ = t_com[j][k];
            double beta_ = sqrt(Gamma_ * Gamma_ - 1) / Gamma_;
            double gamma_c_ = 6 * consts::pi * consts::me * consts::c /
                                  (consts::sigmaT * t_com_ * Bprime_ * Bprime_ * (1 + medium.Y_tilt)) +
                              1;

            double nu_c = syn_nu(gamma_c_, Bprime_);

            // self-absorption frequency
            double I_syn_peak_ = I_syn_peak[j][k];
            double nu_peak = std::min(nu_m, nu_c);
            double gamma_peak = std::min(gamma_m_, gamma_c_);
            double gamma_eos = (4 * Gamma_ + 1) / (3 * Gamma_);
            double kT = (gamma_peak - 1) * consts::me * consts::c2 * (gamma_eos - 1);
            double nu_a = pow(I_syn_peak_ * consts::c2 / pow(nu_peak, 1.0 / 3) / kT / 2,
                              3.0 / 5);  // 2kT(nv_a/c)^2 = I_peak (nu_a/nu_peak)^(1/3)
            if (nu_a > nu_peak) {        // kT = (gamma_a-1) * consts::me * consts::c2*(gamma_eos-1); peak at nu_a
                /*nu_a = pow(I_syn_peak_ / consts::me / 2 / (gamma_eos - 1) /
                               sqrt(4 * consts::pi / 3 * consts::me * consts::c / consts::e / Bprime_),
                           2.0 / 5);*/
                double A = sqrt(4 * consts::pi / 3 * consts::me * consts::c / consts::e / Bprime_);
                double B = I_syn_peak_ / (2 * consts::me * (gamma_eos - 1));
                nu_a = root_bisection([=](double x) -> double { return A * x * x * x * x * x - x * x * x * x - B; },
                                      sqrt(nu_peak), sqrt(nu_M));
                nu_a *= nu_a;
            }

            syn_freq.nu_m[j][k] = nu_m;
            syn_freq.nu_c[j][k] = nu_c;
            syn_freq.nu_a[j][k] = nu_a;
            syn_freq.nu_M[j][k] = nu_M;
        }
    }

    return syn_freq;
}

auto calc_syn_I_nu_obs(Coord const& coord, SynFreq const& syn_freq, MeshGrid const& I_nu_peak,
                       MeshGrid3d const& doppler, double pel, double nu_obs) {
    if (pel < 1) {
        std::cout << "pel should be larger than 1" << std::endl;
        exit(1);
    }

    MeshGrid3d I_nu_syn_obs = createGrid3d(coord.r, coord.theta, coord.phi);

    for (size_t i = 0; i < coord.phi.size(); ++i) {
        for (size_t j = 0; j < coord.theta.size(); ++j) {
            for (size_t k = 0; k < coord.r.size(); ++k) {
                double doppler_ = doppler[i][j][k];
                double I_nu_peak_ = I_nu_peak[j][k];
                double nu_prime = nu_obs / doppler_;

                if (nu_prime >= syn_freq.nu_M[j][k]) {
                    I_nu_syn_obs[i][j][k] = 0;
                } else {
                    I_nu_syn_obs[i][j][k] =
                        doppler_ * doppler_ * doppler_ * I_nu_peak_ *
                        syn_I_nu(pel, nu_prime, syn_freq.nu_m[j][k], syn_freq.nu_c[j][k], syn_freq.nu_a[j][k]);
                }
            }
        }
    }
    return I_nu_syn_obs;
}

auto calc_syn_I_nu_peak(Coord const& coord, MeshGrid const& Gamma, MeshGrid const& Bprime, Medium const& medium,
                        double pel) {
    if (pel < 1) {
        std::cout << "pel should be larger than 1" << std::endl;
        exit(1);
    }

    MeshGrid I_peak = createGrid(coord.r, coord.theta);

    for (size_t j = 0; j < coord.theta.size(); ++j) {
        for (size_t k = 0; k < coord.r.size(); ++k) {
            double r_ = coord.r[k];
            double Gamma_ = Gamma[j][k];
            double Bprime_ = Bprime[j][k];

            I_peak[j][k] = (pel - 1) / 2 * sqrt(3) * consts::e * consts::e * consts::e * Bprime_ /
                           (consts::me * consts::c2) * 4 * Gamma_ * (medium.rho(r_) / consts::mp) * medium.xi *
                           shock_width_com(r_, Gamma_);
        }
    }

    return I_peak;
}

auto get_sorted_EAT_surface(MeshGrid3d const& t_obs) {
    EATsurface eat_s;
    eat_s.reserve(t_obs.size() * t_obs[0].size() * t_obs[0][0].size());

    for (size_t i = 0; i < t_obs.size(); ++i) {
        for (size_t j = 0; j < t_obs[0].size(); ++j) {
            for (size_t k = 0; k < t_obs[0][0].size(); ++k) {
                eat_s.emplace_back(t_obs[i][j][k], Indexes{i, j, k});
            }
        }
    }
    std::sort(eat_s.begin(), eat_s.end());
    return eat_s;
}

Array calc_light_curve(Coord const& coord, MeshGrid3d const& I_nu_obs, EATsurface const& eat_s, Array const& t_bin) {
    Array L_nu = zeros(t_bin.size() - 1);

    for (size_t i = 0, j = 0; i < eat_s.size(); ++i) {
        double t_ = eat_s[i].first;
        size_t phi_idx = eat_s[i].second.i;
        size_t theta_idx = eat_s[i].second.j;
        size_t r_idx = eat_s[i].second.k;

        double dcos = std::cos(coord.theta_b[theta_idx + 1]) - std::cos(coord.theta_b[theta_idx]);
        double dphi = coord.phi_b[phi_idx + 1] - coord.phi_b[phi_idx];
        double dOmega = std::fabs(dphi * dcos);
        double r_ = coord.r[r_idx];

        double dL_nu_ = r_ * r_ * dOmega * I_nu_obs[phi_idx][theta_idx][r_idx];

        if (t_ > t_bin[j + 1]) {
            j++;
        }
        if (j < L_nu.size()) {
            L_nu[j] += dL_nu_;
        } else {
            break;
        }
    }
    for (size_t i = 0; i < L_nu.size(); ++i) {
        L_nu[i] /= (t_bin[i + 1] - t_bin[i]);
    }
    return L_nu;
}

auto afterglow_gen(Coord const& coord, BlastWave const& blast, Medium const& medium, double theta_obs,
                   Array const& nu_obs, size_t light_curve_resolution = 100) {
    // solve dynamics
    auto [Gamma, U, t_eng, t_com] = solveDynamics(blast, medium, coord);

    // calc radiation
    double pel = 2.2;
    auto B = calc_B_field(coord, Gamma, medium);
    auto I_syn_peak = calc_syn_I_nu_peak(coord, Gamma, B, medium, pel);
    auto syn_freq = calc_syn_frequencies(coord, Gamma, t_com, B, I_syn_peak, medium, pel);

    // calc observables
    MeshGrid3d doppler = calc_doppler(coord, Gamma, theta_obs);
    MeshGrid3d t_obs = calc_t_obs(coord, Gamma, theta_obs);
    EATsurface eat_s = get_sorted_EAT_surface(t_obs);

    // calc synchrotron radiation
    std::vector<MeshGrid3d> I_syn_obs;
    std::vector<Array> L_nu;

    for (auto nu_obs_ : nu_obs) {
        I_syn_obs.emplace_back(calc_syn_I_nu_obs(coord, syn_freq, I_syn_peak, doppler, pel, nu_obs_));
    }

    // generate light curve
    Array t_bin = logspace(eat_s.front().first, eat_s.back().first, light_curve_resolution);
    // std::cout << eat_s.front().first << ' ' << eat_s.back().first;
    Array t_lc = boundary2centerlog(t_bin);
    for (auto& I : I_syn_obs) {
        L_nu.emplace_back(calc_light_curve(coord, I, eat_s, t_bin));
    }
    return std::make_tuple(Gamma, U, t_eng, t_com, t_obs, doppler, I_syn_obs, syn_freq, t_lc, L_nu);
}

int main() {
    double E_iso = 1e53 / 9e20 / 2e33;
    double Gamma0 = 300;
    double n_ism = 1 * 1.5e13 * 1.5e13 * 1.5e13;
    double eps_rad = 1;
    double eps_e = 0.1;
    double theta_j = 5.0 / 180 * consts::pi;

    // create model
    auto ism = create_ISM(n_ism, eps_rad, eps_e);
    auto inj = createIsoPowerLawInjection(0 / 9e20 / 2e33, 1000, 1, 2);
    auto blast = createTopHatJet(E_iso, Gamma0, theta_j, inj);
    // auto blast = createPowerLawJet(E_iso, Gamma0, theta_j, 4, inj);
    // auto blast = createGaussianJet(E_iso, Gamma0, theta_j / 6, inj);

    // generate grid
    double M0 = E_iso / (Gamma0 * consts::c * consts::c);
    double R_ES = pow(3 * M0 / (4 * consts::pi * n_ism * consts::mp * Gamma0), 1.0 / 3);
    size_t grid_num = 300;
    Array r = logspace(R_ES / 100, R_ES * 100, grid_num);
    Array theta = linspace(0, theta_j, 120);
    Array phi = linspace(0, 2 * consts::pi, 2);
    Coord coord{r, theta, phi};

    // specify observables
    double nu_obs_R = 1e9 * 500;
    double nu_obs_O = 1e14 * 500;
    double nu_obs_X = 1e20 * 500;
    Array nu_obs{nu_obs_R, nu_obs_O, nu_obs_X};

    double theta_obs = 0.0 / 180 * consts::pi;

    size_t light_curve_resolution = 100;

    auto [Gamma, U, t_eng, t_com, t_obs, doppler, em_syn_obs, syn_freq, t_lc, Lnu] =
        afterglow_gen(coord, blast, ism, theta_obs, nu_obs, light_curve_resolution);

    // write to file
    std::string prefix = "off-axis-tophat/";

    write2file(coord.r, prefix + "r.txt");
    write2file(coord.theta, prefix + "theta.txt");
    write2file(coord.phi, prefix + "phi.txt");
    write2file(Gamma, prefix + "Gamma.txt");
    write2file(t_eng, prefix + "t_eng.txt");
    write2file(t_com, prefix + "t_com.txt");
    write2file(doppler, prefix + "doppler.txt");
    write2file(t_obs, prefix + "t_obs.txt");
    write2file(t_lc, prefix + "t_lc.txt");
    write2file(syn_freq.nu_m, prefix + "nu_m.txt");
    write2file(syn_freq.nu_c, prefix + "nu_c.txt");
    write2file(syn_freq.nu_a, prefix + "nu_a.txt");
    write2file(syn_freq.nu_M, prefix + "nu_Max.txt");

    for (size_t i = 0; i < nu_obs.size(); ++i) {
        write2file(em_syn_obs[i], prefix + "I_nu_" + std::to_string(int(log10(nu_obs[i] / 500))) + ".txt");
        write2file(Lnu[i], prefix + "L_nu_" + std::to_string(int(log10(nu_obs[i] / 500))) + ".txt");
    }
    return 0;
}