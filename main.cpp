#include <algorithm>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

using MeshGrid = std::vector<std::vector<double>>;
using Array = std::vector<double>;
using Profile = std::function<double(double)>;
using Profile2d = std::function<double(double, double)>;

static Profile2d noInjection = Profile2d([](double theta, double t_lab) { return 0; });

namespace consts {
    constexpr double c = 1;
    constexpr double mp = 1.67e-24 / 2e33;
    constexpr double me = mp / 1836;
    // constexpr double kB = 1.38e-16;
    // constexpr double e = 4.8e-10;
    constexpr double pi = 3.14159265358979323846;
}  // namespace consts

MeshGrid createGrid(Array const& r, Array const& theta, double val = 0) {
    return MeshGrid(theta.size(), Array(r.size(), val));
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
    Medium(Profile rho, Profile mass, double eps_rad, double eps_e)
        : rho(rho), mass(mass), eps_rad(eps_rad), eps_e(eps_e){};
    double eps_rad{1};
    double eps_e{0.3};
    Profile rho;
    Profile mass;
};

class DynamicsEqn {
   public:
    DynamicsEqn(Medium medium, BlastWave blast, double theta_lo, double theta_hi, double theta_obs,
                double gamma_eos = 4.0 / 3)
        : medium(medium),
          blast(blast),
          theta_lo(theta_lo),
          theta_hi(theta_hi),
          gamma_eos(gamma_eos),
          theta(0.5 * (theta_lo + theta_hi)),
          theta_obs(theta_obs),
          dOmega(4 * consts::pi * std::fabs(std::cos(theta_hi) - std::cos(theta_lo))){};  // bipolar outflow

    void operator()(Array const& y, Array& dydr, double r) {
        double Gamma = y[0];
        double u = y[1];
        double t_obs = y[2];
        double t_lab = y[3];
        dydr[0] = dGammadr(r, Gamma, u, t_lab);
        dydr[1] = dUdr(r, Gamma, u, t_lab);
        dydr[2] = dt_obsdr(Gamma);
        dydr[3] = dt_labdr(Gamma);
    };

   private:
    inline double dGammadr(double r, double Gamma, double u, double t_lab) {
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
        double E = dOmega * r * r * medium.rho(r) * consts::c * consts::c;
        return (1 - medium.eps_e * medium.eps_rad) * (Gamma - 1) * E -
               (gamma_eos - 1) * (3 / r - dGammadr(r, Gamma, u, t_lab) / Gamma) * u;
    };

    inline double dt_obsdr(double Gamma) {
        double sqr = std::sqrt(Gamma * Gamma - 1);
        return (Gamma - sqr * std::cos(theta_obs - theta)) / (sqr * consts::c);
    };

    inline double dt_labdr(double Gamma) { return Gamma / sqrt(Gamma * Gamma - 1) / consts::c; };

    double gamma_eos{4.0 / 3};
    Medium medium;
    BlastWave blast;
    double theta_lo;
    double theta_hi;
    double theta;
    double theta_obs;
    double dOmega;
};

template <typename Eqn>
void solveShell(Array const& r, Array& Gamma, Array& u, Array& t_obs, Eqn const& eqn) {
    using namespace boost::numeric::odeint;
    double atol = 0;
    double rtol = 1e-9;
    auto stepper = bulirsch_stoer_dense_out<std::vector<double>>{atol, rtol};

    Array state{0, 0, 0};
    double dr = (r[1] - r[0]) / 1000;
    int i = 0;
    stepper.initialize(Array{Gamma[0], u[0], t_obs[0], r[0] / consts::c}, r[0], dr);
    for (; stepper.current_time() <= r.back();) {
        auto [t_last, t_current] = stepper.do_step(eqn);
        for (; stepper.current_time() > r[i + 1] && i + 1 < r.size();) {
            i++;
            stepper.calc_state(r[i], state);
            Gamma[i] = state[0];
            u[i] = state[1];
            t_obs[i] = state[2];
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

auto solveDynamics(BlastWave const& blast, Medium const& medium, Array const& r, Array& theta, double theta_obs) {
    Array theta_c = boundary2center(theta);
    MeshGrid Gamma = createGrid(r, theta_c, 1);
    MeshGrid U = createGrid(r, theta_c, 0);
    MeshGrid t_obs = createGrid(r, theta_c, 0);

    for (size_t i = 0; i < theta_c.size(); ++i) {
        double dOmega = 4 * consts::pi * std::fabs(std::cos(theta[i + 1]) - std::cos(theta[i]));  // bipolar outflow

        Gamma[i][0] = blast.Gamma0(theta_c[i]);
        U[i][0] = medium.mass(r[0]) * dOmega / (4 * consts::pi) * consts::c * consts::c * (Gamma[i][0] - 1);
        t_obs[i][0] = 0;
        auto eqn = DynamicsEqn(medium, blast, theta[i], theta[i + 1], theta_obs);
        solveShell(r, Gamma[i], U[i], t_obs[i], eqn);
    }
    return std::make_tuple(theta_c, Gamma, U, t_obs);
}

void printArray(Array const& arr) {
    for (auto const& a : arr) {
        std::cout << a << " ";
    }
    std::cout << std::endl;
}

void write2file(MeshGrid const& grid, std::string const& filename) {
    std::ofstream file(filename);
    for (size_t j = 0; j < grid[0].size(); ++j) {
        for (size_t i = 0; i < grid.size(); ++i) {
            file << grid[i][j] << " ";
        }
        file << std::endl;
    }
}

void write2file(Array const& arry, std::string const& filename) {
    std::ofstream file(filename);
    for (size_t i = 0; i < arry.size(); ++i) {
        file << arry[i] << " ";
    }
}

Medium create_ISM(double n_ism, double eps_rad, double eps_e) {
    auto rho = [=](double r) { return n_ism * consts::mp; };
    auto m = [=](double r) { return 4 * consts::pi / 3 * r * r * r * n_ism * consts::mp; };
    return Medium{rho, m, eps_rad, eps_e};
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
    auto Gamma = [=](double theta) { return theta < theta_c ? Gamma0 : 1; };
    return BlastWave{dEdOmega, Gamma};
}

BlastWave createPowerLawJet(double E_iso_on_axis, double Gamma0_on_axis, double theta_m, double k,
                            Profile2d inject = noInjection) {
    double e0 = E_iso_on_axis / (4 * consts::pi);
    auto dEdOmega = [=](double theta, double t_lab) {
        return (theta < theta_m ? e0 : e0 * std::pow(theta / theta_m, -k)) + inject(theta, t_lab);
    };
    auto Gamma = [=](double theta) { return Gamma0_on_axis; };
    return BlastWave{dEdOmega, Gamma};
}

BlastWave createGaussianJet(double E_iso_on_axis, double Gamma0_on_axis, double theta_c,
                            Profile2d inject = noInjection) {
    double e0 = E_iso_on_axis / (4 * consts::pi);
    auto dEdOmega = [=](double theta, double t_lab) {
        return e0 * std::exp(-theta * theta / (2 * theta_c * theta_c)) + inject(theta, t_lab);
    };
    auto Gamma = [=](double theta) { return Gamma0_on_axis; };
    return BlastWave{dEdOmega, Gamma};
}

Profile2d createIsoPowerLawInjection(double L0, double t0, double t_wait, double q) {
    double t0q = std::pow(t0, q);
    double tw1q = std::pow(t0, 1 - q);
    double Omega = 4 * consts::pi;
    if (q < 1) {
        return [=](double theta, double t_lab) {
            if (t_lab < t_wait) {
                return 0.0;
            } else {
                return L0 * t0q / (1 - q) * (std::pow(t_lab + t0 - t_wait, 1 - q) - tw1q) / Omega;
            }
        };
    } else if (q > 1) {
        return [=](double theta, double t_lab) {
            if (t_lab < t_wait) {
                return 0.0;
            } else {
                return L0 * t0q / (q - 1) * (tw1q - std::pow(t_lab + t0 - t_wait, 1 - q)) / Omega;
            }
        };
    } else {
        return [=](double theta, double t_lab) {
            if (t_lab < t_wait) {
                return 0.0;
            } else {
                return L0 * t0 * std::log((t_lab + t0 - t_wait) / t_wait) / Omega;
            }
        };
    }
}

void f() {
    double E_iso = 1e51 / 9e20 / 2e33;
    double Gamma0 = 300;
    double n_ism = 1 * 1.5e13 * 1.5e13 * 1.5e13;
    double eps_rad = 1;
    double eps_e = 0.3;
    double theta_j = 35.0 / 180 * consts::pi;

    auto ism = create_ISM(n_ism, eps_rad, eps_e);
    auto inj = createIsoPowerLawInjection(1e51 / 9e20 / 2e33, 1000, 1, 0.5);
    auto blast = createTopHatJet(E_iso, Gamma0, theta_j, inj);

    double M0 = E_iso / (Gamma0 * consts::c * consts::c);

    double R_ES = pow(3 * M0 / (4 * consts::pi * n_ism * consts::mp * Gamma0), 1.0 / 3);

    size_t grid_num = 300;
    Array r = logspace(R_ES / 10, R_ES * 1000, grid_num);
    Array theta = linspace(0, theta_j, 2);

    double theta_obs = 0.0 / 180 * consts::pi;

    auto [theta_c, Gamma, U, t_obs] = solveDynamics(blast, ism, r, theta, theta_obs);

    write2file(r, "r.txt");
    write2file(theta_c, "theta.txt");
    write2file(t_obs, "t_obs.txt");
    write2file(Gamma, "Gamma.txt");
    write2file(U, "U.txt");
}

int main() {
    f();
    return 0;
}