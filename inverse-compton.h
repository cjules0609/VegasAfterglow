#ifndef _INVERSECOMPTON_
#define _INVERSECOMPTON_
#include <array>
#include <cmath>
#include <vector>

#include "macros.h"
#include "medium.h"
#include "mesh.h"
#include "shock.h"
#include "synchrotron.h"
#include "utilities.h"
/*
struct ICPhoton {
    double I_nu_peak{0};
    double nu_E_peak{0};

    double nu_mm{0};
    double nu_mc{0};
    double nu_ma{0};

    double nu_cm{0};
    double nu_cc{0};
    double nu_ca{0};

    double nu_am{0};
    double nu_ac{0};
    double nu_aa{0};

    double p{2.3};

    double I_nu(double nu) const;

   private:
    inline double I_nu_(double nu) const;
};
*/

inline const double IC_x0 = sqrt(2) / 3;

inline double compton_sigma(double nu) {
    double x = con::h * nu / (con::me * con::c2);
    if (x < 1e-2) {
        return con::sigmaT * (1 - 2 * x);
    } else if (x > 1e2) {
        return 3. / 8 * con::sigmaT * (log(2 * x) + 1.0 / 2) / x;
    } else {
        return 0.75 * con::sigmaT *
               ((1 + x) / (x * x * x) * (2 * x * (1 + x) / (1 + 2 * x) - log(1 + 2 * x)) + log(1 + 2 * x) / (2 * x) -
                (1 + 3 * x) / (1 + 2 * x) / (1 + 2 * x));
    }
}

struct IntegratorGrid {
    IntegratorGrid(double x_min, double x_max, double y_min, double y_max)
        : x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max) {
        logspace(x_min, x_max, x_bin);
        logspace(y_min, y_max, y_bin);
        boundary2center(x_bin, x);
        boundary2center(y_bin, y);
    }

    double x_min;
    double x_max;
    double y_min;
    double y_max;
    static constexpr size_t num{80};
    std::array<double, num + 1> x_bin{0};
    std::array<double, num + 1> y_bin{0};
    std::array<double, num> x{0};
    std::array<double, num> y{0};
    std::array<double, num> j_syn{0};
    std::array<double, num> ns{0};
    std::array<std::array<double, num>, num> I0{0};
};

struct ICPhoton {
    ICPhoton() = default;

    template <typename Electrons, typename Photons>
    void gen(Electrons const& e, Photons const& ph, double D_com) {
        double gamma_e_min = std::min(e.gamma_m, std::min(e.gamma_c, e.gamma_a));
        double nu_ph_min = std::min(ph.nu_m, std::min(ph.nu_c, ph.nu_a));

        double nu0_max = ph.nu_M * 10;
        double nu0_min = nu_ph_min / 1e5;

        double gamma_min = e.gamma_N_peak;
        double gamma_max = e.gamma_M * 10;

        IntegratorGrid grid(nu0_min, nu0_max, gamma_min, gamma_max);

        for (size_t i = 0; i < grid.num; i++) {
            grid.j_syn[i] = ph.j_nu(grid.x[i]);
        }

        for (size_t i = 0; i < grid.num; i++) {
            grid.ns[i] = e.n(grid.y[i]);
        }

        for (size_t i = 0; i < grid.num; ++i) {
            double nu0_ = grid.x[i];
            double dnu = grid.x_bin[i + 1] - grid.x_bin[i];
            for (size_t j = 0; j < grid.num; ++j) {
                double gamma_ = grid.y[j];
                double dgamma = grid.y_bin[j + 1] - grid.y_bin[j];
                double dS = fabs(dnu * dgamma);
                double f = 4 * gamma_ * gamma_ * nu0_ * nu0_;
                grid.I0[i][j] = grid.ns[j] * grid.j_syn[i] * dS / f * compton_sigma(nu0_ / gamma_);
            }
        }

        double nu_min = 4 * gamma_min * gamma_min * nu0_min;
        double nu_max = 4 * gamma_max * gamma_max * nu0_max;

        nu_IC_ = logspace(nu_min, nu_max, spectrum_resol);
        j_nu_ = Array(spectrum_resol, 0);

        for (size_t k = 0; k < nu_IC_.size(); ++k) {
            for (size_t i = 0; i < grid.num; ++i) {
                double nu0_ = grid.x[i];
                for (size_t j = 0; j < grid.num; ++j) {
                    double gamma_ = grid.y[j];
                    if (nu0_ < nu_IC_[k] && nu_IC_[k] < 4 * gamma_ * gamma_ * nu0_ * IC_x0) {
                        j_nu_[k] += grid.I0[i][j] * nu_IC_[k];
                    }
                }
            }
        }

        for (size_t k = 0; k < nu_IC_.size(); ++k) {
            j_nu_[k] *= D_com;
        }
    };

    double j_nu(double nu) const;

   public:
    static constexpr size_t spectrum_resol{50};

   private:
    Array j_nu_;
    Array nu_IC_;
};

using ICPhotonArray = std::vector<ICPhoton>;

using ICPhotonMesh = std::vector<std::vector<ICPhoton>>;

ICPhotonMesh create_IC_photon_grid(size_t theta_size, size_t r_size);

ICPhotonMesh gen_IC_photons(SynElectronsMesh const& electron, SynPhotonsMesh const& photon, Shock const& shock);

MeshGrid solve_IC_Y_Thomson(SynElectronsMesh const& electron, Shock const& shock, Medium const& medium);

MeshGrid solve_IC_Y_KN(SynElectronsMesh const& electron, Shock const& shock, Medium const& medium);

#endif