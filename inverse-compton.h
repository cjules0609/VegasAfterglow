#ifndef _INVERSECOMPTON_
#define _INVERSECOMPTON_
#include <vector>

#include "medium.h"
#include "mesh.h"
#include "shock.h"
#include "synchrotron.h"
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

struct ICPhotons_ {
    template <typename Electrons, typename Photons>
    ICPhotons_(Electrons const& e, Photons const& ph);

    double I_nu(double nu) const;

   private:
    Array I_nu_;
    Array nu_IC_;
    double nu_max;
    double nu_min;
    double gamma_max;
    double gamma_min;
    size_t integral_resol{100};
    size_t spectrum_resol{100};
};

using ICPhotonArray = std::vector<ICPhoton>;
using ICPhotonMesh = std::vector<std::vector<ICPhoton>>;

double compton_sigma(double nu);

ICPhotonMesh create_IC_photon_grid(size_t theta_size, size_t r_size);

ICPhotonMesh gen_IC_photons(Coord const& coord, Shock const& shock, SynElectronsMesh const& electron,
                            SynPhotonsMesh const& photon, Medium const& medium);

MeshGrid solve_IC_Y_Thomson(Shock const& shock, SynElectronsMesh const& electron, Medium const& medium);

MeshGrid solve_IC_Y_KN(Shock const& shock, SynElectronsMesh const& electron, Medium const& medium);

double ICPhotons_::I_nu(double nu) const {
    if (nu <= nu_max) {
        return 0.0;
    } else {
        return interp_log_extra_lo(nu, nu_IC_, I_nu_);
    }
}

template <typename Electrons, typename Photons>
ICPhotons_::ICPhotons_(Electrons const& e, Photons const& ph) {
    double gamma_e_min = std::min(e.gamma_m, std::min(e.gamma_c, e.gamma_a));
    double nu_ph_min = std::min(ph.nu_m, std::min(ph.nu_c, ph.nu_a));

    double nu0_max = ph.nu_M * 10;
    double nu0_min = 4 * gamma_e_min * gamma_e_min * nu_ph_min;
    Array nu0_bin = logspace(nu0_min, nu0_max, integral_resol + 1);
    Array nu0 = boundary2center(nu0_bin);

    double gamma_min = 1;
    double gamma_max = e.gamma_M * 10;
    Array gamma_bin = logpace(gamma_min, gamma_max, integral_resol + 1);
    Array gamma = boundary2center(gamma_bin);

    this->nu_min = 4 * gamma_min * gamma_min * nu0_min;
    this->nu_max = 4 * gamma_max * gamma_max * nu0_max;

    MeshGrid I0 = create_grid(nu0.size(), gamma.size(), 0);
    MeshGrid I1 = create_grid(nu0.size(), gamma.size(), 0);
    MeshGrid I2 = create_grid(nu0.size(), gamma.size(), 0);

    for (size_t i = 0; i < nu0.size(); ++i) {
        for (size_t j = 0; j < gamma.size(); ++j) {
            double nu0_ = nu0[i];
            double gamma_ = gamma[j];
            double dS = fabs(nu0_bin[i + 1] - nu0_bin[i]) * fabs(gamma_bin[j + 1] - gamma_bin[j]);
            double f = 4 * gamma_ * gamma_ * nu0_;

            I0[i][j] = compton_sigma(nu0_ / gamma_) * e.dn_dgamma(gamma_) * ph.I_nu(nu0_) * dS / (f * nu0_);
            I1[i][j] = I0[i][j] * (1 - 2 * log(f)) / f;
            I2[i][j] = -2 * I0[i][j] / f * f;
        }
    }

    nu_IC_ = logspace(nu_min, nu_max, spectrum_resol + 1);
    I_nu_ = Array(spectrum_resol, 0);

    for (size_t k = 0; k < nu_IC_.size(); ++k) {
        for (size_t i = 0; i < nu0.size(); ++i) {
            for (size_t j = 0; j < gamma.size(); ++j) {
                double nu0_ = nu0[i];
                double gamma_ = gamma[j];
                if (nu0_ < nu_IC_[k] && nu_IC_[k] < 4 * gamma_ * gamma_ * nu0_) {
                    I_nu_[k] += nu_IC_[k] * I0[i][j];
                    I_nu_[k] += nu_IC_[k] * nu_IC_[k] * I1[i][j];
                    I_nu_[k] += nu_IC_[k] * nu_IC_[k] * nu_IC_[k] * I2[i][j];
                }
            }
        }
    }
}

#endif