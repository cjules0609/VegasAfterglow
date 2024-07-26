#ifndef _SYNCHROTRON_
#define _SYNCHROTRON_
#include <vector>

#include "medium.h"
#include "mesh.h"
#include "shock.h"

struct Y_IC {
    Y_IC(double nu_m, double nu_c, double B, double Y_T);
    Y_IC(double Y_T);
    Y_IC();
    double nu_hat_m{0};
    double nu_hat_c{0};
    double gamma_hat_m{0};
    double gamma_hat_c{0};
    double Y_T{0};
    size_t regime{0};

    double as_nu(double nu, double p) const;
    double as_gamma(double gamma, double p) const;

    static double Y_T_tilt(std::vector<Y_IC> const& Ys);
    static double Y_tilt_at_gamma(std::vector<Y_IC> const& Ys, double gamma, double p);
    static double Y_tilt_at_nu(std::vector<Y_IC> const& Ys, double nu, double p);
};

struct SynElectrons {
    // all in comoving frame
    double I_nu_peak{0};
    double gamma_m{0};
    double gamma_c{0};
    double gamma_a{0};
    double gamma_M{0};
    double p{2.3};
    double N_tot{0};
    double n_tot{0};
    double gamma_N_peak;
    size_t regime{0};
    double N(double gamma) const;
    double n(double gamma) const;

    double Y_c{0};
    std::vector<Y_IC> Ys;

   private:
    inline double gamma_spectrum_(double gamma) const;
};

struct SynPhotons {
    // all in comoving frame
    double L_nu_peak{0};
    double E_nu_peak{0};
    double nu_E_peak{0};
    double nu_m{0};
    double nu_c{0};
    double nu_a{0};
    double nu_M{0};
    double p{2.3};
    size_t regime{0};
    double L_nu(double nu) const;
    double E_nu(double nu) const;

    double Y_c{0};
    std::vector<Y_IC> Ys;

    void update_constant();

   private:
    double a_m_1_3{0};     // a_m_1_3 represents (nu_a / nu_m)^(1/3)
    double c_m_1_2{0};     // c_m_1_2 represents (nu_c / nu_m)^(1/2)
    double m_a_pa4_2{0};   // m_a_pa4_2 represents (nu_m / nu_a)^((p+4)/2)
    double a_m_mpa1_2{0};  // a_m_mpa1_2 represents (nu_a / nu_m)^((-p+1)/2)
    double a_c_1_3{0};     // a_c_1_3 represents (nu_a / nu_c)^(1/3)
    double a_m_1_2{0};     // a_m_1_2 represents (nu_a / nu_m)^(1/2)
    double R4{0};          // R coefficient in case4 in Bing Zhang's Book page 199
    double R5{0};          // R coefficient in case5 in Bing Zhang's Book page 200
    double R6{0};          // R coefficient in case6 in Bing Zhang's Book page 200
    inline double spectrum_(double nu) const;
};

using SynPhotonsArray = std::vector<SynPhotons>;

using SynPhotonsMesh = std::vector<std::vector<SynPhotons>>;

using SynElectronsArray = std::vector<SynElectrons>;

using SynElectronsMesh = std::vector<std::vector<SynElectrons>>;

SynPhotonsMesh create_syn_photons_grid(size_t theta_size, size_t r_size);

SynElectronsMesh create_syn_electrons_grid(size_t theta_size, size_t r_size);

SynElectronsMesh gen_syn_electrons(Coord const& coord, Shock const& shock);

SynPhotonsMesh gen_syn_photons(SynElectronsMesh const& electrons, Coord const& coord, Shock const& shock);

SynPhotonsMesh gen_syn_photons(Coord const& coord, Shock const& shock);

void update_electrons_4_Y(SynElectronsMesh& e, Shock const& shock);

double syn_gamma_c(double t_com, double B, std::vector<Y_IC> const& Ys, double p);

double syn_gamma_N_peak(double gamma_a, double gamma_m, double gamma_c);

double syn_nu(double gamma, double B);
#endif