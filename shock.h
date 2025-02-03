#ifndef _FSDYNAMICS_
#define _FSDYNAMICS_

#include <tuple>

#include "jet.h"
#include "medium.h"
#include "mesh.h"
#include "physics.h"
class Shock {
   public:
    Shock(size_t phi_size, size_t theta_size, size_t r_size, double eps_e, double eps_B);
    Shock() = delete;

    MeshGrid3d t_com;           // comoving time
    MeshGrid3d t_eng;           // engine time
    MeshGrid3d Gamma_rel;       // relative lorentz factor between down stream and up stream
    MeshGrid3d B;               // comoving magnetic field
    MeshGrid3d column_num_den;  // down stream proton column density
    double eps_e{0};
    double eps_B{0};

    auto shape() const { return std::make_tuple(phi_size, theta_size, r_size); }

   private:
    size_t phi_size{0};
    size_t theta_size{0};
    size_t r_size{0};
};

class ForwardShockEqn {
   public:
    using State = std::array<double, 5>;

    ForwardShockEqn(Medium const& medium, Ejecta const& jet, Ejecta const& inject, double phi, double theta,
                    double eps_e);

    Medium const& medium;
    Ejecta const& jet;
    Ejecta const& inject;
    double const phi{0};
    double const theta{0};
    double const eps_e{0};
    double const jet_sigma{0};
    double gamma4{1};
    double spreading_factor{1};

    void operator()(State const& y, State& dydr, double r);

   private:
    inline double dGammadr(double r, double Gamma, double u, double t_eng, double ad_idx, double rho, double dtdr);
    inline double dUdr(double r, double Gamma, double u, double t_eng, double ad_idx, double rho, double dGdr);

    double const jet_Gamma0{0};
    double const inj_Gamma0{0};
    double const inj_sigma{0};
    double const dM0{0};
};

class FRShockEqn {
   public:
    using State = std::array<double, 5>;

    FRShockEqn(Medium const& medium, Ejecta const& jet, Ejecta const& inject, double phi, double theta);

    Medium const& medium;
    Ejecta const& jet;
    Ejecta const& inject;
    double const phi{0};
    double const theta{0};
    double const jet_sigma{0};
    double gamma4{1};

    void operator()(State const& y, State& dydr, double r);

   private:
    inline double dN3drPerOmega(double r, double n1, double n4, double gamma3);
};

Shock genForwardShock(Coord const& coord, Ejecta const& jet, Medium const& medium, double eps_e, double eps_B,
                      Ejecta const& inject = inject::none);

Shock genForwardShock3D(Coord const& coord, Ejecta const& jet, Medium const& medium, double eps_e, double eps_B,
                        Ejecta const& inject = inject::none);

std::pair<Shock, Shock> genFRShocks(Coord const& coord, Ejecta const& jet, Medium const& medium, double eps_e,
                                    double eps_B, Ejecta const& inject = inject::none);

std::pair<Shock, Shock> genFRShocks3D(Coord const& coord, Ejecta const& jet, Medium const& medium, double eps_e,
                                      double eps_B, Ejecta const& inject = inject::none);

double find_r_max(ForwardShockEqn& eqn, double r_min, double t_max);
#endif