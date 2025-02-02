#ifndef _FSDYNAMICS_
#define _FSDYNAMICS_

#include <tuple>

#include "jet.h"
#include "medium.h"
#include "mesh.h"
#include "physics.h"
#include "shock.h"
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
    double phi{0};
    double theta{0};
    double eps_e{0};
    double jet_sigma0{0};
    double inject_sigma0{0};
    double gamma4{1};
    double spreading_factor{1};

    void operator()(State const& y, State& dydr, double r);

   private:
    inline double dGammadr(double r, double Gamma, double u, double t_eng);
    inline double dUdr(double r, double Gamma, double u, double t_eng);

    double jet_Gamma0{0};
    double inject_Gamma0{0};
};

class FRShockEqn {
   public:
    using State = std::array<double, 5>;

    FRShockEqn(Medium const& medium, Ejecta const& jet, Ejecta const& inject, double phi, double theta);

    Medium const& medium;
    Ejecta const& jet;
    Ejecta const& inject;
    double phi{0};
    double theta{0};
    double sigma{0};
    double gamma4{1};

    void operator()(State const& y, State& dydr, double r);

   private:
    inline double dN3drPerOmega(double r, double n1, double n4, double gamma3);
};

Shock genForwardShock(Coord const& coord, Ejecta const& jet, Medium const& medium, double eps_e, double eps_B,
                      Ejecta const& inject = Ejecta());

Shock genForwardShock3D(Coord const& coord, Ejecta const& jet, Medium const& medium, double eps_e, double eps_B,
                        Ejecta const& inject = Ejecta());

std::pair<Shock, Shock> genFRShocks(Coord const& coord, Ejecta const& jet, Medium const& medium, double eps_e,
                                    double eps_B, Ejecta const& inject = Ejecta());

std::pair<Shock, Shock> genFRShocks3D(Coord const& coord, Ejecta const& jet, Medium const& medium, double eps_e,
                                      double eps_B, Ejecta const& inject = Ejecta());
#endif