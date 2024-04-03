#ifndef _INVERSECOMPTON_
#define _INVERSECOMPTON_
#include <vector>

#include "medium.h"
#include "mesh.h"
#include "synchrotron.h"
struct ICRad {
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

    double pel{2.3};

    double I_nu(double nu) const;

   private:
    inline double I_nu_(double nu) const;
};
using ICRadArray = std::vector<ICRad>;
using ICRadMesh = std::vector<std::vector<ICRad>>;

double compton_sigma(double gamma, double nu);

ICRadMesh createICRadGrid(size_t theta_size, size_t r_size, ICRad val = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.3});

ICRadMesh calc_IC_radiation(Coord const& coord, MeshGrid const& Gamma, MeshGrid B, SynRadMesh const& electron,
                            SynRadMesh const& photon, Medium const& medium);

MeshGrid IC_cooling_noKN(MeshGrid const& Gamma, MeshGrid const& t_com, MeshGrid const& B, SynRadMesh& electron,
                    ICRadMesh const& photon, Medium const& medium, size_t order = 1);

MeshGrid IC_cooling(MeshGrid const& Gamma, MeshGrid const& t_com, MeshGrid const& B, SynRadMesh& electron,
                       ICRadMesh const& photon, Medium const& medium, size_t order = 1);
#endif