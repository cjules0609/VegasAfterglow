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

    double pel{2.3};

    double I_nu(double nu) const;

   private:
    inline double I_nu_(double nu) const;
};
using ICPhotonArray = std::vector<ICPhoton>;
using ICPhotonMesh = std::vector<std::vector<ICPhoton>>;

double compton_sigma(double gamma, double nu);

ICPhotonMesh create_IC_photon_grid(size_t theta_size, size_t r_size,
                                   ICPhoton val = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.3});

ICPhotonMesh gen_IC_photons(Coord const& coord, Shock const& shock, SynElectronMesh const& electron,
                            SynElectronMesh const& photon, Medium const& medium);

MeshGrid IC_cooling_Thomson(Shock const& shock, SynElectronMesh& electron, ICPhotonMesh const& photon,
                            Medium const& medium);

MeshGrid IC_cooling_KN(Shock const& shock, SynElectronMesh& electron, ICPhotonMesh const& photon, Medium const& medium);
#endif