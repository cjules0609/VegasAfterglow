#ifndef _SHOCK_
#define _SHOCK_
#include "mesh.h"
class Shock {
   public:
    Shock(Coord const& coord);
    MeshGrid t_com;
    MeshGrid Gamma;
    MeshGrid B;
};

inline double co_moving_shock_width(double r, double Gamma) { return r / Gamma / 12; }
#endif