#ifndef _SHOCK_
#define _SHOCK_
#include "mesh.h"
class Shock {
   public:
    Shock(Coord const& coord);
    MeshGrid t_com;
    MeshGrid Gamma;
    MeshGrid B;
    MeshGrid D_com;  // co-moving shock width
};
#endif