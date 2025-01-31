#ifndef _PROMPT_
#define _PROMPT_

#include "jet.h"
#include "mesh.h"
#include "physics.h"

struct PromptPhotons {
    double E_nu_peak{0};
    double nu_0{0};
    double alpha{0};

    double I_nu(double nu) const;
};

using PromptPhotonsArray = std::vector<PromptPhotons>;
using PromptPhotonsMesh = std::vector<std::vector<PromptPhotons>>;

PromptPhotonsMesh genPromptPhotons(Coord const& coord, Jet const& jet, double R0, double nu_0, double alpha);

#endif