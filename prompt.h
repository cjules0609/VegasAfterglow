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

using PromptPhotonsGrid = boost::multi_array<PromptPhotons, 3>;
PromptPhotonsGrid createPromptPhotonsGrid(size_t phi_size, size_t theta_size, size_t r_size);
PromptPhotonsGrid genPromptPhotons(Coord const& coord, Ejecta const& jet, double R0, double nu_0, double alpha);

#endif