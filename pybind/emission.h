#pragma once

#include "afterglow.h"
#include "macros.h"
#include "mesh.h"
#include "utilities.h"
#include "setup.h"

class Emission {
public:
    Emission(Real nu, const std::vector<Real>& t);
    Emission(const std::vector<Real>& nu, Real t);

    void generate(const Params& param, const ConfigParams& config);

    // Vectors
    Array t;
    Array nu;

    // Result
    Array F_nu;
};

