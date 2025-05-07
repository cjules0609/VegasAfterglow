#pragma once

#include "afterglow.h"
#include "macros.h"
#include "mesh.h"
#include "utilities.h"
#include "struct.h"

class Emission {
public:
    Emission(const Params& param, const ConfigParams& config);

    std::vector<Real> generate_lc(Real nu, const std::vector<Real>& t);
    std::vector<Real> generate_spec(const std::vector<Real>& nu, Real t);

    Array t;
    Array nu;
    Array F_nu;

    Ejecta jet;
    Medium medium;
    Observer obs;

private:
    Params param;
    ConfigParams config;

    void observe(const Params& param, const ConfigParams& config, const Array& t);

    SynElectronGrid electrons;
    SynPhotonGrid photons;
};
