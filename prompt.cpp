#include "prompt.h"

PromptPhotonsGrid createPromptPhotonsGrid(size_t phi_size, size_t theta_size, size_t r_size) {
    return PromptPhotonsGrid(boost::extents[phi_size][theta_size][r_size]);
}

double PromptPhotons::I_nu(double nu) const { return E_nu_peak * std::pow(nu / nu_0, -alpha); }