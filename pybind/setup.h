#pragma once

#include "afterglow.h"
#include "macros.h"
#include "mesh.h"
#include "utilities.h"

struct Params {
    double E_iso{1e52};
    double Gamma0{300};
    double theta_c{0.1};
    double theta_v{0};
    double theta_w{con::pi / 2};
    double p{2.3};
    double eps_e{0.1};
    double eps_B{0.01};
    double n_ism{1};
    double A_star{0.01};
    double xi{1};
    double k_jet{2};
};

struct ConfigParams {
    double lumi_dist{1e26};
    double z{0};
    std::string medium{"ism"};
    std::string jet{"tophat"};
    size_t t_grid{24};
    size_t phi_grid{24};
    size_t theta_grid{24};
    double rtol{1e-5};
};