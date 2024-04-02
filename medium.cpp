#include "medium.h"

#include "macros.h"

Medium create_ISM(double n_ism, double eps_e, double eps_B, double pel) {
    auto rho = [=](double r) { return n_ism * con::mp; };
    auto m = [=](double r) { return 4 * con::pi / 3 * r * r * r * n_ism * con::mp; };
    return Medium{rho, m, eps_e, eps_B, pel};
}