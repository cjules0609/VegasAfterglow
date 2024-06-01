#include "medium.h"

#include <cmath>

#include "macros.h"

Medium create_ISM(double n_ism, double cs) {
    auto rho = [=](double r) { return n_ism * con::mp; };
    auto m = [=](double r) { return 4 * con::pi / 3 * r * r * r * n_ism * con::mp; };
    return Medium{rho, m, cs};
}