#include "medium.h"

#include <cmath>

#include "macros.h"

Medium::Medium(double k, double n_c, double r_c) : k(k), n_c(n_c), r_c(r_c) {};
Medium::Medium(double n_c) : k(0), n_c(n_c), r_c(con::cm) {};


Medium createISM(double n_ism) { return Medium(n_ism); }