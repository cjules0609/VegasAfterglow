#ifndef _RELATIVITY_H_
#define _RELATIVITY_H_
#include <cmath>

#include "macros.h"

inline double gamma_to_beta(double gamma) { return sqrt(1 - 1 / (gamma * gamma)); }

inline double adiabatic_index(double gamma) { return (4 * gamma + 1) / (3 * gamma); }

#endif  // _RELATIVITY_H_
