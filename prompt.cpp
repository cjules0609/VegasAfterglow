//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "prompt.h"

PromptPhotonsGrid createPromptPhotonsGrid(size_t phi_size, size_t theta_size, size_t r_size) {
    return PromptPhotonsGrid(boost::extents[phi_size][theta_size][r_size]);
}

double PromptPhotons::I_nu(double nu) const { return E_nu_peak * std::pow(nu / nu_0, -alpha); }