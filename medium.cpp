//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "medium.h"

#include <cmath>

#include "macros.h"

Medium::Medium(double k, double n_c, double r_c) : k(k), n_c(n_c), r_c(r_c) {};
Medium::Medium(double n_c) : k(0), n_c(n_c), r_c(con::cm) {};
Medium createISM(double n_ism) { return Medium(n_ism); }
