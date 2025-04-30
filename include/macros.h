//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <limits>
#include <numeric>
/********************************************************************************************************************
 * NAMESPACE: con
 * DESCRIPTION: Contains constant definitions used throughout the simulation. These constants define unit
 *              conversions (e.g., seconds, centimeters, grams), physical constants (e.g., speed of light, masses,
 *              Planck's constant), energy units, and cosmological parameters.
 ********************************************************************************************************************/
using Real = double;
constexpr Real operator"" _r(long double x) { return static_cast<Real>(x); }

// Alternatively, if Real might be float and you want to support integer literals too:
constexpr Real operator"" _r(unsigned long long x) { return static_cast<Real>(x); }

namespace unit {
    constexpr Real len = 1.5e13;
    constexpr Real cm = 1 / len;
    constexpr Real sec = 3e10 / len;
    constexpr Real cm2 = cm * cm;       // Square centimeters
    constexpr Real cm3 = cm * cm * cm;  // Cubic centimeters
    constexpr Real g = 1 / 2e33;
    constexpr Real kg = 1000 * g;
    constexpr Real Gauss = 8.66e-11 / sec;  // sqrt(g/cm)/sec
    constexpr Real Hz = 1 / sec;
    constexpr Real erg = g * cm * cm / sec / sec;
    constexpr Real M_sun = 2e33 * g;
    constexpr Real eV = 1.60218e-12 * erg;
    constexpr Real keV = 1e3 * eV;
    constexpr Real MeV = 1e6 * eV;
    constexpr Real GeV = 1e9 * eV;
    constexpr Real TeV = 1e12 * eV;
    constexpr Real deg = 3.14159265358979323846 / 180;
    constexpr Real flux_cgs = erg / cm2 / sec;
    constexpr Real flux_den_cgs = erg / cm2 / sec / Hz;
    constexpr Real Jy = 1e-23 * erg / cm2 / sec / Hz;
    constexpr Real mJy = 1e-3 * Jy;
    constexpr Real uJy = 1e-6 * mJy;
    constexpr Real m = 100 * cm;
    constexpr Real km = 1000 * m;
    constexpr Real au = 1.5e13 * cm;
    constexpr Real pc = 2.06265e5 * au;
    constexpr Real kpc = 1000 * pc;
    constexpr Real Mpc = 1e6 * pc;
    constexpr Real hr = 3600 * sec;
    constexpr Real day = 24 * hr;
    constexpr Real yr = 365.2425 * day;
    constexpr Real Myr = 1e6 * yr;
    constexpr Real Gyr = 1e9 * yr;
}  // namespace unit
namespace con {
    constexpr Real c = 1;
    constexpr Real c2 = c * c;
    constexpr Real mp = 1.67e-24 * unit::g;
    constexpr Real me = mp / 1836;
    constexpr Real h = 6.63e-27 * unit::erg * unit::sec;
    constexpr Real e = 4.8e-10 / 4.472136e16 / 5.809475e19 / unit::sec;  // Units: M^(1/2)L^(3/2)/T
    constexpr Real e2 = e * e;
    constexpr Real e3 = e2 * e;
    constexpr Real pi = 3.14159265358979323846;
    constexpr Real sigmaT = 6.65e-25 * unit::cm * unit::cm;  // Thomson cross-section in cm².
    constexpr Real Omega_m = 0.27;
    constexpr Real Omega_L = 0.73;
    constexpr Real H0 = 67.66 * unit::km / unit::sec / unit::Mpc;
    constexpr Real Gamma_cut = 1 + 1e-5;
    constexpr Real inf = std::numeric_limits<Real>::infinity();
}  // namespace con
