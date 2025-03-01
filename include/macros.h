//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _CONSTS_
#define _CONSTS_
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
namespace con {
    constexpr double len = 1.5e13;
    // Length unit: 1 cm is defined as 1 / len (arbitrary unit conversion)
    constexpr double cm = 1 / len;

    // Time unit: 1 second is defined as c/len (arbitrary unit conversion)
    constexpr double sec = 3e10 / len;

    constexpr double cm2 = cm * cm;       // Square centimeters
    constexpr double cm3 = cm * cm * cm;  // Cubic centimeters
    // Mass unit: 1 gram is defined as 1 / 2e33 (arbitrary unit conversion)
    constexpr double g = 1 / 2e33;
    // Frequency unit: 1 Hz is defined as 1 / sec
    constexpr double Hz = 1 / sec;
    // Energy unit: 1 erg, derived from g, cm, and sec (erg = g * cm^2 / sec^2)
    constexpr double erg = g * cm * cm / sec / sec;
    // Speed of light and its square (set to 1 in these natural units)
    constexpr double c = 1;
    constexpr double c2 = c * c;
    // Proton mass (in grams) using the defined g
    constexpr double mp = 1.67e-24 * g;
    // Electron mass (set relative to the proton mass)
    constexpr double me = mp / 1836;
    // Planck's constant (in erg*sec)
    constexpr double h = 6.63e-27 * erg * sec;
    // Boltzmann constant is commented out (could be defined if needed)
    // constexpr double kB = 1.38e-16;
    // Elementary charge (e) with a complex conversion factor.
    constexpr double e = 4.8e-10 / 4.472136e16 / 5.809475e19 / sec;  // Units: M^(1/2)L^(3/2)/T
    // Energy units in electronvolts
    constexpr double eV = 1.60218e-12 * erg;
    constexpr double keV = 1e3 * eV;
    constexpr double MeV = 1e6 * eV;
    constexpr double GeV = 1e9 * eV;
    constexpr double TeV = 1e12 * eV;
    // Powers of the elementary charge.
    constexpr double e2 = e * e;
    constexpr double e3 = e2 * e;
    // Mathematical constant π.
    constexpr double pi = 3.14159265358979323846;
    // Thomson cross-section in cm².
    constexpr double sigmaT = 6.65e-25 * cm * cm;
    // Conversion factor from degrees to radians.
    constexpr double deg = pi / 180;

    // Astronomical unit (AU), set to 1 in these units.
    constexpr double au = 1;
    // Meter, defined as 100 centimeters.
    constexpr double m = 100 * cm;
    // Kilometer, defined as 1000 meters.
    constexpr double km = 1000 * m;
    // Parsec (pc) in AU; note: conversion factor provided.
    constexpr double pc = 2.06265e5 * au;
    // Kiloparsec (kpc) and Megaparsec (Mpc).
    constexpr double kpc = 1000 * pc;
    constexpr double Mpc = 1e6 * pc;
    // Time units: hour, day, year, and their multiples.
    constexpr double hr = 3600 * sec;
    constexpr double day = 24 * hr;
    constexpr double yr = 365.2425 * day;
    constexpr double Myr = 1e6 * yr;
    constexpr double Gyr = 1e9 * yr;
    // Cosmological density parameters.
    constexpr double Omega_m = 0.27;
    constexpr double Omega_L = 0.73;
    // Hubble constant H0 in km/s/Mpc.
    constexpr double H0 = 67.66 * km / sec / Mpc;

    // Gamma cut-off, defined very close to 1.
    constexpr double Gamma_cut = 1 + 1e-6;
    constexpr Real inf = std::numeric_limits<Real>::infinity();
}  // namespace con
#endif
