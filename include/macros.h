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

namespace con {
    constexpr Real len = 1.5e13;
    // Length unit: 1 cm is defined as 1 / len (arbitrary unit conversion)
    constexpr Real cm = 1 / len;

    // Time unit: 1 second is defined as c/len (arbitrary unit conversion)
    constexpr Real sec = 3e10 / len;

    constexpr Real cm2 = cm * cm;       // Square centimeters
    constexpr Real cm3 = cm * cm * cm;  // Cubic centimeters
    // Mass unit: 1 gram is defined as 1 / 2e33 (arbitrary unit conversion)
    constexpr Real g = 1 / 2e33;

    constexpr Real Gauss = 8.66e-11 / sec;  // sqrt(g/cm)/sec
    // Frequency unit: 1 Hz is defined as 1 / sec
    constexpr Real Hz = 1 / sec;
    // Energy unit: 1 erg, derived from g, cm, and sec (erg = g * cm^2 / sec^2)
    constexpr Real erg = g * cm * cm / sec / sec;
    // Speed of light and its square (set to 1 in these natural units)
    constexpr Real c = 1;
    constexpr Real c2 = c * c;
    // Proton mass (in grams) using the defined g
    constexpr Real mp = 1.67e-24 * g;
    // Electron mass (set relative to the proton mass)
    constexpr Real me = mp / 1836;
    constexpr Real M_sun = 2e33 * g;
    // Planck's constant (in erg*sec)
    constexpr Real h = 6.63e-27 * erg * sec;
    // Boltzmann constant is commented out (could be defined if needed)
    // constexpr double kB = 1.38e-16;
    // Elementary charge (e) with a complex conversion factor.
    constexpr Real e = 4.8e-10 / 4.472136e16 / 5.809475e19 / sec;  // Units: M^(1/2)L^(3/2)/T
    // Energy units in electronvolts
    constexpr Real eV = 1.60218e-12 * erg;
    constexpr Real keV = 1e3 * eV;
    constexpr Real MeV = 1e6 * eV;
    constexpr Real GeV = 1e9 * eV;
    constexpr Real TeV = 1e12 * eV;
    // Powers of the elementary charge.
    constexpr Real e2 = e * e;
    constexpr Real e3 = e2 * e;
    // Mathematical constant π.
    constexpr Real pi = 3.14159265358979323846;
    // Thomson cross-section in cm².
    constexpr Real sigmaT = 6.65e-25 * cm * cm;
    // Conversion factor from degrees to radians.
    constexpr Real deg = pi / 180;

    constexpr Real Jy = 1e-23 * erg / cm2 / sec / Hz;

    constexpr Real mJy = 1e-3 * Jy;

    constexpr Real uJy = 1e-6 * mJy;

    // Astronomical unit (AU), set to 1 in these units.
    constexpr Real au = 1;
    // Meter, defined as 100 centimeters.
    constexpr Real m = 100 * cm;
    // Kilometer, defined as 1000 meters.
    constexpr Real km = 1000 * m;
    // Parsec (pc) in AU; note: conversion factor provided.
    constexpr Real pc = 2.06265e5 * au;
    // Kiloparsec (kpc) and Megaparsec (Mpc).
    constexpr Real kpc = 1000 * pc;
    constexpr Real Mpc = 1e6 * pc;
    // Time units: hour, day, year, and their multiples.
    constexpr Real hr = 3600 * sec;
    constexpr Real day = 24 * hr;
    constexpr Real yr = 365.2425 * day;
    constexpr Real Myr = 1e6 * yr;
    constexpr Real Gyr = 1e9 * yr;
    // Cosmological density parameters.
    constexpr Real Omega_m = 0.27;
    constexpr Real Omega_L = 0.73;
    // Hubble constant H0 in km/s/Mpc.
    constexpr Real H0 = 67.66 * km / sec / Mpc;

    // Gamma cut-off, defined very close to 1.
    constexpr Real Gamma_cut = 1 + 1e-6;
    constexpr Real inf = std::numeric_limits<Real>::infinity();
}  // namespace con
