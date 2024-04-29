#ifndef _CONSTS_
#define _CONSTS_

namespace con {

    constexpr double sec = 1.0 / 500;
    constexpr double cm = 1 / 1.5e13;
    constexpr double g = 1 / 2e33;
    constexpr double Hz = 1 / sec;
    constexpr double erg = g * cm * cm / sec / sec;
    constexpr double c = 1;
    constexpr double c2 = c * c;
    constexpr double mp = 1.67e-24 * g;
    constexpr double me = mp / 1836;
    constexpr double h = 6.63e-27 * erg * sec;
    // constexpr double kB = 1.38e-16;
    constexpr double e = 4.8e-10 / 4.472136e16 / 5.809475e19 / sec;  // M^(1/2)L^(3/2)/T
    constexpr double eV = 1.60218e-12 * erg;
    constexpr double keV = 1e3 * eV;
    constexpr double MeV = 1e6 * eV;
    constexpr double GeV = 1e9 * eV;
    constexpr double TeV = 1e12 * eV;
    constexpr double e2 = e * e;
    constexpr double e3 = e2 * e;
    constexpr double pi = 3.14159265358979323846;
    constexpr double sigmaT = 6.65e-25 * cm * cm;
    constexpr double deg = pi / 180;

    constexpr double au = 1;
    constexpr double m = 100 * cm;
    constexpr double km = 1000 * m;
    constexpr double pc = 2.06265e5 * au;
    constexpr double kpc = 1000 * pc;
    constexpr double Mpc = 1e6 * pc;
    constexpr double hr = 3600 * sec;
    constexpr double day = 24 * hr;
    constexpr double yr = 365.2425 * day;
    constexpr double Myr = 1e6 * yr;
    constexpr double Gyr = 1e9 * yr;
    constexpr double Omega_m = 0.27;
    constexpr double Omega_L = 0.73;
    constexpr double H0 = 67.66 * km / sec / Mpc;
}  // namespace con

#endif