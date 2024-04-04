#ifndef _CONSTS_
#define _CONSTS_

namespace con {
    constexpr double Hz = 500;
    constexpr double sec = 1.0 / 500;
    constexpr double cm = 1 / 1.5e13;
    constexpr double g = 1 / 2e33;
    constexpr double erg = 1 / 9e20 / 2e33;
    constexpr double c = 1;
    constexpr double c2 = c * c;
    constexpr double mp = 1.67e-24 * g;
    constexpr double me = mp / 1836;
    constexpr double h = 6.63e-27 * erg * sec;
    // constexpr double kB = 1.38e-16;
    constexpr double e = 4.8e-10 / 4.472136e16 / 5.809475e19 / sec;  // M^(1/2)L^(3/2)/T
    constexpr double e2 = e * e;
    constexpr double e3 = e2 * e;
    constexpr double pi = 3.14159265358979323846;
    constexpr double sigmaT = 6.65e-25 * cm * cm;
    constexpr double deg = pi / 180;
}  // namespace con

#endif