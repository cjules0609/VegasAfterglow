#include <boost/test/unit_test.hpp>

#include <array>
#include <cmath>
#include <limits>

#include "radiation/synchrotron-kernel.h"

BOOST_AUTO_TEST_SUITE(SynchrotronKernel)

BOOST_AUTO_TEST_CASE(reference_values) {
    constexpr std::array<std::pair<Real, Real>, 16> reference = {{
        {1e-8, 0.004631000072771483},
        {1e-6, 0.021493468615984595},
        {1e-4, 0.09959088308506682},
        {0.01, 0.44497250411421074},
        {0.05, 0.70157192969349103},
        {0.1, 0.81818553487285406},
        {0.2, 0.90338599391440655},
        {0.285812, 0.9180123331838651},
        {0.5, 0.8708191468754688},
        {1, 0.65142281535536328},
        {2, 0.30163590285073949},
        {5, 0.021248129774981972},
        {10, 0.00019223826430086864},
        {20, 1.1968634456097954e-8},
        {50, 1.7347852040130392e-21},
        {100, 4.697593665922195e-43},
    }};

    for (const auto& [x, expected] : reference) {
        const Real actual = synchrotron_kernel(x);
        BOOST_CHECK_SMALL(actual / expected - 1, 1e-4);
        BOOST_CHECK_SMALL(std::exp2(log2_synchrotron_kernel(std::log2(x))) / expected - 1, 1e-4);
    }
}

BOOST_AUTO_TEST_CASE(absorption_kernel_reference_values) {
    constexpr std::array<std::pair<Real, Real>, 16> reference = {{
        {1e-8, 0.003087345473843415},
        {1e-6, 0.014330188276891179},
        {1e-4, 0.066514841640181899},
        {0.01, 0.30872298702966244},
        {0.05, 0.52743969618307729},
        {0.1, 0.66272663682545552},
        {0.2, 0.82631773600156555},
        {0.285812, 0.9180121098184425},
        {0.5, 1.051211411708572},
        {1, 1.0977307162471457},
        {2, 0.7990836518196639},
        {5, 0.11885137496465183},
        {10, 0.0020296099946992691},
        {20, 2.4575004213034363e-7},
        {50, 8.7631516887388006e-20},
        {100, 4.7214289764710537e-41},
    }};

    for (const auto& [x, expected] : reference) {
        BOOST_CHECK_SMALL(synchrotron_absorption_kernel(x) / expected - 1, 1e-4);
        BOOST_CHECK_SMALL(std::exp2(log2_synchrotron_absorption_kernel(std::log2(x))) / expected - 1, 1e-4);
    }
}

BOOST_AUTO_TEST_CASE(low_x_asymptote) {
    constexpr Real coefficient = 2.1495282415344786;
    const Real x1 = 1e-12;
    const Real x2 = 1e-9;
    BOOST_CHECK_SMALL(synchrotron_kernel(x1) / (coefficient * std::cbrt(x1)) - 1, 1e-7);
    const Real slope = std::log(synchrotron_kernel(x2) / synchrotron_kernel(x1)) / std::log(x2 / x1);
    BOOST_CHECK_SMALL(slope - 1.0 / 3.0, 1e-5);

    constexpr Real absorption_coefficient = 1.4330188276896527;
    BOOST_CHECK_SMALL(synchrotron_absorption_kernel(x1) /
                              (absorption_coefficient * std::cbrt(x1)) -
                          1,
                      1e-10);
}

BOOST_AUTO_TEST_CASE(high_x_asymptote) {
    const Real x = 100;
    const Real leading = std::sqrt(con::pi * x / 2) * std::exp(-x);
    const Real corrected = leading * (1 + 55.0 / (72 * x) - 10151.0 / (10368 * x * x));
    BOOST_CHECK_SMALL(synchrotron_kernel(x) / corrected - 1, 1e-12);

    const Real absorption_corrected = std::sqrt(con::pi / 2) * std::pow(x, 1.5) * std::exp(-x) *
                                      (1 + 91.0 / (72 * x) + 1729.0 / (10368 * x * x));
    BOOST_CHECK_SMALL(synchrotron_absorption_kernel(x) / absorption_corrected - 1, 1e-8);
}

BOOST_AUTO_TEST_CASE(peak) {
    Real x_peak = 0;
    Real f_peak = 0;
    for (size_t i = 0; i <= 20000; ++i) {
        const Real log_x = std::log(0.1) + static_cast<Real>(i) / 20000 * std::log(10.0);
        const Real x = std::exp(log_x);
        const Real value = synchrotron_kernel(x);
        if (value > f_peak) {
            x_peak = x;
            f_peak = value;
        }
    }
    BOOST_CHECK_CLOSE(x_peak, 0.2858, 0.1);
    BOOST_CHECK_CLOSE(f_peak, 0.9180, 0.02);
}

BOOST_AUTO_TEST_CASE(asymptotic_boundaries_are_continuous) {
    for (Real boundary : {1e-2, 50.0}) {
        const Real below = synchrotron_kernel(boundary * (1 - 1e-10));
        const Real above = synchrotron_kernel(boundary * (1 + 1e-10));
        BOOST_CHECK_SMALL(above / below - 1, 3e-5);
        const Real absorption_below = synchrotron_absorption_kernel(boundary * (1 - 1e-10));
        const Real absorption_above = synchrotron_absorption_kernel(boundary * (1 + 1e-10));
        BOOST_CHECK_SMALL(absorption_above / absorption_below - 1, 3e-5);
    }
}

BOOST_AUTO_TEST_CASE(invalid_and_underflow_inputs) {
    BOOST_CHECK_EQUAL(synchrotron_kernel(0), 0);
    BOOST_CHECK_EQUAL(synchrotron_kernel(-1), 0);
    BOOST_CHECK_EQUAL(synchrotron_kernel(con::inf), 0);
    BOOST_CHECK_EQUAL(synchrotron_kernel(std::numeric_limits<Real>::quiet_NaN()), 0);
    BOOST_CHECK_EQUAL(synchrotron_kernel(1e4), 0);
    BOOST_CHECK_EQUAL(log2_synchrotron_kernel(-con::inf), -con::inf);
    BOOST_CHECK_EQUAL(synchrotron_absorption_kernel(0), 0);
    BOOST_CHECK_EQUAL(synchrotron_absorption_kernel(-1), 0);
    BOOST_CHECK_EQUAL(synchrotron_absorption_kernel(con::inf), 0);
    BOOST_CHECK_EQUAL(log2_synchrotron_absorption_kernel(-con::inf), -con::inf);
}

BOOST_AUTO_TEST_SUITE_END()
