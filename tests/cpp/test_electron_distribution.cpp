#include <boost/test/unit_test.hpp>

#include <cmath>
#include <limits>

#include "radiation/electron-distribution.h"
#include "util/macros.h"

BOOST_AUTO_TEST_SUITE(ElectronDistribution)

BOOST_AUTO_TEST_CASE(built_in_shapes) {
    FlatElectronShape flat(10, 100);
    BOOST_CHECK_EQUAL(flat.evaluate(9), 0);
    BOOST_CHECK_EQUAL(flat.evaluate(10), 1);
    BOOST_CHECK_EQUAL(flat.evaluate(50), 1);
    BOOST_CHECK_EQUAL(flat.evaluate(101), 0);

    PowerLawElectronShape power_law(10, 100, 2.5);
    BOOST_CHECK_CLOSE(power_law.evaluate(20), std::pow(20.0, -2.5), 1e-10);
    BOOST_CHECK_EQUAL(power_law.evaluate(1), 0);

    CutoffPowerLawElectronShape cutoff(10, 100, 2.5, 40);
    BOOST_CHECK_CLOSE(cutoff.evaluate(20), std::pow(20.0, -2.5) * std::exp(-0.5), 1e-10);
    BOOST_CHECK_LT(cutoff.evaluate(80), power_law.evaluate(80));
}

BOOST_AUTO_TEST_CASE(parameter_validation) {
    BOOST_CHECK_THROW(FlatElectronShape(0.99, 100), std::invalid_argument);
    BOOST_CHECK_THROW(FlatElectronShape(10, 10), std::invalid_argument);
    BOOST_CHECK_THROW(FlatElectronShape(10, con::inf), std::invalid_argument);
    BOOST_CHECK_THROW(PowerLawElectronShape(10, 100, con::inf), std::invalid_argument);
    BOOST_CHECK_THROW(CutoffPowerLawElectronShape(10, 100, 2.3, 0), std::invalid_argument);
    BOOST_CHECK_THROW(CutoffPowerLawElectronShape(10, 100, 2.3, con::inf), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(flat_integrals) {
    constexpr Real a = 1;
    constexpr Real b = 1e5;
    auto sampled = sample_electron_shape(FlatElectronShape(a, b), 32);
    const Real number_exact = b - a;
    const Real energy_exact = 0.5 * (b * b - a * a) - (b - a);
    BOOST_CHECK_SMALL(sampled->number_integral / number_exact - 1, 1e-10);
    BOOST_CHECK_SMALL(sampled->energy_integral / energy_exact - 1, 3e-10);
}

BOOST_AUTO_TEST_CASE(power_law_integrals) {
    constexpr Real a = 10;
    constexpr Real b = 1e7;
    constexpr Real p = 2.3;
    auto sampled = sample_electron_shape(PowerLawElectronShape(a, b, p), 32);
    const Real number_exact = (std::pow(b, 1 - p) - std::pow(a, 1 - p)) / (1 - p);
    const Real gamma_moment = (std::pow(b, 2 - p) - std::pow(a, 2 - p)) / (2 - p);
    const Real energy_exact = gamma_moment - number_exact;
    BOOST_CHECK_SMALL(sampled->number_integral / number_exact - 1, 2e-8);
    BOOST_CHECK_SMALL(sampled->energy_integral / energy_exact - 1, 2e-8);
}

BOOST_AUTO_TEST_CASE(cutoff_integral_convergence) {
    NumericalElectronShape shape = CutoffPowerLawElectronShape(10, 1e8, 2.3, 1e5);
    auto coarse = sample_electron_shape(shape, 8);
    auto medium = sample_electron_shape(shape, 16);
    auto fine = sample_electron_shape(shape, 32);
    const Real coarse_error = std::abs(coarse->energy_integral / fine->energy_integral - 1);
    const Real medium_error = std::abs(medium->energy_integral / fine->energy_integral - 1);
    BOOST_CHECK_LT(medium_error, coarse_error);
    BOOST_CHECK_LT(medium_error, 2e-5);
}

BOOST_AUTO_TEST_CASE(pre_sampled_native_shape_uses_the_builtin_table_path) {
    Array gamma = make_electron_gamma_grid(10, 1e6, 32);
    Array values = xt::ones<Real>({gamma.size()});
    auto sampled = sample_electron_shape(gamma, values);
    auto built_in = sample_electron_shape(FlatElectronShape(10, 1e6), 32);

    BOOST_CHECK_EQUAL_COLLECTIONS(sampled->gamma.begin(), sampled->gamma.end(), built_in->gamma.begin(),
                                  built_in->gamma.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(sampled->shape.begin(), sampled->shape.end(), built_in->shape.begin(),
                                  built_in->shape.end());
    BOOST_CHECK_EQUAL(sampled->number_integral, built_in->number_integral);
    BOOST_CHECK_EQUAL(sampled->energy_integral, built_in->energy_integral);

    Array wrong_size = xt::ones<Real>({gamma.size() - 1});
    BOOST_CHECK_THROW((void)sample_electron_shape(gamma, wrong_size), std::invalid_argument);

    Array invalid = values;
    invalid(3) = -1;
    BOOST_CHECK_THROW((void)sample_electron_shape(gamma, invalid), std::invalid_argument);
    invalid(3) = std::numeric_limits<Real>::quiet_NaN();
    BOOST_CHECK_THROW((void)sample_electron_shape(gamma, invalid), std::invalid_argument);

    Array zero = xt::zeros<Real>({gamma.size()});
    BOOST_CHECK_THROW((void)sample_electron_shape(gamma, zero), std::invalid_argument);

    Array nonuniform = gamma;
    nonuniform(3) *= 1.01;
    BOOST_CHECK_THROW((void)sample_electron_shape(nonuniform, values), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(energy_normalization) {
    auto shape = sample_electron_shape(FlatElectronShape(10, 1000));
    constexpr Real proton_column = 2e20;
    constexpr Real Gamma_th = 3;
    constexpr Real eps_e = 0.1;
    NumericalElectronDistribution electrons(shape, ElectronNormalization::energy, proton_column, Gamma_th, eps_e,
                                            0.3);

    BOOST_CHECK_SMALL(electrons.implied_eps_e / eps_e - 1, 1e-12);
    BOOST_CHECK_SMALL(electrons.number_column / proton_column - electrons.implied_xi_e, 1e-12);
    BOOST_CHECK_GT(electrons.compute_column_den(100), 0);
    BOOST_CHECK_EQUAL(electrons.compute_column_den(9), 0);
    BOOST_CHECK_EQUAL(electrons.compute_column_den(1001), 0);
}

BOOST_AUTO_TEST_CASE(number_normalization_and_unclipped_fractions) {
    auto shape = sample_electron_shape(FlatElectronShape(1e4, 1e6));
    constexpr Real proton_column = 1e20;
    constexpr Real xi_e = 2.0;
    NumericalElectronDistribution electrons(shape, ElectronNormalization::number, proton_column, 1.01, 0.1, xi_e);

    BOOST_CHECK_SMALL(electrons.number_column / (xi_e * proton_column) - 1, 1e-12);
    BOOST_CHECK_SMALL(electrons.implied_xi_e / xi_e - 1, 1e-12);
    BOOST_CHECK_GT(electrons.implied_eps_e, 1);
}

BOOST_AUTO_TEST_CASE(distribution_validation) {
    auto shape = sample_electron_shape(FlatElectronShape(10, 100));
    BOOST_CHECK_THROW(NumericalElectronDistribution(nullptr, ElectronNormalization::energy, 1, 2, 0.1, 1),
                      std::invalid_argument);
    BOOST_CHECK_THROW(NumericalElectronDistribution(shape, ElectronNormalization::energy, -1, 2, 0.1, 1),
                      std::invalid_argument);
    BOOST_CHECK_THROW(NumericalElectronDistribution(shape, ElectronNormalization::energy, 1, 0.9, 0.1, 1),
                      std::invalid_argument);
    BOOST_CHECK_THROW(NumericalElectronDistribution(shape, ElectronNormalization::energy, 1, 2, con::inf, 1),
                      std::invalid_argument);
    BOOST_CHECK_THROW(NumericalElectronDistribution(shape, ElectronNormalization::number, 1, 2, 0.1, con::inf),
                      std::invalid_argument);
    BOOST_CHECK_THROW(NumericalElectronDistribution(shape, static_cast<ElectronNormalization>(99), 1, 2, 0.1, 1),
                      std::invalid_argument);
    BOOST_CHECK_THROW((void)sample_electron_shape(FlatElectronShape(10, 100), 0), std::invalid_argument);

    NumericalElectronShape overflowing = PowerLawElectronShape(1, 1e200, -2);
    BOOST_CHECK_THROW((void)sample_electron_shape(overflowing), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(number_normalization_at_zero_thermal_energy_is_visible) {
    auto shape = sample_electron_shape(FlatElectronShape(10, 100));
    NumericalElectronDistribution electrons(shape, ElectronNormalization::number, 1e20, 1, 0.1, 1);
    BOOST_CHECK(std::isinf(electrons.implied_eps_e));
}

BOOST_AUTO_TEST_SUITE_END()
