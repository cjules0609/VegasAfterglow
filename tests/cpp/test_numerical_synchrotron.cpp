#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>

#include "radiation/electron-distribution.h"
#include "radiation/numerical-synchrotron.h"
#include "radiation/synchrotron.h"
#include "util/macros.h"

namespace {
    template <typename T>
    concept HasStandardBreaks = requires(T value) {
        value.gamma_m;
        value.gamma_c;
        value.gamma_a;
    };

    Real peak_frequency(NumericalSynchrotron const& photons, Real nu_min, Real nu_max) {
        Real peak_nu = nu_min;
        Real peak_I = 0;
        for (size_t i = 0; i <= 4000; ++i) {
            const Real nu = std::exp(std::log(nu_min) + static_cast<Real>(i) / 4000 * std::log(nu_max / nu_min));
            const Real intensity = photons.compute_I_nu(nu);
            if (intensity > peak_I) {
                peak_I = intensity;
                peak_nu = nu;
            }
        }
        return peak_nu;
    }
} // namespace

static_assert(!HasStandardBreaks<NumericalElectronDistribution>);
static_assert(!HasStandardBreaks<NumericalSynchrotron>);

BOOST_AUTO_TEST_SUITE(NumericalSynchrotronTests)

BOOST_AUTO_TEST_CASE(column_and_magnetic_scaling) {
    auto shape = sample_electron_shape(PowerLawElectronShape(10, 1e6, 2.3));
    auto table = build_numerical_synchrotron_table(shape);
    NumericalElectronDistribution e1(shape, ElectronNormalization::number, 1e20, 2, 0.1, 0.2);
    NumericalElectronDistribution e2(shape, ElectronNormalization::number, 1e20, 2, 0.1, 0.6);
    const Real B1 = 0.1 * unit::Gauss;
    const Real B2 = 0.4 * unit::Gauss;
    NumericalSynchrotron p1(table, e1, B1);
    NumericalSynchrotron p2(table, e2, B1);
    NumericalSynchrotron pB(table, e1, B2);

    const Real u = 1e7;
    const Real nu1 = compute_syn_freq(1, B1) * u;
    const Real nu2 = compute_syn_freq(1, B2) * u;
    BOOST_CHECK_SMALL(p2.compute_I_nu(nu1) / p1.compute_I_nu(nu1) - 3, 2e-12);
    BOOST_CHECK_SMALL(pB.compute_I_nu(nu2) / p1.compute_I_nu(nu1) - B2 / B1, 2e-12);
}

BOOST_AUTO_TEST_CASE(characteristic_frequency_scales_as_B_gamma_squared) {
    const Real B1 = 0.1 * unit::Gauss;
    const Real B2 = 0.3 * unit::Gauss;
    auto shape1 = sample_electron_shape(FlatElectronShape(100, 101), 64);
    auto shape2 = sample_electron_shape(FlatElectronShape(200, 202), 64);
    auto table1 = build_numerical_synchrotron_table(shape1, 48);
    auto table2 = build_numerical_synchrotron_table(shape2, 48);
    NumericalElectronDistribution e1(shape1, ElectronNormalization::number, 1e20, 2, 0.1, 1);
    NumericalElectronDistribution e2(shape2, ElectronNormalization::number, 1e20, 2, 0.1, 1);
    NumericalSynchrotron p1(table1, e1, B1);
    NumericalSynchrotron p2(table2, e2, B2);

    const Real peak1 = peak_frequency(p1, compute_syn_freq(100, B1) * 0.03, compute_syn_freq(101, B1) * 3);
    const Real peak2 = peak_frequency(p2, compute_syn_freq(200, B2) * 0.03, compute_syn_freq(202, B2) * 3);
    BOOST_CHECK_CLOSE(peak2 / peak1, (B2 / B1) * 4, 0.3);
}

BOOST_AUTO_TEST_CASE(low_frequency_slope) {
    auto shape = sample_electron_shape(FlatElectronShape(100, 1e4));
    auto table = build_numerical_synchrotron_table(shape);
    NumericalElectronDistribution electrons(shape, ElectronNormalization::number, 1e20, 2, 0.1, 1);
    const Real B = unit::Gauss;
    NumericalSynchrotron photons(table, electrons, B);
    const Real nu1 = compute_syn_freq(1, B) * std::exp2(table->log2_u_min - 6);
    const Real nu2 = nu1 * 8;
    const Real slope = std::log(photons.compute_I_nu(nu2) / photons.compute_I_nu(nu1)) / std::log(nu2 / nu1);
    BOOST_CHECK_SMALL(slope - 1.0 / 3.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(power_law_asymptotic_slope) {
    constexpr Real p = 2.4;
    auto shape = sample_electron_shape(PowerLawElectronShape(10, 1e8, p), 40);
    auto table = build_numerical_synchrotron_table(shape, 32);
    NumericalElectronDistribution electrons(shape, ElectronNormalization::number, 1e20, 2, 0.1, 1);
    const Real B = unit::Gauss;
    NumericalSynchrotron photons(table, electrons, B);
    const Real nu1 = compute_syn_freq(1e3, B);
    const Real nu2 = compute_syn_freq(1e6, B);
    const Real slope = std::log(photons.compute_I_nu(nu2) / photons.compute_I_nu(nu1)) / std::log(nu2 / nu1);
    BOOST_TEST_MESSAGE("power-law synchrotron slope error = " << slope + (p - 1) / 2);
    BOOST_CHECK_SMALL(slope + (p - 1) / 2, 3e-3);
}

BOOST_AUTO_TEST_CASE(flat_distribution_emits_without_standard_breaks) {
    auto shape = sample_electron_shape(FlatElectronShape(100, 1e5));
    auto table = build_numerical_synchrotron_table(shape);
    NumericalElectronDistribution electrons(shape, ElectronNormalization::energy, 1e20, 3, 0.1, 1);
    NumericalSynchrotron photons(table, electrons, 0.2 * unit::Gauss);

    BOOST_CHECK_GT(photons.compute_I_nu(compute_syn_freq(1000, 0.2 * unit::Gauss)), 0);
    BOOST_CHECK_GT(photons.compute_I_nu(compute_syn_freq(1e4, 0.2 * unit::Gauss)), 0);
    // Neither object has gamma_m or gamma_c: the static assertions above make accidental standard-break coupling
    // a compile-time failure rather than a spectral-regression convention.
}

BOOST_AUTO_TEST_CASE(cutoff_power_law_suppresses_high_frequency_emission) {
    constexpr Real p = 2.3;
    auto plain_shape = sample_electron_shape(PowerLawElectronShape(10, 1e8, p), 40);
    auto cutoff_shape = sample_electron_shape(CutoffPowerLawElectronShape(10, 1e8, p, 1e5), 40);
    auto plain_table = build_numerical_synchrotron_table(plain_shape, 32);
    auto cutoff_table = build_numerical_synchrotron_table(cutoff_shape, 32);
    NumericalElectronDistribution plain_e(plain_shape, ElectronNormalization::number, 1e20, 2, 0.1, 1);
    NumericalElectronDistribution cutoff_e(cutoff_shape, ElectronNormalization::number, 1e20, 2, 0.1, 1);
    const Real B = unit::Gauss;
    NumericalSynchrotron plain(plain_table, plain_e, B);
    NumericalSynchrotron cutoff(cutoff_table, cutoff_e, B);
    const Real high_nu = compute_syn_freq(1e6, B);
    BOOST_CHECK_LT(cutoff.compute_I_nu(high_nu), 0.2 * plain.compute_I_nu(high_nu));
}

BOOST_AUTO_TEST_CASE(grid_refinement_converges) {
    NumericalElectronShape model = CutoffPowerLawElectronShape(10, 1e7, 2.3, 2e5);
    auto coarse_shape = sample_electron_shape(model, 8);
    auto medium_shape = sample_electron_shape(model, 16);
    auto fine_shape = sample_electron_shape(model, 32);
    auto coarse_table = build_numerical_synchrotron_table(coarse_shape, 8);
    auto medium_table = build_numerical_synchrotron_table(medium_shape, 16);
    auto fine_table = build_numerical_synchrotron_table(fine_shape, 32);
    NumericalElectronDistribution coarse_e(coarse_shape, ElectronNormalization::number, 1e20, 2, 0.1, 1);
    NumericalElectronDistribution medium_e(medium_shape, ElectronNormalization::number, 1e20, 2, 0.1, 1);
    NumericalElectronDistribution fine_e(fine_shape, ElectronNormalization::number, 1e20, 2, 0.1, 1);
    const Real B = unit::Gauss;
    NumericalSynchrotron coarse(coarse_table, coarse_e, B);
    NumericalSynchrotron medium(medium_table, medium_e, B);
    NumericalSynchrotron fine(fine_table, fine_e, B);
    const Real nu = compute_syn_freq(2e5, B);
    const Real fine_I = fine.compute_I_nu(nu);
    const Real coarse_error = std::abs(coarse.compute_I_nu(nu) / fine_I - 1);
    const Real medium_error = std::abs(medium.compute_I_nu(nu) / fine_I - 1);
    BOOST_TEST_MESSAGE("synchrotron refinement errors: coarse=" << coarse_error << ", medium=" << medium_error);
    BOOST_CHECK_LT(medium_error, coarse_error);
    BOOST_CHECK_LT(medium_error, 5e-3);
}

BOOST_AUTO_TEST_CASE(table_build_performance_smoke) {
    auto shape = sample_electron_shape(CutoffPowerLawElectronShape(10, 1e8, 2.3, 1e6), 32);
    const auto start = std::chrono::steady_clock::now();
    auto table = build_numerical_synchrotron_table(shape, 24);
    const auto elapsed = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
    BOOST_TEST_MESSAGE("one numerical synchrotron shape table: " << elapsed << " ms (" << shape->gamma.size()
                                                                  << " gamma x " << table->log2_J.size()
                                                                  << " frequency nodes)");
    BOOST_CHECK_GT(table->log2_J.size(), 1u);
}

BOOST_AUTO_TEST_CASE(zero_and_invalid_inputs) {
    auto shape = sample_electron_shape(FlatElectronShape(10, 100));
    auto table = build_numerical_synchrotron_table(shape);
    NumericalElectronDistribution electrons(shape, ElectronNormalization::number, 1e20, 2, 0.1, 1);
    NumericalSynchrotron zero_B(table, electrons, 0);
    BOOST_CHECK_EQUAL(zero_B.compute_I_nu(1), 0);
    BOOST_CHECK_EQUAL(zero_B.compute_log2_I_nu(0), -con::inf);
    BOOST_CHECK_THROW(NumericalSynchrotron(table, electrons, -1), std::invalid_argument);
    BOOST_CHECK_THROW((void)build_numerical_synchrotron_table(shape, 0), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
