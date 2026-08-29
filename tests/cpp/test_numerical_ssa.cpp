#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <future>
#include <vector>

#include "radiation/electron-distribution.h"
#include "radiation/numerical-synchrotron.h"
#include "radiation/synchrotron-kernel.h"
#include "radiation/synchrotron.h"
#include "util/macros.h"

namespace {
    Real slope(Real x1, Real y1, Real x2, Real y2) {
        return std::log(y2 / y1) / std::log(x2 / x1);
    }

    Real interpolated_K(NumericalSynchrotronTable const& table, Real u) {
        const Array& log2_K = ensure_numerical_synchrotron_absorption_table(table);
        const Real log2_u = std::log2(u);
        if (log2_u < table.log2_u_min) {
            return std::exp2(log2_K(0) + (log2_u - table.log2_u_min) / 3);
        }
        const Real position = (log2_u - table.log2_u_min) / table.log2_u_step;
        const size_t i = std::min(static_cast<size_t>(position), log2_K.size() - 2);
        const Real fraction = position - static_cast<Real>(i);
        return std::exp2(log2_K(i) + fraction * (log2_K(i + 1) - log2_K(i)));
    }
} // namespace

BOOST_AUTO_TEST_SUITE(NumericalSynchrotronSSA)

BOOST_AUTO_TEST_CASE(all_builtin_shapes_start_without_absorption) {
    const std::array<NumericalElectronShape, 3> shapes = {
        FlatElectronShape(10, 1e6),
        PowerLawElectronShape(10, 1e7, 2.3),
        CutoffPowerLawElectronShape(10, 1e7, 2.3, 1e5),
    };
    for (auto const& shape : shapes) {
        auto sampled = sample_electron_shape(shape);
        auto table = build_numerical_synchrotron_table(sampled);
        BOOST_CHECK(!numerical_synchrotron_absorption_ready(*table));
    }
}

BOOST_AUTO_TEST_CASE(absorption_table_is_lazy_shared_and_constructed_once) {
    auto shape = sample_electron_shape(PowerLawElectronShape(10, 1e7, 2.3), 48);
    auto table = build_numerical_synchrotron_table(shape, 40);
    BOOST_CHECK(!numerical_synchrotron_absorption_ready(*table));

    NumericalElectronDistribution electrons(shape, ElectronNormalization::number, 1e20 / unit::cm2, 2, 0.1, 1);
    NumericalSynchrotron thin(table, electrons, unit::Gauss, false);
    const Real nu = compute_syn_freq(1e4, unit::Gauss);
    for (size_t i = 0; i < 100; ++i) {
        BOOST_CHECK(thin.compute_I_nu(nu) > 0);
    }
    BOOST_CHECK(!numerical_synchrotron_absorption_ready(*table));

    const Real tau = thin.compute_optical_depth(nu);
    BOOST_CHECK(tau > 0);
    BOOST_CHECK(numerical_synchrotron_absorption_ready(*table));
    auto const* cached = &ensure_numerical_synchrotron_absorption_table(*table);
    for (size_t i = 0; i < 100; ++i) {
        BOOST_CHECK_EQUAL(thin.compute_optical_depth(nu), tau);
        BOOST_CHECK_EQUAL(&ensure_numerical_synchrotron_absorption_table(*table), cached);
    }

    // Copying table/configuration values retains the same non-copyable synchronization state behind shared
    // ownership; copying a Radiation or Model configuration therefore cannot duplicate K(u).
    NumericalSynchrotronTable copied = *table;
    BOOST_CHECK_EQUAL(copied.absorption.get(), table->absorption.get());
    BOOST_CHECK_EQUAL(&ensure_numerical_synchrotron_absorption_table(copied), cached);
}

BOOST_AUTO_TEST_CASE(concurrent_first_absorption_use_is_thread_safe) {
    auto shape = sample_electron_shape(PowerLawElectronShape(10, 1e7, 2.3), 48);
    auto table = build_numerical_synchrotron_table(shape, 40);
    NumericalElectronDistribution electrons(shape, ElectronNormalization::number, 1e20 / unit::cm2, 2, 0.1, 1);
    const Real nu = compute_syn_freq(1e4, unit::Gauss);

    std::vector<std::future<Real>> futures;
    for (size_t i = 0; i < 8; ++i) {
        futures.emplace_back(std::async(std::launch::async, [table, &electrons, nu] {
            NumericalSynchrotron photons(table, electrons, unit::Gauss, true);
            return photons.compute_optical_depth(nu);
        }));
    }
    const Real reference = futures.front().get();
    BOOST_CHECK(reference > 0);
    for (size_t i = 1; i < futures.size(); ++i) {
        BOOST_CHECK_EQUAL(futures[i].get(), reference);
    }
    BOOST_CHECK(numerical_synchrotron_absorption_ready(*table));
}

BOOST_AUTO_TEST_CASE(flat_hard_boundaries_match_distributional_derivative) {
    constexpr std::array<std::array<Real, 3>, 5> cases = {{{1e-4, 2, 20},
                                                           {1, 10, 1e3},
                                                           {1e4, 3, 50},
                                                           {1e8, 100, 1e5},
                                                           {1e12, 1e2, 1e8}}};
    for (auto const& test : cases) {
        const Real u = test[0];
        const Real gamma_min = test[1];
        const Real gamma_max = test[2];
        auto shape = sample_electron_shape(FlatElectronShape(gamma_min, gamma_max), 256);
        auto table = build_numerical_synchrotron_table(shape, 96);

        Real interior = 0;
        for (size_t i = 0; i < shape->gamma.size(); ++i) {
            const Real gamma = shape->gamma(i);
            interior += 2 * shape->quadrature_weight(i) / gamma * synchrotron_kernel(u / (gamma * gamma));
        }
        // q=f/gamma^2 has jumps +delta(gamma-a)/a^2 and -delta(gamma-b)/b^2.
        const Real boundary = -synchrotron_kernel(u / (gamma_min * gamma_min)) +
                              synchrotron_kernel(u / (gamma_max * gamma_max));
        const Real direct = interior + boundary;
        const Real by_parts = interpolated_K(*table, u);
        BOOST_TEST_MESSAGE("flat-boundary check u=" << u << " gamma=[" << gamma_min << ", " << gamma_max
                                                     << "] relative error=" << by_parts / direct - 1);
        // The independently evaluated high-precision check agrees to 4e-14; this sampled-grid comparison also
        // includes the production quadrature error near an exponentially sharp endpoint.
        BOOST_CHECK_SMALL(by_parts / direct - 1, 3e-4);
    }
}

BOOST_AUTO_TEST_CASE(optical_depth_column_and_fixed_u_magnetic_scaling) {
    auto shape = sample_electron_shape(PowerLawElectronShape(10, 1e7, 2.3), 48);
    auto table = build_numerical_synchrotron_table(shape, 40);
    NumericalElectronDistribution e1(shape, ElectronNormalization::number, 1e20 / unit::cm2, 2, 0.1, 0.2);
    NumericalElectronDistribution e2(shape, ElectronNormalization::number, 1e20 / unit::cm2, 2, 0.1, 0.6);
    const Real B1 = 0.2 * unit::Gauss;
    const Real B2 = 0.8 * unit::Gauss;
    NumericalSynchrotron p1(table, e1, B1, true);
    NumericalSynchrotron p2(table, e2, B1, true);
    NumericalSynchrotron pB(table, e1, B2, true);
    const Real u = 1e8;
    const Real tau1 = p1.compute_optical_depth(compute_syn_freq(1, B1) * u);
    const Real tau2 = p2.compute_optical_depth(compute_syn_freq(1, B1) * u);
    const Real tauB = pB.compute_optical_depth(compute_syn_freq(1, B2) * u);
    BOOST_CHECK_SMALL(tau2 / tau1 - 3, 3e-12);
    BOOST_CHECK_SMALL(tauB / tau1 - B1 / B2, 3e-12);
}

BOOST_AUTO_TEST_CASE(power_law_absorption_asymptotes) {
    constexpr Real p = 2.4;
    auto shape = sample_electron_shape(PowerLawElectronShape(10, 1e9, p), 64);
    auto table = build_numerical_synchrotron_table(shape, 48);
    NumericalElectronDistribution electrons(shape, ElectronNormalization::number, 1e20 / unit::cm2, 2, 0.1, 1);
    const Real B1 = 0.1 * unit::Gauss;
    const Real B2 = 0.4 * unit::Gauss;
    NumericalSynchrotron ph1(table, electrons, B1, true);
    NumericalSynchrotron ph2(table, electrons, B2, true);
    const Real nu1 = compute_syn_freq(1e3, B1);
    const Real nu2 = compute_syn_freq(1e6, B1);
    const Real nu_slope = slope(nu1, ph1.compute_optical_depth(nu1), nu2, ph1.compute_optical_depth(nu2));
    BOOST_CHECK_SMALL(nu_slope + (p + 4) / 2, 3e-3);

    const Real nu_fixed = compute_syn_freq(1e5, B1);
    const Real magnetic_ratio = ph2.compute_optical_depth(nu_fixed) / ph1.compute_optical_depth(nu_fixed);
    BOOST_CHECK_SMALL(magnetic_ratio / std::pow(B2 / B1, (p + 2) / 2) - 1, 3e-3);
}

BOOST_AUTO_TEST_CASE(transfer_limits_are_stable) {
    auto shape = sample_electron_shape(PowerLawElectronShape(10, 1e7, 2.3), 48);
    auto table = build_numerical_synchrotron_table(shape, 40);
    const Real B = unit::Gauss;
    const Real nu = compute_syn_freq(1e4, B);

    NumericalElectronDistribution trial(shape, ElectronNormalization::number, 1e20 / unit::cm2, 2, 0.1, 1);
    NumericalSynchrotron trial_ph(table, trial, B, true);
    const Real trial_tau = trial_ph.compute_optical_depth(nu);

    const Real thin_column = (1e20 / unit::cm2) * 1e-10 / trial_tau;
    NumericalElectronDistribution thin_e(shape, ElectronNormalization::number, thin_column, 2, 0.1, 1);
    NumericalSynchrotron thin(table, thin_e, B, false);
    NumericalSynchrotron thin_ssa(table, thin_e, B, true);
    const Real thin_tau = thin_ssa.compute_optical_depth(nu);
    BOOST_CHECK_CLOSE(thin_tau, 1e-10, 1e-8);
    BOOST_CHECK_SMALL(thin_ssa.compute_I_nu(nu) / thin.compute_I_nu(nu) - 1, 1e-9);

    const Real thick_column = (1e20 / unit::cm2) * 1e8 / trial_tau;
    NumericalElectronDistribution thick_e(shape, ElectronNormalization::number, thick_column, 2, 0.1, 1);
    NumericalSynchrotron thick(table, thick_e, B, false);
    NumericalSynchrotron thick_ssa(table, thick_e, B, true);
    const Real thick_tau = thick_ssa.compute_optical_depth(nu);
    BOOST_CHECK_CLOSE(thick_tau, 1e8, 1e-8);
    BOOST_CHECK_SMALL(thick_ssa.compute_I_nu(nu) / (thick.compute_I_nu(nu) / thick_tau) - 1, 2e-12);
}

BOOST_AUTO_TEST_CASE(optically_thick_power_law_slope_is_five_halves) {
    constexpr Real p = 2.4;
    auto shape = sample_electron_shape(PowerLawElectronShape(10, 1e9, p), 64);
    auto table = build_numerical_synchrotron_table(shape, 48);
    const Real B = unit::Gauss;
    const Real nu1 = compute_syn_freq(1e3, B);
    const Real nu2 = compute_syn_freq(1e5, B);
    NumericalElectronDistribution trial(shape, ElectronNormalization::number, 1e20 / unit::cm2, 2, 0.1, 1);
    NumericalSynchrotron trial_ph(table, trial, B, true);
    const Real min_tau = std::min(trial_ph.compute_optical_depth(nu1), trial_ph.compute_optical_depth(nu2));
    NumericalElectronDistribution thick(shape, ElectronNormalization::number,
                                         (1e20 / unit::cm2) * 1e6 / min_tau, 2, 0.1, 1);
    NumericalSynchrotron photons(table, thick, B, true);
    const Real spectral_slope = slope(nu1, photons.compute_I_nu(nu1), nu2, photons.compute_I_nu(nu2));
    BOOST_CHECK_SMALL(spectral_slope - 2.5, 4e-3);
}

BOOST_AUTO_TEST_CASE(deep_low_frequency_thick_slope_is_two) {
    auto shape = sample_electron_shape(FlatElectronShape(100, 1e5), 64);
    auto table = build_numerical_synchrotron_table(shape, 48);
    const Real B = unit::Gauss;
    const Real nu1 = compute_syn_freq(1, B) * std::exp2(table->log2_u_min - 12);
    const Real nu2 = 64 * nu1;
    NumericalElectronDistribution trial(shape, ElectronNormalization::number, 1e20 / unit::cm2, 2, 0.1, 1);
    NumericalSynchrotron trial_ph(table, trial, B, true);
    const Real min_tau = std::min(trial_ph.compute_optical_depth(nu1), trial_ph.compute_optical_depth(nu2));
    NumericalElectronDistribution thick(shape, ElectronNormalization::number,
                                         (1e20 / unit::cm2) * 1e6 / min_tau, 2, 0.1, 1);
    NumericalSynchrotron photons(table, thick, B, true);
    BOOST_CHECK_SMALL(slope(nu1, photons.compute_I_nu(nu1), nu2, photons.compute_I_nu(nu2)) - 2, 2e-12);
}

BOOST_AUTO_TEST_SUITE_END()
