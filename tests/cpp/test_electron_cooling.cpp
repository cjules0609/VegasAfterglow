#include <boost/test/unit_test.hpp>

#include <cmath>
#include <numeric>

#include "radiation/electron-cooling.h"
#include "util/macros.h"

BOOST_AUTO_TEST_SUITE(ElectronCooling)

BOOST_AUTO_TEST_CASE(single_electron_synchrotron_cooling_tracks_characteristic) {
    Array gamma = make_electron_gamma_grid(1, 1e6, 96);
    Array initial = xt::zeros<Real>({gamma.size()});
    constexpr Real gamma0 = 1e4;
    size_t peak = 0;
    for (size_t i = 1; i < gamma.size(); ++i) {
        if (std::abs(std::log(gamma(i) / gamma0)) < std::abs(std::log(gamma(peak) / gamma0))) {
            peak = i;
        }
    }
    initial(peak) = 1;
    ElectronKineticState state(gamma, initial);
    const Real number0 = state.number();

    ElectronCoolingConfig config;
    config.max_step_fraction = 0.2;
    constexpr Real B = 1 * unit::Gauss;
    const Real a = con::sigmaT * B * B / (6 * con::pi * con::me * con::c);
    const Real dt = 0.5 / (a * gamma0);
    advance_electron_distribution(state, dt, B, config);

    Real mean_gamma = 0;
    for (size_t i = 0; i < gamma.size(); ++i) {
        mean_gamma += gamma(i) * state.bin_number(i);
    }
    mean_gamma /= state.number();
    const Real exact = 1 / std::tanh(a * dt + std::atanh(1 / gamma0));
    BOOST_CHECK_SMALL(state.number() / number0 - 1, 2e-12);
    BOOST_CHECK_SMALL(mean_gamma / exact - 1, 0.025);
    BOOST_CHECK_LT(state.kinetic_energy_moment(), (gamma0 - 1) * number0 * 1.02);
}

BOOST_AUTO_TEST_CASE(particle_number_is_conserved_with_closed_lower_boundary) {
    auto shape = sample_electron_shape(PowerLawElectronShape(1, 1e7, 2.3), 32);
    ElectronKineticState state(shape);
    const Real number0 = state.number();
    ElectronCoolingConfig config;
    advance_electron_distribution(state, 1e7 * unit::sec, 10 * unit::Gauss, config);
    BOOST_CHECK_SMALL(state.number() / number0 - 1, 2e-12);
    BOOST_CHECK_EQUAL(state.escaped_lower, 0);
}

BOOST_AUTO_TEST_CASE(continuous_injection_adds_the_integrated_particle_rate) {
    auto shape = sample_electron_shape(FlatElectronShape(1, 1e3), 24);
    ElectronKineticState state(shape);
    const Real number0 = state.number();
    Array injection = xt::ones<Real>({state.gamma.size()});
    const Real rate = std::accumulate(state.bin_width.begin(), state.bin_width.end(), 0.0);
    ElectronCoolingConfig config;
    config.synchrotron = false;
    constexpr Real dt = 3;
    advance_electron_distribution(state, dt, 0, config, 0, 0, injection);
    BOOST_CHECK_SMALL(state.number() / (number0 + dt * rate) - 1, 2e-12);
}

BOOST_AUTO_TEST_CASE(continuously_injected_power_law_steepens_by_one) {
    constexpr Real p = 2.3;
    Array gamma = make_electron_gamma_grid(1, 1e6, 64);
    Array source = xt::zeros<Real>({gamma.size()});
    Array seed = xt::zeros<Real>({gamma.size()});
    for (size_t i = 0; i < gamma.size(); ++i) {
        if (gamma(i) >= 10) {
            source(i) = std::pow(gamma(i), -p);
            seed(i) = 1e-6 * source(i);
        }
    }
    ElectronKineticState state(gamma, seed);
    ElectronCoolingConfig config;
    config.max_step_fraction = 0.3;
    advance_electron_distribution(state, 1e5 * unit::sec, 10 * unit::Gauss, config, 0, 0,
                                  source / unit::sec);

    Array cooled = state.distribution();
    Real sx = 0;
    Real sy = 0;
    Real sxx = 0;
    Real sxy = 0;
    size_t count = 0;
    for (size_t i = 0; i < gamma.size(); ++i) {
        if (gamma(i) >= 300 && gamma(i) <= 3e4 && cooled(i) > 0) {
            const Real x = std::log(gamma(i));
            const Real y = std::log(cooled(i));
            sx += x;
            sy += y;
            sxx += x * x;
            sxy += x * y;
            ++count;
        }
    }
    const Real slope = (count * sxy - sx * sy) / (count * sxx - sx * sx);
    BOOST_CHECK_CLOSE(slope, -(p + 1), 5.0);
}

BOOST_AUTO_TEST_CASE(ic_and_adiabatic_terms_are_independently_configurable) {
    ElectronCoolingConfig none;
    none.synchrotron = false;
    BOOST_CHECK_EQUAL(electron_cooling_rate(100, 1 * unit::Gauss, none), 0);

    ElectronCoolingConfig ic = none;
    ic.inverse_compton = true;
    ic.compton_y = 2;
    BOOST_CHECK_LT(electron_cooling_rate(100, 1 * unit::Gauss, ic), 0);

    ElectronCoolingConfig ad = none;
    ad.adiabatic = true;
    BOOST_CHECK_LT(electron_cooling_rate(100, 0, ad, 1 / unit::sec), 0);
}

BOOST_AUTO_TEST_SUITE_END()
