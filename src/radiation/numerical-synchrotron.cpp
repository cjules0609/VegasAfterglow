//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "numerical-synchrotron.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

#include "../util/fast-math.h"
#include "synchrotron-kernel.h"
#include "synchrotron.h"

std::shared_ptr<NumericalSynchrotronTable const>
build_numerical_synchrotron_table(std::shared_ptr<ElectronShapeTable const> electrons, size_t points_per_decade,
                                  Real low_x, Real high_x) {
    if (!electrons) {
        throw std::invalid_argument("electron shape table must not be null");
    }
    if (points_per_decade == 0) {
        throw std::invalid_argument("points_per_decade must be greater than zero");
    }
    if (!(std::isfinite(low_x) && low_x > 0 && std::isfinite(high_x) && high_x > low_x)) {
        throw std::invalid_argument("synchrotron kernel margins must be finite, positive, and ordered");
    }

    auto table = std::make_shared<NumericalSynchrotronTable>();
    table->electrons = std::move(electrons);
    table->log2_u_min = fast_log2(low_x) + 2 * fast_log2(table->electrons->gamma_min);
    table->log2_u_max = fast_log2(high_x) + 2 * fast_log2(table->electrons->gamma_max);
    const Real decades = (table->log2_u_max - table->log2_u_min) / std::log2(10.0);
    const size_t intervals =
        std::max<size_t>(1, static_cast<size_t>(std::ceil(decades * static_cast<Real>(points_per_decade))));
    const size_t size = intervals + 1;
    table->log2_u_step = (table->log2_u_max - table->log2_u_min) / static_cast<Real>(intervals);
    table->log2_J = Array::from_shape({size});
    table->absorption = std::make_shared<NumericalSynchrotronAbsorptionCache>();

    const size_t gamma_size = table->electrons->gamma.size();
    static thread_local Array log2_emission_terms;
    log2_emission_terms.resize({gamma_size});
    for (size_t j = 0; j < size; ++j) {
        const Real log2_u = table->log2_u_min + table->log2_u_step * static_cast<Real>(j);
        Real max_emission_term = -con::inf;
        for (size_t i = 0; i < gamma_size; ++i) {
            if (!std::isfinite(table->electrons->log2_shape(i))) {
                continue;
            }
            const Real log2_x = log2_u - 2 * table->electrons->log2_gamma(i);
            const Real common = table->electrons->log2_shape(i) + table->electrons->log2_quadrature_weight(i);
            const Real log2_emission_term = common + log2_synchrotron_kernel(log2_x);
            log2_emission_terms(i) = log2_emission_term;
            max_emission_term = std::max(max_emission_term, log2_emission_term);
        }

        if (!std::isfinite(max_emission_term)) {
            table->log2_J(j) = -con::inf;
        } else {
            Real scaled_sum = 0;
            for (size_t i = 0; i < gamma_size; ++i) {
                if (std::isfinite(table->electrons->log2_shape(i))) {
                    scaled_sum += fast_exp2(log2_emission_terms(i) - max_emission_term);
                }
            }
            table->log2_J(j) = max_emission_term + fast_log2(scaled_sum);
        }
    }
    return table;
}

Array const& ensure_numerical_synchrotron_absorption_table(NumericalSynchrotronTable const& table) {
    if (!table.electrons || !table.absorption) {
        throw std::invalid_argument("numerical synchrotron table has no retained electron shape or absorption cache");
    }
    std::call_once(table.absorption->once, [&table] {
        const size_t size = table.log2_J.size();
        const size_t gamma_size = table.electrons->gamma.size();
        Array log2_K = Array::from_shape({size});
        static thread_local Array log2_absorption_terms;
        log2_absorption_terms.resize({gamma_size});

        for (size_t j = 0; j < size; ++j) {
            const Real log2_u = table.log2_u_min + table.log2_u_step * static_cast<Real>(j);
            Real max_absorption_term = -con::inf;
            for (size_t i = 0; i < gamma_size; ++i) {
                if (!std::isfinite(table.electrons->log2_shape(i))) {
                    continue;
                }
                const Real log2_x = log2_u - 2 * table.electrons->log2_gamma(i);
                const Real common =
                    table.electrons->log2_shape(i) + table.electrons->log2_quadrature_weight(i);
                // The leading one is log2(2), from K(u) = 2 integral f(gamma) G(x) dgamma/gamma.
                const Real term =
                    1 + common - table.electrons->log2_gamma(i) + log2_synchrotron_absorption_kernel(log2_x);
                log2_absorption_terms(i) = term;
                max_absorption_term = std::max(max_absorption_term, term);
            }

            if (!std::isfinite(max_absorption_term)) {
                log2_K(j) = -con::inf;
            } else {
                Real scaled_sum = 0;
                for (size_t i = 0; i < gamma_size; ++i) {
                    if (std::isfinite(table.electrons->log2_shape(i))) {
                        scaled_sum += fast_exp2(log2_absorption_terms(i) - max_absorption_term);
                    }
                }
                log2_K(j) = max_absorption_term + fast_log2(scaled_sum);
            }
        }

        table.absorption->log2_K = std::move(log2_K);
        table.absorption->ready.store(true, std::memory_order_release);
    });
    return table.absorption->log2_K;
}

bool numerical_synchrotron_absorption_ready(NumericalSynchrotronTable const& table) noexcept {
    return table.absorption && table.absorption->ready.load(std::memory_order_acquire);
}

NumericalSynchrotron::NumericalSynchrotron(std::shared_ptr<NumericalSynchrotronTable const> table,
                                           NumericalElectronDistribution const& electrons, Real B, bool ssa_enabled)
    : table(std::move(table)), ssa_enabled(ssa_enabled) {
    if (!this->table) {
        throw std::invalid_argument("numerical synchrotron table must not be null");
    }
    if (this->table->electrons != electrons.shape) {
        throw std::invalid_argument("numerical synchrotron table and electron distribution must share one shape");
    }
    if (!(std::isfinite(B) && B >= 0)) {
        throw std::invalid_argument("B must be finite and non-negative");
    }
    if (ssa_enabled && !numerical_synchrotron_absorption_ready(*this->table)) {
        (void)ensure_numerical_synchrotron_absorption_table(*this->table);
    }
    if (B == 0 || electrons.column_scale == 0) {
        return;
    }

    // Keep the same characteristic-frequency definition as compute_syn_freq(gamma, B).
    log2_frequency_scale = fast_log2(compute_syn_freq(1, B));

    // Match compute_single_elec_P_nu_max: effective <sin(alpha)> = pi/4 and F_max = 0.92.
    // Here the full kernel F(x), rather than F_max, is carried by the shared J(u) table.
    constexpr Real sin_angle_ave = con::pi / 4;
    const Real power_prefactor =
        std::numbers::sqrt3 * con::e3 * B * sin_angle_ave / (con::me * con::c2);
    const Real intensity_scale = electrons.column_scale * power_prefactor / (4 * con::pi);
    log2_intensity_scale = fast_log2(intensity_scale);

    // Column form of synchrotron absorption:
    // tau_nu = A_cell P0 B K(u) / (8 pi m_e nu^2), with the factor two carried by K.
    const Real optical_depth_scale = electrons.column_scale * power_prefactor / (8 * con::pi * con::me);
    log2_optical_depth_scale = fast_log2(optical_depth_scale);
}

Real NumericalSynchrotron::compute_I_nu(Real nu) const noexcept {
    if (!(nu > 0) || !std::isfinite(nu)) {
        return 0;
    }
    const Real log2_I_nu = compute_log2_I_nu(fast_log2(nu));
    constexpr Real log2_smallest = -1074;
    return (log2_I_nu > log2_smallest) ? fast_exp2(log2_I_nu) : 0;
}

Real NumericalSynchrotron::compute_log2_I_nu(Real log2_nu) const noexcept {
    const Real log2_I_thin = compute_log2_thin_I_nu(log2_nu);
    if (!ssa_enabled || !std::isfinite(log2_I_thin)) {
        return log2_I_thin;
    }

    // An SSA-enabled object initializes the shared cache in its constructor, so the hot spectrum path needs no
    // synchronization or readiness check.
    const Real log2_tau = compute_log2_optical_depth_from_table(log2_nu, table->absorption->log2_K);
    if (!std::isfinite(log2_tau)) {
        return log2_I_thin;
    }

    // Below this point 1 - R(tau) is below double-precision resolution.
    if (log2_tau < -52) {
        return log2_I_thin;
    }
    // exp(-tau) is negligible above tau=64; staying in log space also avoids overflow.
    if (log2_tau > 6) {
        return log2_I_thin - log2_tau;
    }

    const Real tau = fast_exp2(log2_tau);
    const Real transfer = -std::expm1(-tau) / tau;
    return log2_I_thin + fast_log2(transfer);
}

Real NumericalSynchrotron::compute_log2_thin_I_nu(Real log2_nu) const noexcept {
    if (!table || !std::isfinite(log2_nu) || !std::isfinite(log2_frequency_scale) ||
        !std::isfinite(log2_intensity_scale)) {
        return -con::inf;
    }

    const Real log2_u = log2_nu - log2_frequency_scale;
    if (log2_u < table->log2_u_min) {
        // Preserve the historical optically-thin operation ordering exactly when SSA is disabled.
        return log2_intensity_scale + table->log2_J(0) + (log2_u - table->log2_u_min) / 3;
    }
    if (log2_u > table->log2_u_max) {
        return -con::inf;
    }
    const size_t size = table->log2_J.size();
    if (log2_u == table->log2_u_max) {
        return log2_intensity_scale + table->log2_J(size - 1);
    }
    const Real position = (log2_u - table->log2_u_min) / table->log2_u_step;
    const size_t i = std::min(static_cast<size_t>(position), size - 2);
    const Real fraction = position - static_cast<Real>(i);
    return log2_intensity_scale + table->log2_J(i) +
           fraction * (table->log2_J(i + 1) - table->log2_J(i));
}

Real NumericalSynchrotron::compute_optical_depth(Real nu) const {
    if (!(nu > 0) || !std::isfinite(nu)) {
        return 0;
    }
    const Real log2_tau = compute_log2_optical_depth(fast_log2(nu));
    if (log2_tau > 1023) {
        return con::inf;
    }
    constexpr Real log2_smallest = -1074;
    return (log2_tau > log2_smallest) ? fast_exp2(log2_tau) : 0;
}

Real NumericalSynchrotron::compute_log2_optical_depth(Real log2_nu) const {
    if (!table || !std::isfinite(log2_nu) || !std::isfinite(log2_frequency_scale) ||
        !std::isfinite(log2_optical_depth_scale)) {
        return -con::inf;
    }
    const Array& absorption = numerical_synchrotron_absorption_ready(*table)
                                  ? table->absorption->log2_K
                                  : ensure_numerical_synchrotron_absorption_table(*table);
    return compute_log2_optical_depth_from_table(log2_nu, absorption);
}

Real NumericalSynchrotron::compute_log2_optical_depth_from_table(Real log2_nu,
                                                                 Array const& absorption) const noexcept {
    const Real log2_u = log2_nu - log2_frequency_scale;
    const Real log2_K = interpolate_log2_table(absorption, log2_u, 1.0 / 3.0);
    return std::isfinite(log2_K) ? log2_optical_depth_scale - 2 * log2_nu + log2_K : -con::inf;
}

Real NumericalSynchrotron::interpolate_log2_table(Array const& values, Real log2_u, Real low_slope) const noexcept {
    if (log2_u < table->log2_u_min) {
        return values(0) + low_slope * (log2_u - table->log2_u_min);
    }
    if (log2_u > table->log2_u_max) {
        return -con::inf;
    }

    const size_t size = values.size();
    if (log2_u == table->log2_u_max) {
        return values(size - 1);
    }
    const Real position = (log2_u - table->log2_u_min) / table->log2_u_step;
    const size_t i = std::min(static_cast<size_t>(position), size - 2);
    const Real f = position - static_cast<Real>(i);
    return values(i) + f * (values(i + 1) - values(i));
}
