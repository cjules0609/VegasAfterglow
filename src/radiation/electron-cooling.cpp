#include "electron-cooling.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace {
    void validate_config(ElectronCoolingConfig const& config) {
        if (!(std::isfinite(config.compton_y) && config.compton_y >= 0)) {
            throw std::invalid_argument("compton_y must be finite and non-negative");
        }
        if (!(std::isfinite(config.max_step_fraction) && config.max_step_fraction > 0)) {
            throw std::invalid_argument("max_step_fraction must be finite and positive");
        }
        if (config.max_substeps == 0) {
            throw std::invalid_argument("max_substeps must be positive");
        }
        if (config.inverse_compton && config.compton_y == 0) {
            // An external photon energy density may still be supplied at advance time.
            return;
        }
    }

    Array make_bin_edges(Array const& gamma) {
        if (gamma.dimension() != 1 || gamma.size() < 2) {
            throw std::invalid_argument("electron gamma grid must be one-dimensional with at least two samples");
        }
        Array edges = Array::from_shape({gamma.size() + 1});
        for (size_t i = 0; i < gamma.size(); ++i) {
            if (!(std::isfinite(gamma(i)) && gamma(i) >= 1 && (i == 0 || gamma(i) > gamma(i - 1)))) {
                throw std::invalid_argument("electron gamma grid must be finite, strictly increasing, and >= 1");
            }
        }
        edges(0) = gamma(0);
        for (size_t i = 1; i < gamma.size(); ++i) {
            edges(i) = std::sqrt(gamma(i - 1) * gamma(i));
        }
        edges(gamma.size()) = gamma(gamma.size() - 1);
        return edges;
    }

    Real stable_substep(ElectronKineticState const& state, Real magnetic_field,
                        ElectronCoolingConfig const& config, Real expansion_rate,
                        Real photon_energy_density) {
        Real result = con::inf;
        for (size_t i = 0; i < state.gamma.size(); ++i) {
            const Real speed = std::abs(electron_cooling_rate(state.gamma(i), magnetic_field, config,
                                                               expansion_rate, photon_energy_density));
            if (speed > 0) {
                result = std::min(result, config.max_step_fraction * state.bin_width(i) / speed);
            }
        }
        return result;
    }
} // namespace

ElectronKineticState::ElectronKineticState(std::shared_ptr<ElectronShapeTable const> distribution)
    : ElectronKineticState(distribution ? distribution->gamma : Array{},
                           distribution ? distribution->shape : Array{}) {
    if (!distribution) {
        throw std::invalid_argument("electron distribution must not be null");
    }
}
ElectronKineticState::ElectronKineticState(Array const& gamma, Array const& distribution)
    : gamma(gamma), bin_edges(make_bin_edges(gamma)), bin_width(Array::from_shape({gamma.size()})),
      bin_number(Array::from_shape({gamma.size()})) {
    if (distribution.dimension() != 1 || distribution.size() != gamma.size()) {
        throw std::invalid_argument("electron distribution must match the one-dimensional gamma grid");
    }
    Real total = 0;
    for (size_t i = 0; i < gamma.size(); ++i) {
        if (!(std::isfinite(distribution(i)) && distribution(i) >= 0)) {
            throw std::invalid_argument("electron distribution must be finite and non-negative");
        }
        bin_width(i) = bin_edges(i + 1) - bin_edges(i);
        bin_number(i) = distribution(i) * bin_width(i);
        total += bin_number(i);
    }
    if (!(std::isfinite(total) && total > 0)) {
        throw std::invalid_argument("electron distribution must contain a finite positive particle number");
    }
}

Array ElectronKineticState::distribution() const {
    Array result = Array::from_shape({gamma.size()});
    for (size_t i = 0; i < gamma.size(); ++i) {
        result(i) = bin_number(i) / bin_width(i);
    }
    return result;
}

Real ElectronKineticState::number() const noexcept {
    return std::accumulate(bin_number.begin(), bin_number.end(), 0.0);
}

Real ElectronKineticState::kinetic_energy_moment() const noexcept {
    Real result = 0;
    for (size_t i = 0; i < gamma.size(); ++i) {
        result += (gamma(i) - 1) * bin_number(i);
    }
    return result;
}

Real electron_cooling_rate(Real gamma, Real magnetic_field, ElectronCoolingConfig const& config,
                           Real expansion_rate, Real photon_energy_density) noexcept {
    if (!(gamma >= 1) || !std::isfinite(gamma)) {
        return 0;
    }
    const Real momentum2 = (gamma - 1) * (gamma + 1);
    Real rate = 0;
    if (config.synchrotron && magnetic_field > 0 && std::isfinite(magnetic_field)) {
        const Real a_syn = con::sigmaT * magnetic_field * magnetic_field / (6 * con::pi * con::me * con::c);
        rate -= a_syn * momentum2;
    }
    if (config.inverse_compton) {
        if (config.compton_y > 0 && magnetic_field > 0 && std::isfinite(magnetic_field)) {
            const Real a_syn = con::sigmaT * magnetic_field * magnetic_field / (6 * con::pi * con::me * con::c);
            rate -= config.compton_y * a_syn * momentum2;
        }
        if (photon_energy_density > 0 && std::isfinite(photon_energy_density)) {
            const Real a_ic = 4 * con::sigmaT * photon_energy_density / (3 * con::me * con::c);
            rate -= a_ic * momentum2;
        }
    }
    if (config.adiabatic && expansion_rate > 0 && std::isfinite(expansion_rate)) {
        rate -= expansion_rate * momentum2 / (3 * gamma);
    }
    return rate;
}

void advance_electron_distribution(ElectronKineticState& state, Real dt, Real magnetic_field,
                                   ElectronCoolingConfig const& config, Real expansion_rate,
                                   Real photon_energy_density, Array const& injection, Real escape_time) {
    validate_config(config);
    if (!(std::isfinite(dt) && dt >= 0)) {
        throw std::invalid_argument("cooling timestep must be finite and non-negative");
    }
    if (!(std::isfinite(magnetic_field) && magnetic_field >= 0)) {
        throw std::invalid_argument("magnetic field must be finite and non-negative");
    }
    if (!(std::isfinite(expansion_rate) && expansion_rate >= 0)) {
        throw std::invalid_argument("expansion_rate must be finite and non-negative");
    }
    if (!(std::isfinite(photon_energy_density) && photon_energy_density >= 0)) {
        throw std::invalid_argument("photon_energy_density must be finite and non-negative");
    }
    if (!(escape_time > 0) || std::isnan(escape_time)) {
        throw std::invalid_argument("escape_time must be positive or infinity");
    }
    const bool has_injection = injection.size() != 0;
    if (has_injection && (injection.dimension() != 1 || injection.size() != state.gamma.size())) {
        throw std::invalid_argument("injection must be empty or match the electron gamma grid");
    }
    if (dt == 0) {
        return;
    }
    if (has_injection) {
        for (Real value : injection) {
            if (!(std::isfinite(value) && value >= 0)) {
                throw std::invalid_argument("injection must be finite and non-negative");
            }
        }
    }

    const Real suggested = stable_substep(state, magnetic_field, config, expansion_rate, photon_energy_density);
    size_t substeps = 1;
    if (std::isfinite(suggested) && suggested > 0) {
        substeps = std::max<size_t>(1, static_cast<size_t>(std::ceil(dt / suggested)));
        substeps = std::min(substeps, config.max_substeps);
    }
    const Real h = dt / static_cast<Real>(substeps);
    const size_t n = state.gamma.size();
    Array velocity = Array::from_shape({n + 1});
    Array next = Array::from_shape({n});

    for (size_t step = 0; step < substeps; ++step) {
        for (size_t edge = 0; edge <= n; ++edge) {
            velocity(edge) = electron_cooling_rate(state.bin_edges(edge), magnetic_field, config,
                                                   expansion_rate, photon_energy_density);
        }
        if (config.accumulate_at_gamma_min) {
            velocity(0) = 0;
        }

        // Backward-Euler donor-cell transport. Cooling velocities are non-positive,
        // so the matrix is upper bidiagonal and can be solved from high to low gamma.
        for (size_t reverse = 0; reverse < n; ++reverse) {
            const size_t i = n - 1 - reverse;
            Real rhs = state.bin_number(i);
            if (has_injection) {
                rhs += h * injection(i) * state.bin_width(i);
            }
            const Real escape = std::isfinite(escape_time) ? h / escape_time : 0;
            const Real diagonal = 1 + escape - h * velocity(i) / state.bin_width(i);
            const Real upper = (i + 1 < n) ? h * velocity(i + 1) / state.bin_width(i + 1) : 0;
            const Real coupled = (i + 1 < n) ? upper * next(i + 1) : 0;
            next(i) = std::max(0.0, (rhs - coupled) / diagonal);
        }
        if (!config.accumulate_at_gamma_min && velocity(0) < 0) {
            state.escaped_lower += -h * velocity(0) * next(0) / state.bin_width(0);
        }
        state.bin_number = next;
    }
    state.time += dt;
}
