// Conservative kinetic evolution of isotropic electron distributions.
#pragma once

#include <cstddef>
#include <memory>

#include "electron-distribution.h"

/** Configuration for comoving one-zone electron energy losses. */
struct ElectronCoolingConfig {
    bool synchrotron{true};
    bool inverse_compton{false};
    bool adiabatic{false};
    Real compton_y{0};
    Real max_step_fraction{0.5};
    size_t max_substeps{4096};
    bool accumulate_at_gamma_min{true};
};

/**
 * Finite-volume state for N(gamma)=dN/dgamma on a fixed logarithmic grid.
 *
 * `bin_number` contains bin-integrated particle numbers.  Consequently the
 * transport step is particle-conservative (apart from explicit escape and an
 * open lower boundary) and never renormalizes the population after cooling.
 */
struct ElectronKineticState {
    explicit ElectronKineticState(std::shared_ptr<ElectronShapeTable const> distribution);
    ElectronKineticState(Array const& gamma, Array const& distribution);

    Array gamma;
    Array bin_edges;
    Array bin_width;
    Array bin_number;
    Real time{0};
    Real escaped_lower{0};

    [[nodiscard]] Array distribution() const;
    [[nodiscard]] Real number() const noexcept;
    [[nodiscard]] Real kinetic_energy_moment() const noexcept;
};

/** Return dgamma/dt' in internal units at one Lorentz factor. */
[[nodiscard]] Real electron_cooling_rate(Real gamma, Real magnetic_field,
                                         ElectronCoolingConfig const& config,
                                         Real expansion_rate = 0,
                                         Real photon_energy_density = 0) noexcept;

/**
 * Advance one comoving interval with fields held constant over the interval.
 * `injection` is Q(gamma)=dN/(dgamma dt') sampled at the state centers.
 * `escape_time` is a comoving time; infinity disables escape.
 */
void advance_electron_distribution(ElectronKineticState& state, Real dt, Real magnetic_field,
                                   ElectronCoolingConfig const& config,
                                   Real expansion_rate = 0,
                                   Real photon_energy_density = 0,
                                   Array const& injection = Array{},
                                   Real escape_time = con::inf);
