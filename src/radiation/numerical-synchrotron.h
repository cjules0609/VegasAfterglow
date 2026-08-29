//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <atomic>
#include <memory>
#include <mutex>

#include "electron-distribution.h"

/**
 * <!-- ************************************************************************************** -->
 * @struct NumericalSynchrotronAbsorptionCache
 * @brief Thread-safe lazy storage for a dimensionless synchrotron absorption table.
 * @details This state remains behind shared ownership so NumericalSynchrotronTable stays copyable. The once flag
 *          serializes the first K(u) construction, while ready permits inexpensive state inspection and avoids
 *          repeated call_once operations during per-cell construction.
 * <!-- ************************************************************************************** -->
 */
struct NumericalSynchrotronAbsorptionCache {
    std::once_flag once;
    std::atomic<bool> ready{false};
    Array log2_K;
};

/**
 * <!-- ************************************************************************************** -->
 * @struct NumericalSynchrotronTable
 * @brief Shared dimensionless synchrotron emission table and lazy absorption state for one electron shape.
 * @details Stores J(u) = integral f(gamma) F(u/gamma^2) dgamma eagerly. The matching
 *          K(u) = 2 integral f(gamma) G(u/gamma^2) dgamma/gamma is constructed on first use from the retained
 *          immutable electron table. The default domain extends from u = 10^-6 gamma_min^2 to
 *          u = 50 gamma_max^2, covering the low-frequency power-law tail and the exponentially suppressed
 *          high-frequency tail of the single-electron kernel.
 * <!-- ************************************************************************************** -->
 */
struct NumericalSynchrotronTable {
    std::shared_ptr<ElectronShapeTable const> electrons;
    Array log2_J;
    std::shared_ptr<NumericalSynchrotronAbsorptionCache> absorption;
    Real log2_u_min{0};
    Real log2_u_step{0};
    Real log2_u_max{0};
};

/**
 * <!-- ************************************************************************************** -->
 * @brief Build a shared dimensionless synchrotron spectrum for an electron shape.
 * @param electrons Immutable sampled electron shape
 * @param points_per_decade Approximate dimensionless-frequency nodes per decade
 * @param low_x Lowest single-electron kernel ratio retained at gamma_min
 * @param high_x Highest single-electron kernel ratio retained at gamma_max
 * @return Shared immutable numerical synchrotron table
 * <!-- ************************************************************************************** -->
 */
[[nodiscard]] std::shared_ptr<NumericalSynchrotronTable const>
build_numerical_synchrotron_table(std::shared_ptr<ElectronShapeTable const> electrons,
                                  size_t points_per_decade = 24, Real low_x = 1e-6, Real high_x = 50);

/** @brief Construct K(u) once from the retained native shape, then return the cached table. */
[[nodiscard]] Array const& ensure_numerical_synchrotron_absorption_table(NumericalSynchrotronTable const& table);

/** @brief Return whether K(u) has already been constructed without triggering construction. */
[[nodiscard]] bool numerical_synchrotron_absorption_ready(NumericalSynchrotronTable const& table) noexcept;

/**
 * <!-- ************************************************************************************** -->
 * @struct NumericalSynchrotron
 * @brief Lightweight per-cell numerical synchrotron photon spectrum.
 * @details Holds a shared dimensionless spectrum, the magnetic frequency shift, and cell-specific intensity and
 *          optical-depth scales. Self-absorption is explicit radiative transfer through the electron column; it
 *          introduces no gamma_a or standard afterglow break. It intentionally has no electron index, cooling
 *          regime, or inverse-Compton state.
 * <!-- ************************************************************************************** -->
 */
struct NumericalSynchrotron {
    NumericalSynchrotron() noexcept = default;

    NumericalSynchrotron(std::shared_ptr<NumericalSynchrotronTable const> table,
                         NumericalElectronDistribution const& electrons, Real B, bool ssa_enabled = false);

    std::shared_ptr<NumericalSynchrotronTable const> table;
    Real log2_frequency_scale{-con::inf};
    Real log2_intensity_scale{-con::inf};
    Real log2_optical_depth_scale{-con::inf};
    bool ssa_enabled{false};

    [[nodiscard]] Real compute_I_nu(Real nu) const noexcept;
    [[nodiscard]] Real compute_log2_I_nu(Real log2_nu) const noexcept;
    [[nodiscard]] Real compute_optical_depth(Real nu) const;
    [[nodiscard]] Real compute_log2_optical_depth(Real log2_nu) const;

  private:
    [[nodiscard]] Real compute_log2_thin_I_nu(Real log2_nu) const noexcept;
    [[nodiscard]] Real compute_log2_optical_depth_from_table(Real log2_nu, Array const& absorption) const noexcept;
    [[nodiscard]] Real interpolate_log2_table(Array const& values, Real log2_u, Real low_slope) const noexcept;
};
