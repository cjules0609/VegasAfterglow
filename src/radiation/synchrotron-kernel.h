//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include "../util/macros.h"

/**
 * <!-- ************************************************************************************** -->
 * @brief Evaluate the single-electron synchrotron kernel F(x).
 * @details Computes F(x) = x integral_x^infinity K_{5/3}(z) dz from an immutable logarithmic lookup table,
 *          using analytic asymptotes outside the tabulated interval. Invalid, non-positive, and exponentially
 *          underflowed arguments return zero.
 * @param x Dimensionless frequency ratio
 * @return Synchrotron kernel value
 * <!-- ************************************************************************************** -->
 */
[[nodiscard]] Real synchrotron_kernel(Real x) noexcept;

/**
 * <!-- ************************************************************************************** -->
 * @brief Evaluate log2(F(x)) directly from log2(x).
 * @details Avoids underflow in numerical synchrotron integrations. Invalid inputs return negative infinity,
 *          matching the photon-spectrum convention for zero intensity.
 * @param log2_x Base-2 logarithm of the dimensionless frequency ratio
 * @return Base-2 logarithm of the synchrotron kernel
 * <!-- ************************************************************************************** -->
 */
[[nodiscard]] Real log2_synchrotron_kernel(Real log2_x) noexcept;

/**
 * <!-- ************************************************************************************** -->
 * @brief Evaluate the absorption kernel G(x) = F(x) - x F'(x) = x^2 K_{5/3}(x).
 * @details Uses the same immutable logarithmic lookup-table strategy as synchrotron_kernel. This positive
 *          kernel permits integration by parts in the general SSA expression without differentiating a sampled
 *          electron distribution.
 * @param x Dimensionless frequency ratio
 * @return Synchrotron absorption kernel value
 * <!-- ************************************************************************************** -->
 */
[[nodiscard]] Real synchrotron_absorption_kernel(Real x) noexcept;

/**
 * <!-- ************************************************************************************** -->
 * @brief Evaluate log2(G(x)) directly from log2(x).
 * @param log2_x Base-2 logarithm of the dimensionless frequency ratio
 * @return Base-2 logarithm of the absorption kernel, or negative infinity for invalid/underflowed input
 * <!-- ************************************************************************************** -->
 */
[[nodiscard]] Real log2_synchrotron_absorption_kernel(Real log2_x) noexcept;
