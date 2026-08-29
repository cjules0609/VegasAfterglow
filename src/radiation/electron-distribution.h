//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <memory>
#include <variant>

#include "../core/mesh.h"

/**
 * <!-- ************************************************************************************** -->
 * @enum ElectronNormalization
 * @brief Selects which shock constraint normalizes an instantaneous numerical electron distribution.
 * @details Fixed distribution bounds generally prevent eps_e and xi_e from being satisfied simultaneously.
 *          Every numerical distribution must therefore select energy or number normalization explicitly.
 * <!-- ************************************************************************************** -->
 */
enum class ElectronNormalization { energy, number };

/**
 * <!-- ************************************************************************************** -->
 * @struct FlatElectronShape
 * @brief Bounded dimensionless electron shape f(gamma) = 1.
 * <!-- ************************************************************************************** -->
 */
struct FlatElectronShape {
    FlatElectronShape(Real gamma_min, Real gamma_max);

    Real gamma_min;
    Real gamma_max;

    [[nodiscard]] Real evaluate(Real gamma) const noexcept;
};

/**
 * <!-- ************************************************************************************** -->
 * @struct PowerLawElectronShape
 * @brief Bounded dimensionless electron shape f(gamma) = gamma^(-p).
 * <!-- ************************************************************************************** -->
 */
struct PowerLawElectronShape {
    PowerLawElectronShape(Real gamma_min, Real gamma_max, Real p);

    Real gamma_min;
    Real gamma_max;
    Real p;

    [[nodiscard]] Real evaluate(Real gamma) const noexcept;
};

/**
 * <!-- ************************************************************************************** -->
 * @struct CutoffPowerLawElectronShape
 * @brief Bounded dimensionless electron shape f(gamma) = gamma^(-p) exp(-gamma/gamma_cut).
 * <!-- ************************************************************************************** -->
 */
struct CutoffPowerLawElectronShape {
    CutoffPowerLawElectronShape(Real gamma_min, Real gamma_max, Real p, Real gamma_cut);

    Real gamma_min;
    Real gamma_max;
    Real p;
    Real gamma_cut;

    [[nodiscard]] Real evaluate(Real gamma) const noexcept;
};

using NumericalElectronShape =
    std::variant<FlatElectronShape, PowerLawElectronShape, CutoffPowerLawElectronShape>;

/**
 * <!-- ************************************************************************************** -->
 * @brief Construct the logarithmic gamma grid used to sample a numerical electron shape.
 * @param gamma_min Lower Lorentz-factor bound, at least one
 * @param gamma_max Upper Lorentz-factor bound, greater than gamma_min
 * @param points_per_decade Approximate number of samples per gamma decade
 * @return Logarithmic grid with a number of intervals compatible with the existing quadrature
 * <!-- ************************************************************************************** -->
 */
[[nodiscard]] Array make_electron_gamma_grid(Real gamma_min, Real gamma_max, size_t points_per_decade = 32);

/**
 * <!-- ************************************************************************************** -->
 * @struct ElectronShapeTable
 * @brief Immutable logarithmic sampling of a bounded dimensionless electron shape.
 * @details Quadrature weights include dgamma. number_integral is integral f dgamma and energy_integral is
 *          integral (gamma-1) f dgamma. The table represents an instantaneous emitting population shape,
 *          not an injection function and not a cooled evolution of one.
 * <!-- ************************************************************************************** -->
 */
struct ElectronShapeTable {
    Array gamma;
    Array log2_gamma;
    Array shape;
    Array log2_shape;
    Array quadrature_weight;
    Array log2_quadrature_weight;
    Real number_integral{0};
    Real energy_integral{0};
    Real gamma_min{1};
    Real gamma_max{1};
    Real log2_gamma_min{0};
    Real log2_gamma_step{0};
};

/**
 * <!-- ************************************************************************************** -->
 * @brief Sample a bounded electron shape on a logarithmic gamma grid.
 * @param shape Built-in dimensionless electron shape
 * @param points_per_decade Approximate number of gamma nodes per decade
 * @return Shared immutable shape table
 * <!-- ************************************************************************************** -->
 */
[[nodiscard]] std::shared_ptr<ElectronShapeTable const>
sample_electron_shape(NumericalElectronShape const& shape, size_t points_per_decade = 32);

/**
 * <!-- ************************************************************************************** -->
 * @brief Finalize an electron shape that has already been sampled on a compatible logarithmic grid.
 * @details This is the common native-array boundary used by built-in shapes and one-shot external samplers.
 *          Values are copied into an immutable table and must be finite, non-negative, and not identically zero.
 * @param gamma Logarithmic gamma sampling grid from make_electron_gamma_grid
 * @param values Dimensionless unnormalized shape values evaluated on gamma
 * @return Shared immutable sampled shape table
 * <!-- ************************************************************************************** -->
 */
[[nodiscard]] std::shared_ptr<ElectronShapeTable const>
sample_electron_shape(Array const& gamma, Array const& values);

/**
 * <!-- ************************************************************************************** -->
 * @struct NumericalElectronDistribution
 * @brief Physical instantaneous electron column distribution dSigma_e/dgamma for one shock cell.
 * @details The supplied shape is authoritative. No gamma_m, gamma_c, absorption break, or cooling operation is
 *          applied. The normalization choice is mandatory because fixed shapes cannot generally satisfy both
 *          eps_e and xi_e.
 * <!-- ************************************************************************************** -->
 */
struct NumericalElectronDistribution {
    NumericalElectronDistribution(std::shared_ptr<ElectronShapeTable const> shape,
                                  ElectronNormalization normalization, Real proton_column, Real Gamma_th,
                                  Real eps_e, Real xi_e);

    std::shared_ptr<ElectronShapeTable const> shape;
    ElectronNormalization normalization;
    Real column_scale{0};       ///< A_cell in dSigma_e/dgamma = A_cell f(gamma)
    Real number_column{0};      ///< Integral dSigma_e/dgamma dgamma
    Real energy_column{0};      ///< Electron kinetic-energy column in internal energy/area units
    Real implied_xi_e{0};       ///< Electron number column divided by proton column
    Real implied_eps_e{0};      ///< Electron energy divided by the local shock thermal-energy column

    [[nodiscard]] Real gamma_min() const noexcept;
    [[nodiscard]] Real gamma_max() const noexcept;
    [[nodiscard]] Real compute_column_den(Real gamma) const noexcept;
};
