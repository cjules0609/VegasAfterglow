//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "electron-distribution.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "../core/quadrature.h"
#include "../util/fast-math.h"
#include "../util/macros.h"

namespace {
    void validate_bounds(Real gamma_min, Real gamma_max) {
        if (!(std::isfinite(gamma_min) && gamma_min >= 1)) {
            throw std::invalid_argument("gamma_min must be finite and >= 1");
        }
        if (!(std::isfinite(gamma_max) && gamma_max > gamma_min)) {
            throw std::invalid_argument("gamma_max must be finite and greater than gamma_min");
        }
    }

    [[nodiscard]] std::pair<Real, Real> shape_bounds(NumericalElectronShape const& shape) {
        return std::visit([](auto const& model) { return std::pair{model.gamma_min, model.gamma_max}; }, shape);
    }
} // namespace

FlatElectronShape::FlatElectronShape(Real gamma_min, Real gamma_max) : gamma_min(gamma_min), gamma_max(gamma_max) {
    validate_bounds(gamma_min, gamma_max);
}

Real FlatElectronShape::evaluate(Real gamma) const noexcept {
    return (gamma >= gamma_min && gamma <= gamma_max) ? 1 : 0;
}

PowerLawElectronShape::PowerLawElectronShape(Real gamma_min, Real gamma_max, Real p)
    : gamma_min(gamma_min), gamma_max(gamma_max), p(p) {
    validate_bounds(gamma_min, gamma_max);
    if (!std::isfinite(p)) {
        throw std::invalid_argument("p must be finite");
    }
}

Real PowerLawElectronShape::evaluate(Real gamma) const noexcept {
    if (!(gamma >= gamma_min && gamma <= gamma_max)) {
        return 0;
    }
    return std::exp(-p * std::log(gamma));
}

CutoffPowerLawElectronShape::CutoffPowerLawElectronShape(Real gamma_min, Real gamma_max, Real p, Real gamma_cut)
    : gamma_min(gamma_min), gamma_max(gamma_max), p(p), gamma_cut(gamma_cut) {
    validate_bounds(gamma_min, gamma_max);
    if (!std::isfinite(p)) {
        throw std::invalid_argument("p must be finite");
    }
    if (!(std::isfinite(gamma_cut) && gamma_cut > 0)) {
        throw std::invalid_argument("gamma_cut must be finite and greater than zero");
    }
}

Real CutoffPowerLawElectronShape::evaluate(Real gamma) const noexcept {
    if (!(gamma >= gamma_min && gamma <= gamma_max)) {
        return 0;
    }
    return std::exp(-p * std::log(gamma) - gamma / gamma_cut);
}

Array make_electron_gamma_grid(Real gamma_min, Real gamma_max, size_t points_per_decade) {
    validate_bounds(gamma_min, gamma_max);
    if (points_per_decade == 0) {
        throw std::invalid_argument("points_per_decade must be greater than zero");
    }

    const Real decades = std::log10(gamma_max / gamma_min);
    size_t intervals = std::max<size_t>(6, static_cast<size_t>(std::ceil(decades * points_per_decade)));
    intervals += (6 - intervals % 6) % 6; // align with the 7-point Newton-Cotes groups
    const size_t size = intervals + 1;

    Array gamma = Array::from_shape({size});
    const Real log2_gamma_min = fast_log2(gamma_min);
    const Real log2_gamma_step = fast_log2(gamma_max / gamma_min) / static_cast<Real>(intervals);
    for (size_t i = 0; i < size; ++i) {
        gamma(i) = fast_exp2(log2_gamma_min + log2_gamma_step * static_cast<Real>(i));
    }
    gamma(0) = gamma_min;
    gamma(size - 1) = gamma_max;
    return gamma;
}

std::shared_ptr<ElectronShapeTable const> sample_electron_shape(Array const& gamma, Array const& values) {
    if (gamma.dimension() != 1 || values.dimension() != 1 || gamma.size() != values.size()) {
        throw std::invalid_argument("electron gamma and shape samples must be one-dimensional and equal length");
    }
    if (gamma.size() < 7 || (gamma.size() - 1) % 6 != 0) {
        throw std::invalid_argument("electron gamma grid must contain 6n + 1 samples for quadrature");
    }

    const size_t size = gamma.size();
    const Real gamma_min = gamma(0);
    const Real gamma_max = gamma(size - 1);
    validate_bounds(gamma_min, gamma_max);

    auto table = std::make_shared<ElectronShapeTable>();
    table->gamma = gamma;
    table->log2_gamma = Array::from_shape({size});
    table->shape = values;
    table->log2_shape = Array::from_shape({size});
    table->gamma_min = gamma_min;
    table->gamma_max = gamma_max;
    table->log2_gamma_min = fast_log2(gamma_min);
    table->log2_gamma_step = fast_log2(gamma_max / gamma_min) / static_cast<Real>(size - 1);

    for (size_t i = 0; i < size; ++i) {
        const Real expected = fast_exp2(table->log2_gamma_min + table->log2_gamma_step * static_cast<Real>(i));
        const Real relative_error = std::abs(table->gamma(i) / expected - 1);
        const bool log_uniform = (i == 0 || i + 1 == size || relative_error < 1e-10);
        if (!(std::isfinite(table->gamma(i)) && table->gamma(i) >= 1 && log_uniform)) {
            throw std::invalid_argument("electron gamma samples must form a finite uniform logarithmic grid");
        }
        table->log2_gamma(i) = fast_log2(table->gamma(i));
    }

    compute_nc7_weights(table->gamma, table->quadrature_weight);
    table->log2_quadrature_weight = Array::from_shape({size});
    Real number_integral = 0;
    Real energy_integral = 0;
    Real number_compensation = 0;
    Real energy_compensation = 0;
    for (size_t i = 0; i < size; ++i) {
        const Real value = table->shape(i);
        if (!(std::isfinite(value) && value >= 0)) {
            throw std::invalid_argument("electron shape must be finite and non-negative throughout its bounds");
        }
        table->log2_shape(i) = (value > 0) ? fast_log2(value) : -con::inf;
        table->log2_quadrature_weight(i) = fast_log2(table->quadrature_weight(i));

        // Kahan summation: both moments are positive, but may span many decades.
        const Real number_term = value * table->quadrature_weight(i);
        const Real number_y = number_term - number_compensation;
        const Real number_next = number_integral + number_y;
        number_compensation = (number_next - number_integral) - number_y;
        number_integral = number_next;

        const Real energy_term = (table->gamma(i) - 1) * number_term;
        const Real energy_y = energy_term - energy_compensation;
        const Real energy_next = energy_integral + energy_y;
        energy_compensation = (energy_next - energy_integral) - energy_y;
        energy_integral = energy_next;
    }
    if (!(std::isfinite(number_integral) && number_integral > 0)) {
        throw std::invalid_argument("electron shape has no finite positive number integral");
    }
    if (!(std::isfinite(energy_integral) && energy_integral > 0)) {
        throw std::invalid_argument("electron shape has no finite positive energy integral");
    }
    table->number_integral = number_integral;
    table->energy_integral = energy_integral;
    return table;
}

std::shared_ptr<ElectronShapeTable const> sample_electron_shape(NumericalElectronShape const& model,
                                                                size_t points_per_decade) {
    const auto [gamma_min, gamma_max] = shape_bounds(model);
    Array gamma = make_electron_gamma_grid(gamma_min, gamma_max, points_per_decade);
    Array values = Array::from_shape(gamma.shape());
    std::visit(
        [&](auto const& shape) {
            for (size_t i = 0; i < gamma.size(); ++i) {
                values(i) = shape.evaluate(gamma(i));
            }
        },
        model);
    return sample_electron_shape(gamma, values);
}

NumericalElectronDistribution::NumericalElectronDistribution(std::shared_ptr<ElectronShapeTable const> shape,
                                                             ElectronNormalization normalization, Real proton_column,
                                                             Real Gamma_th, Real eps_e, Real xi_e)
    : shape(std::move(shape)), normalization(normalization) {
    if (!this->shape) {
        throw std::invalid_argument("electron shape table must not be null");
    }
    if (!(std::isfinite(proton_column) && proton_column >= 0)) {
        throw std::invalid_argument("proton_column must be finite and non-negative");
    }
    if (!(std::isfinite(Gamma_th) && Gamma_th >= 1)) {
        throw std::invalid_argument("Gamma_th must be finite and >= 1");
    }
    if (!(std::isfinite(eps_e) && eps_e >= 0)) {
        throw std::invalid_argument("eps_e must be finite and non-negative");
    }
    if (!(std::isfinite(xi_e) && xi_e >= 0)) {
        throw std::invalid_argument("xi_e must be finite and non-negative");
    }

    const Real thermal_gamma = Gamma_th - 1;
    if (normalization == ElectronNormalization::energy) {
        column_scale = eps_e * thermal_gamma * (con::mp / con::me) * proton_column / this->shape->energy_integral;
    } else if (normalization == ElectronNormalization::number) {
        column_scale = xi_e * proton_column / this->shape->number_integral;
    } else {
        throw std::invalid_argument("unknown electron normalization mode");
    }
    if (!(std::isfinite(column_scale) && column_scale >= 0)) {
        throw std::invalid_argument("electron normalization is not finite");
    }

    number_column = column_scale * this->shape->number_integral;
    energy_column = con::mec2 * column_scale * this->shape->energy_integral;
    implied_xi_e = (proton_column > 0) ? number_column / proton_column : 0;

    const Real thermal_energy_column = thermal_gamma * con::mpc2 * proton_column;
    if (thermal_energy_column > 0) {
        implied_eps_e = energy_column / thermal_energy_column;
    } else {
        implied_eps_e = (energy_column > 0) ? con::inf : 0;
    }
}

Real NumericalElectronDistribution::gamma_min() const noexcept {
    return shape ? shape->gamma_min : 1;
}

Real NumericalElectronDistribution::gamma_max() const noexcept {
    return shape ? shape->gamma_max : 1;
}

Real NumericalElectronDistribution::compute_column_den(Real gamma) const noexcept {
    if (!shape || !(gamma >= shape->gamma_min && gamma <= shape->gamma_max) || column_scale == 0) {
        return 0;
    }
    const size_t size = shape->gamma.size();
    if (gamma == shape->gamma_max) {
        return column_scale * shape->shape(size - 1);
    }

    const Real position = (fast_log2(gamma) - shape->log2_gamma_min) / shape->log2_gamma_step;
    const size_t i = std::min(static_cast<size_t>(std::max(position, 0.0)), size - 2);
    const Real f = position - static_cast<Real>(i);
    const Real lo = shape->log2_shape(i);
    const Real hi = shape->log2_shape(i + 1);
    if (std::isfinite(lo) && std::isfinite(hi)) {
        return column_scale * fast_exp2(lo + f * (hi - lo));
    }
    return column_scale * ((1 - f) * shape->shape(i) + f * shape->shape(i + 1));
}
