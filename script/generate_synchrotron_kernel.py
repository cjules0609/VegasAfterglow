"""Generate and validate the checked-in synchrotron-kernel lookup tables.

This developer-only script requires mpmath. Production builds consume only the
generated numeric data and have no Python or special-function dependency. Run
with ``--validate`` to check the absorption kernel and the finite-boundary SSA
identity independently of the production quadrature.
"""

from __future__ import annotations

import argparse
import math

import mpmath as mp


X_MIN = 1e-2
X_MAX = 50.0
POINTS_PER_DECADE = 48
GUARD_POINTS = 2
LOG2_G_LOW_LEADING = 0.51905756454496719
LOG2_PI_OVER_2 = 0.6514961294723187
INV_LN2 = 1.4426950408889634


def synchrotron_kernel(x: float) -> float:
    """Evaluate x integral_x^infinity K_5/3(z) dz at high precision."""
    x_mp = mp.mpf(x)
    return float(x_mp * mp.quad(lambda z: mp.besselk(mp.mpf(5) / 3, z), [x_mp, mp.inf]))


def synchrotron_absorption_kernel(x: float) -> float:
    """Evaluate G(x) = x^2 K_5/3(x) at high precision."""
    x_mp = mp.mpf(x)
    return float(x_mp * x_mp * mp.besselk(mp.mpf(5) / 3, x_mp))


def table_geometry() -> tuple[int, float, float, int]:
    intervals = math.ceil(math.log10(X_MAX / X_MIN) * POINTS_PER_DECADE)
    log2_step = (math.log2(X_MAX) - math.log2(X_MIN)) / intervals
    log2_start = math.log2(X_MIN) - GUARD_POINTS * log2_step
    size = intervals + 2 * GUARD_POINTS + 1
    return intervals, log2_step, log2_start, size


def absorption_table() -> list[float]:
    _, log2_step, log2_start, size = table_geometry()
    return [
        math.log2(synchrotron_absorption_kernel(math.exp2(log2_start + i * log2_step)))
        for i in range(size)
    ]


def interpolated_log2_absorption(log2_x: float, values: list[float]) -> float:
    """Mirror the production G(x) asymptotes and four-point log-space interpolation."""
    if log2_x < math.log2(X_MIN):
        x2 = math.exp2(2 * log2_x)
        return LOG2_G_LOW_LEADING + log2_x / 3 + math.log2(1 - 3 * x2 / 8)
    if log2_x > math.log2(X_MAX):
        x = math.exp2(log2_x)
        inv_x = 1 / x
        correction = 1 + 91 * inv_x / 72 + 1729 * inv_x * inv_x / 10368
        return 0.5 * LOG2_PI_OVER_2 + 1.5 * log2_x - x * INV_LN2 + math.log2(correction)

    _, step, start, _ = table_geometry()
    position = (log2_x - start) / step
    i = int(position)
    fraction = position - i
    y0, y1, y2, y3 = values[i - 1 : i + 3]
    return y1 + 0.5 * fraction * (
        y2
        - y0
        + fraction
        * (2 * y0 - 5 * y1 + 4 * y2 - y3 + fraction * (3 * (y1 - y2) + y3 - y0))
    )


def _logarithmic_samples(lower: float, upper: float, count: int) -> list[float]:
    return [math.exp(math.log(lower) + i * math.log(upper / lower) / (count - 1)) for i in range(count)]


def validate_absorption_kernel() -> tuple[float, float, str]:
    """Compare the complete production-style G(x) evaluator with mpmath."""
    values = absorption_table()
    regions = {
        "low-x asymptote": _logarithmic_samples(1e-20, X_MIN, 401),
        "lookup table": _logarithmic_samples(X_MIN, X_MAX, 2001),
        "high-x asymptote": _logarithmic_samples(X_MAX, 700, 401),
    }
    worst_error = -1.0
    worst_x = math.nan
    worst_region = ""
    for region, samples in regions.items():
        for x in samples:
            approximate_log2 = interpolated_log2_absorption(math.log2(x), values)
            x_mp = mp.mpf(x)
            reference_log2 = mp.log(x_mp * x_mp * mp.besselk(mp.mpf(5) / 3, x_mp), 2)
            relative_error = float(abs(mp.power(2, mp.mpf(approximate_log2) - reference_log2) - 1))
            if relative_error > worst_error:
                worst_error = relative_error
                worst_x = x
                worst_region = region
    if worst_error > 1e-4:
        raise RuntimeError(f"G(x) validation failed: maximum relative error {worst_error:.6e}")
    return worst_error, worst_x, worst_region


def _integrate_bessel_representation(function, x_min: mp.mpf) -> mp.mpf:
    """Integrate over the K_nu exponential representation on a finite partition."""
    upper = max(mp.mpf(12), mp.log(2 / x_min) + 6)
    points = [mp.mpf(0)]
    point = mp.mpf(1)
    while point < upper:
        points.append(point)
        point *= 2
    points.append(upper)
    return mp.quad(function, points)


def validate_finite_boundaries() -> mp.mpf:
    """Check the top-hat derivative form, including endpoints, against the G form."""
    cases = ((1e-4, 2, 20), (1, 10, 1e3), (1e4, 3, 50), (1e8, 100, 1e5), (1e12, 1e2, 1e8))
    order = mp.mpf(5) / 3
    worst_error = mp.mpf(0)
    for u, gamma_min, gamma_max in cases:
        x_a = mp.mpf(u) / mp.mpf(gamma_min) ** 2
        x_b = mp.mpf(u) / mp.mpf(gamma_max) ** 2

        def f_kernel(x: mp.mpf) -> mp.mpf:
            return x * _integrate_bessel_representation(
                lambda t: mp.exp(-x * mp.cosh(t)) * mp.cosh(order * t) / mp.cosh(t), x
            )

        # Transform the interior gamma integral to x, while evaluating F from
        # the independent exponential representation of K_5/3.
        interior = _integrate_bessel_representation(
            lambda t: mp.cosh(order * t)
            / mp.cosh(t) ** 2
            * (mp.exp(-x_b * mp.cosh(t)) - mp.exp(-x_a * mp.cosh(t))),
            x_b,
        )
        direct = interior - f_kernel(x_a) + f_kernel(x_b)

        # This is integral G(x)/x dx after the same gamma-to-x change of variables.
        by_parts = _integrate_bessel_representation(
            lambda t: mp.cosh(order * t)
            / mp.cosh(t) ** 2
            * (
                (x_b * mp.cosh(t) + 1) * mp.exp(-x_b * mp.cosh(t))
                - (x_a * mp.cosh(t) + 1) * mp.exp(-x_a * mp.cosh(t))
            ),
            x_b,
        )
        worst_error = max(worst_error, abs(direct / by_parts - 1))
    if worst_error > mp.mpf("1e-30"):
        raise RuntimeError(f"finite-boundary validation failed: relative error {mp.nstr(worst_error, 8)}")
    return worst_error


def validate() -> None:
    mp.mp.dps = 50
    kernel_error, worst_x, worst_region = validate_absorption_kernel()
    boundary_error = validate_finite_boundaries()
    _, step, _, size = table_geometry()
    print(f"table range: [{X_MIN:g}, {X_MAX:g}]")
    print(f"samples per decade: {POINTS_PER_DECADE}")
    print(f"stored samples: {size} ({GUARD_POINTS} guard points per side)")
    print(f"log2 spacing: {step:.17g}")
    print("interpolation: four-point cubic interpolation of log2(G) on a uniform log2(x) grid")
    print("validation range: [1e-20, 700]")
    print(f"maximum G relative error: {kernel_error:.8e}")
    print(f"worst case: x={worst_x:.8e} ({worst_region})")
    print(f"maximum finite-boundary identity relative error: {mp.nstr(boundary_error, 8)}")


def generate() -> None:
    mp.mp.dps = 50
    _, log2_step, log2_start, size = table_geometry()

    print(f"// size={size}, log2_start={log2_start:.17g}, log2_step={log2_step:.17g}")
    print("constexpr std::array<Real, %d> log2_kernel = {" % size)
    values = []
    for i in range(size):
        x = math.exp2(log2_start + i * log2_step)
        values.append(math.log2(synchrotron_kernel(x)))
    for i in range(0, size, 4):
        print("    " + ", ".join(f"{v:.17g}" for v in values[i : i + 4]) + ",")
    print("};")

    print("constexpr std::array<Real, %d> log2_absorption_kernel = {" % size)
    values = absorption_table()
    for i in range(0, size, 4):
        print("    " + ", ".join(f"{v:.17g}" for v in values[i : i + 4]) + ",")
    print("};")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validate", action="store_true", help="validate G(x) and the finite-boundary SSA identity")
    arguments = parser.parse_args()
    if arguments.validate:
        validate()
    else:
        generate()


if __name__ == "__main__":
    main()
