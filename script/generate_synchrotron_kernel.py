"""Generate the checked-in synchrotron-kernel lookup table.

This developer-only script requires mpmath. Production builds consume only the
generated numeric data and have no Python or special-function dependency.
"""

from __future__ import annotations

import math

import mpmath as mp


X_MIN = 1e-2
X_MAX = 50.0
POINTS_PER_DECADE = 48
GUARD_POINTS = 2


def synchrotron_kernel(x: float) -> float:
    """Evaluate x integral_x^infinity K_5/3(z) dz at high precision."""
    x_mp = mp.mpf(x)
    return float(x_mp * mp.quad(lambda z: mp.besselk(mp.mpf(5) / 3, z), [x_mp, mp.inf]))


def synchrotron_absorption_kernel(x: float) -> float:
    """Evaluate G(x) = x^2 K_5/3(x) at high precision."""
    x_mp = mp.mpf(x)
    return float(x_mp * x_mp * mp.besselk(mp.mpf(5) / 3, x_mp))


def main() -> None:
    mp.mp.dps = 50
    intervals = math.ceil(math.log10(X_MAX / X_MIN) * POINTS_PER_DECADE)
    log2_step = (math.log2(X_MAX) - math.log2(X_MIN)) / intervals
    log2_start = math.log2(X_MIN) - GUARD_POINTS * log2_step
    size = intervals + 2 * GUARD_POINTS + 1

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
    values = []
    for i in range(size):
        x = math.exp2(log2_start + i * log2_step)
        values.append(math.log2(synchrotron_absorption_kernel(x)))
    for i in range(0, size, 4):
        print("    " + ", ".join(f"{v:.17g}" for v in values[i : i + 4]) + ",")
    print("};")


if __name__ == "__main__":
    main()
