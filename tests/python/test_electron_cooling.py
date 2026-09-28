"""Validated one-zone cooling for arbitrary instantaneous electron shapes."""

import numpy as np

from VegasAfterglow import ElectronCooling, ElectronDistribution, Radiation


def _narrow_population(gamma0=1e4):
    return ElectronDistribution(
        function=lambda gamma: np.exp(-0.5 * (np.log(gamma / gamma0) / 0.015) ** 2),
        gamma_min=1.0,
        gamma_max=1e6,
        normalization="number",
        samples_per_decade=128,
    )


def test_synchrotron_cooling_conserves_number_and_reduces_energy():
    electrons = _narrow_population()
    cooling = ElectronCooling(synchrotron=True, max_step_fraction=0.2)
    result = cooling.evolve(electrons, np.array([0.0, 2e4, 4e4]), 10.0)

    assert result.distribution.shape == (3, result.gamma.size)
    np.testing.assert_allclose(result.number, result.number[0], rtol=2e-12)
    assert np.all(np.diff(result.kinetic_energy) < 0)
    assert result.escaped_lower == 0


def test_continuous_injection_and_final_distribution_round_trip():
    electrons = ElectronDistribution(
        function=lambda gamma: gamma**-2.3,
        gamma_min=1.0,
        gamma_max=1e5,
        normalization="number",
        samples_per_decade=32,
    )
    cooling = ElectronCooling(synchrotron=False)
    source = np.ones_like(electrons.gamma) * 1e-5
    result = cooling.evolve(
        electrons,
        np.array([0.0, 2.0, 5.0]),
        magnetic_field=0.0,
        injection=source,
    )

    assert result.number[-1] > result.number[0]
    final = result.final_electrons()
    assert final.normalization == "number"
    assert np.all(np.isfinite(final.shape))
    assert np.all(final.shape >= 0)
    Radiation(eps_e=0.1, eps_B=1e-3, p=2.3, electrons=final)


def test_optional_loss_channels_are_independent():
    electrons = _narrow_population(1e3)
    times = np.array([0.0, 1e3])
    unchanged = ElectronCooling(synchrotron=False).evolve(electrons, times, 10.0)
    adiabatic = ElectronCooling(synchrotron=False, adiabatic=True).evolve(
        electrons, times, 10.0, expansion_rate=1e-3
    )
    ic = ElectronCooling(synchrotron=False, inverse_compton=True).evolve(
        electrons, times, 0.0, photon_energy_density=1.0
    )

    np.testing.assert_allclose(unchanged.kinetic_energy, unchanged.kinetic_energy[0])
    assert adiabatic.kinetic_energy[-1] < adiabatic.kinetic_energy[0]
    assert ic.kinetic_energy[-1] < ic.kinetic_energy[0]
