"""Runtime integration tests for bounded instantaneous electron populations."""

import numpy as np
import pytest

from VegasAfterglow import (
    CutoffPowerLawElectrons,
    FlatElectrons,
    ISM,
    Model,
    Observer,
    PowerLawElectrons,
    Radiation,
    TophatJet,
)


def _radiation(electrons=None, *, eps_e=0.1, xi_e=1.0):
    return Radiation(eps_e=eps_e, eps_B=1e-3, p=2.3, xi_e=xi_e, electrons=electrons)


def _model(radiation, rvs_rad=None):
    return Model(
        jet=TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=100),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=1e28, z=0.1, theta_obs=0.0),
        fwd_rad=radiation,
        rvs_rad=rvs_rad,
        resolutions=(0.1, 0.3, 8),
    )


def _cell_spectrum(model, electrons, t_min=1e3, t_max=1e5):
    details = model.details(t_min, t_max)
    grid = details.fwd.sync_spectrum
    assert grid is not None
    return details, grid[0, 0, 2]


def test_typed_electron_construction_and_radiation_properties():
    flat = FlatElectrons(gamma_min=1e2, gamma_max=1e7, normalization="energy")
    powerlaw = PowerLawElectrons(p=2.3, gamma_min=1e2, gamma_max=1e8, normalization="number")
    cutoff = CutoffPowerLawElectrons(
        p=2.3, gamma_min=1e2, gamma_max=1e8, gamma_cut=1e6, normalization="energy"
    )

    for electrons, name in ((flat, "flat"), (powerlaw, "powerlaw"), (cutoff, "cutoff_powerlaw")):
        rad = _radiation(electrons)
        assert rad.electron_model == name
        assert rad.electrons is not None
        assert rad.supports_ssa
        assert not rad.ssa_enabled
        assert not rad.supports_ssc

    standard = _radiation()
    assert standard.electron_model == "standard"
    assert standard.electrons is None
    assert standard.supports_ssa
    assert standard.ssa_enabled
    assert standard.supports_ssc


def test_default_and_explicit_none_are_bitwise_identical():
    historical = Radiation(eps_e=0.1, eps_B=1e-3, p=2.3)
    explicit_none = Radiation(eps_e=0.1, eps_B=1e-3, p=2.3, electrons=None)
    time = np.logspace(3, 6, 12)
    frequency = np.full_like(time, 1e14)

    old_flux = _model(historical).flux_density(time, frequency)
    none_flux = _model(explicit_none).flux_density(time, frequency)
    assert np.array_equal(old_flux.total, none_flux.total)
    assert np.array_equal(old_flux.fwd.sync, none_flux.fwd.sync)


def test_flat_distribution_produces_finite_light_curve_without_standard_breaks():
    electrons = FlatElectrons(gamma_min=1e2, gamma_max=1e7, normalization="energy")
    model = _model(_radiation(electrons))
    time = np.logspace(3, 6, 12)
    flux = np.asarray(model.flux_density(time, np.full_like(time, 1e14)).total)
    assert np.all(np.isfinite(flux))
    assert np.all(flux > 0)

    details, spectrum = _cell_spectrum(model, electrons)
    assert np.all(np.isnan(details.fwd.gamma_m))
    assert np.all(np.isnan(details.fwd.gamma_c))
    assert spectrum(1e14) > 0


def test_numerical_powerlaw_has_optically_thin_asymptotic_slope():
    p = 2.3
    electrons = PowerLawElectrons(p=p, gamma_min=1e2, gamma_max=1e8, normalization="energy")
    details, spectrum = _cell_spectrum(_model(_radiation(electrons)), electrons)
    B = float(details.fwd.B_comv[0, 0, 2])
    frequencies = np.logspace(np.log10(4e6 * B * 1e6), np.log10(4e6 * B * 1e12), 80)
    intensity = np.asarray(spectrum(frequencies))
    slope = np.polyfit(np.log(frequencies), np.log(intensity), 1)[0]
    assert slope == pytest.approx(-(p - 1) / 2, abs=0.03)


def test_cutoff_powerlaw_suppresses_high_frequency_emission():
    common = dict(p=2.3, gamma_min=1e2, gamma_max=1e8, normalization="energy")
    powerlaw = PowerLawElectrons(**common)
    cutoff = CutoffPowerLawElectrons(**common, gamma_cut=1e5)
    power_details, power_spectrum = _cell_spectrum(_model(_radiation(powerlaw)), powerlaw)
    _, cutoff_spectrum = _cell_spectrum(_model(_radiation(cutoff)), cutoff)
    B = float(power_details.fwd.B_comv[0, 0, 2])
    low = 4e6 * B * 1e6
    high = 4e6 * B * 1e12
    power_ratio = power_spectrum(high) / power_spectrum(low)
    cutoff_ratio = cutoff_spectrum(high) / cutoff_spectrum(low)
    assert cutoff_ratio < 0.25 * power_ratio


@pytest.mark.parametrize("normalization", ["energy", "number"])
def test_normalization_modes_and_diagnostics(normalization):
    electrons = FlatElectrons(gamma_min=1e2, gamma_max=1e5, normalization=normalization)
    eps_e = 0.07
    xi_e = 0.2
    details = _model(_radiation(electrons, eps_e=eps_e, xi_e=xi_e)).details(1e3, 1e5)
    shock = details.fwd

    assert np.all(shock.electron_gamma_min == 1e2)
    assert np.all(shock.electron_gamma_max == 1e5)
    if normalization == "energy":
        assert np.allclose(shock.electron_energy_fraction, eps_e, rtol=1e-12, atol=0)
    else:
        assert np.allclose(shock.electron_number_fraction, xi_e, rtol=1e-12, atol=0)

    distribution = shock.electron_column_distribution[0, 0, 2]
    gamma = np.array([10.0, 1e2, 1e4, 1e5, 1e6])
    values = np.asarray(distribution(gamma))
    assert values[0] == 0
    assert values[-1] == 0
    assert np.all(values[1:-1] > 0)
    assert distribution(1e4) == pytest.approx(values[2])


def test_numerical_electrons_support_reverse_shock_and_sky_image():
    electrons = FlatElectrons(gamma_min=10, gamma_max=1e6, normalization="number")
    rad = _radiation(electrons, xi_e=0.3)
    model = _model(rad, rvs_rad=rad)
    time = np.logspace(2, 5, 8)
    flux = model.flux_density(time, np.full_like(time, 1e12))
    assert np.any(np.asarray(flux.rvs.sync) > 0)

    details = model.details(1e2, 1e5)
    assert details.rvs.electron_column_distribution is not None
    assert details.rvs.sync_spectrum is not None

    image = model.sky_image(np.array([1e4]), nu_obs=1e12, fov=1e-8, npixel=16)
    assert image.image.shape == (1, 16, 16)
    assert np.all(np.isfinite(image.image))
    assert np.any(image.image > 0)


@pytest.mark.parametrize(
    "factory",
    [
        lambda: FlatElectrons(gamma_min=0.9, gamma_max=10, normalization="energy"),
        lambda: FlatElectrons(gamma_min=10, gamma_max=10, normalization="energy"),
        lambda: FlatElectrons(gamma_min=np.nan, gamma_max=10, normalization="energy"),
        lambda: PowerLawElectrons(p=np.inf, gamma_min=1, gamma_max=10, normalization="number"),
        lambda: CutoffPowerLawElectrons(
            p=2.3, gamma_min=1, gamma_max=10, gamma_cut=0, normalization="energy"
        ),
        lambda: FlatElectrons(gamma_min=1, gamma_max=10, normalization="automatic"),
    ],
)
def test_invalid_electron_configuration(factory):
    with pytest.raises((ValueError, RuntimeError)):
        factory()


def test_unsupported_radiation_combinations_fail_clearly():
    electrons = FlatElectrons(gamma_min=1, gamma_max=10, normalization="energy")
    with pytest.raises(ValueError, match="ssc=True is not currently supported"):
        Radiation(eps_e=0.1, eps_B=1e-3, p=2.3, ssc=True, electrons=electrons)
    with pytest.raises(ValueError, match="kn=True is not currently supported"):
        Radiation(eps_e=0.1, eps_B=1e-3, p=2.3, kn=True, electrons=electrons)
    with pytest.raises(TypeError, match="electrons must be None"):
        Radiation(eps_e=0.1, eps_B=1e-3, p=2.3, electrons=object())
