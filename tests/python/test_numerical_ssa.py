"""Numerical synchrotron self-absorption for instantaneous electron distributions."""

import numpy as np
import pytest

from VegasAfterglow import (
    CutoffPowerLawElectrons,
    ElectronDistribution,
    FlatElectrons,
    ISM,
    Model,
    Observer,
    PowerLawElectrons,
    Radiation,
    TophatJet,
)


def _model(electrons, *, ssa=None, reverse=False):
    kwargs = {} if ssa is None else {"ssa": ssa}
    radiation = Radiation(
        eps_e=0.1,
        eps_B=1e-2,
        p=2.3,
        xi_e=1,
        electrons=electrons,
        **kwargs,
    )
    return Model(
        jet=TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=100),
        medium=ISM(n_ism=10),
        observer=Observer(lumi_dist=1e28, z=0.1, theta_obs=0),
        fwd_rad=radiation,
        rvs_rad=radiation if reverse else None,
        resolutions=(0.1, 0.3, 8),
    )


def test_ssa_metadata_and_legacy_defaults():
    electrons = FlatElectrons(gamma_min=10, gamma_max=1e5, normalization="number")
    legacy = Radiation(eps_e=0.1, eps_B=1e-2, p=2.3, electrons=electrons)
    disabled = Radiation(eps_e=0.1, eps_B=1e-2, p=2.3, electrons=electrons, ssa=False)
    enabled = Radiation(eps_e=0.1, eps_B=1e-2, p=2.3, electrons=electrons, ssa=True)
    standard = Radiation(eps_e=0.1, eps_B=1e-2, p=2.3)

    assert legacy.supports_ssa and not legacy.ssa_enabled
    assert disabled.supports_ssa and not disabled.ssa_enabled
    assert enabled.supports_ssa and enabled.ssa_enabled
    assert standard.supports_ssa and standard.ssa_enabled
    with pytest.raises(ValueError, match="ssa=False is not supported with the standard"):
        Radiation(eps_e=0.1, eps_B=1e-2, p=2.3, ssa=False)

    time = np.logspace(3, 6, 10)
    frequency = np.full_like(time, 1e8)
    legacy_flux = _model(electrons).flux_density(time, frequency).total
    disabled_flux = _model(electrons, ssa=False).flux_density(time, frequency).total
    assert np.array_equal(legacy_flux, disabled_flux)


def test_optical_depth_diagnostic_and_transferred_spectrum():
    electrons = PowerLawElectrons(p=2.3, gamma_min=10, gamma_max=1e8, normalization="number")
    absorbed = _model(electrons, ssa=True).details(1e2, 1e6)
    thin = _model(electrons, ssa=False).details(1e2, 1e6)
    tau = absorbed.fwd.sync_optical_depth[0, 0, 2]
    absorbed_spectrum = absorbed.fwd.sync_spectrum[0, 0, 2]
    thin_spectrum = thin.fwd.sync_spectrum[0, 0, 2]
    frequency = np.logspace(5, 15, 100)

    optical_depth = np.asarray(tau(frequency))
    absorbed_intensity = np.asarray(absorbed_spectrum(frequency))
    thin_intensity = np.asarray(thin_spectrum(frequency))
    expected_transfer = -np.expm1(-optical_depth) / optical_depth
    assert np.all(np.isfinite(optical_depth))
    assert np.all(optical_depth >= 0)
    assert np.allclose(absorbed_intensity / thin_intensity, expected_transfer, rtol=3e-12, atol=0)
    assert absorbed_intensity[0] < 1e-5 * thin_intensity[0]
    assert absorbed_intensity[-1] == pytest.approx(thin_intensity[-1], rel=2e-12)
    assert np.all(np.isnan(absorbed.fwd.gamma_a))
    assert np.all(np.isnan(absorbed.fwd.nu_a))


@pytest.mark.parametrize(
    "electrons",
    [
        FlatElectrons(gamma_min=10, gamma_max=1e6, normalization="number"),
        PowerLawElectrons(p=2.3, gamma_min=10, gamma_max=1e8, normalization="number"),
        CutoffPowerLawElectrons(
            p=2.3, gamma_min=10, gamma_max=1e8, gamma_cut=1e5, normalization="number"
        ),
        ElectronDistribution(
            function=lambda gamma: np.exp(-0.5 * (np.log(gamma / 1e4) / 0.3) ** 2),
            gamma_min=10,
            gamma_max=1e7,
            normalization="number",
        ),
        ElectronDistribution(
            function=lambda gamma: np.where(gamma < 1e4, (gamma / 1e4) ** -1.5, (gamma / 1e4) ** -3),
            gamma_min=10,
            gamma_max=1e8,
            normalization="number",
        ),
        ElectronDistribution(
            function=lambda gamma: gamma**2 * np.exp(-gamma / 300),
            gamma_min=1,
            gamma_max=1e5,
            normalization="number",
        ),
        ElectronDistribution(
            function=lambda gamma: np.where((gamma > 1e3) & (gamma < 1e4), 0, gamma**-2),
            gamma_min=10,
            gamma_max=1e7,
            normalization="number",
        ),
    ],
    ids=["flat", "powerlaw", "cutoff", "pileup", "broken", "quasi_thermal", "internal_zero"],
)
def test_supported_shapes_have_finite_nonnegative_ssa(electrons):
    details = _model(electrons, ssa=True).details(1e2, 1e6)
    frequency = np.logspace(5, 16, 80)
    tau = np.asarray(details.fwd.sync_optical_depth[0, 0, 2](frequency))
    intensity = np.asarray(details.fwd.sync_spectrum[0, 0, 2](frequency))
    assert np.all(np.isfinite(tau))
    assert np.all(tau >= 0)
    assert np.all(np.isfinite(intensity))
    assert np.any(intensity > 0)


def test_numerical_ssa_runs_light_curve_reverse_shock_and_sky_image():
    electrons = FlatElectrons(gamma_min=10, gamma_max=1e5, normalization="number")
    model = _model(electrons, ssa=True, reverse=True)
    time = np.logspace(2, 5, 8)
    frequency = np.full_like(time, 1e8)
    flux = model.flux_density(time, frequency)
    assert np.any(np.asarray(flux.fwd.sync) > 0)
    assert np.any(np.asarray(flux.rvs.sync) > 0)

    details = model.details(1e2, 1e5)
    assert details.fwd.sync_optical_depth is not None
    assert details.rvs.sync_optical_depth is not None
    assert details.rvs.sync_optical_depth[0, 0, 2](1e8) >= 0

    image = model.sky_image(np.array([1e4]), nu_obs=1e8, fov=1e-8, npixel=16)
    assert image.image.shape == (1, 16, 16)
    assert np.all(np.isfinite(image.image))
    assert np.any(image.image > 0)


def test_custom_callable_is_not_revisited_when_lazy_ssa_is_constructed():
    calls = 0

    def shape(gamma):
        nonlocal calls
        calls += 1
        return gamma**-2.3

    electrons = ElectronDistribution(
        function=shape,
        gamma_min=1e2,
        gamma_max=1e8,
        normalization="energy",
    )
    assert calls == 1

    thin_model = _model(electrons, ssa=False)
    thin_model.flux_density(np.array([1e5]), np.array([1e10]))
    assert calls == 1

    absorbed_model = _model(electrons, ssa=True, reverse=True)
    absorbed_model.flux_density(np.array([1e5]), np.array([1e10]))
    details = absorbed_model.details(1e2, 1e5)
    details.fwd.sync_optical_depth[0, 0, 2](1e10)
    absorbed_model.sky_image(np.array([1e4]), nu_obs=1e8, fov=1e-8, npixel=8)
    assert calls == 1
