"""User-defined instantaneous electron distributions sampled once at construction."""

import gc
import weakref

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


def _model(electrons, *, p=2.3, radiative_fireball=False):
    return Model(
        jet=TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=100),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=1e28, z=0.1, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=0.1, eps_B=1e-3, p=p, electrons=electrons),
        resolutions=(0.1, 0.3, 8),
        radiative_fireball=radiative_fireball,
    )


def _spectrum(electrons, *, p=2.3):
    details = _model(electrons, p=p).details(1e3, 1e5)
    return details, details.fwd.sync_spectrum[0, 0, 2]


@pytest.mark.parametrize(
    ("custom_function", "built_in"),
    [
        (
            lambda gamma: np.ones_like(gamma),
            FlatElectrons(gamma_min=1e2, gamma_max=1e8, normalization="energy"),
        ),
        (
            lambda gamma: gamma**-2.3,
            PowerLawElectrons(p=2.3, gamma_min=1e2, gamma_max=1e8, normalization="energy"),
        ),
        (
            lambda gamma: gamma**-2.3 * np.exp(-gamma / 1e6),
            CutoffPowerLawElectrons(
                p=2.3, gamma_min=1e2, gamma_max=1e8, gamma_cut=1e6, normalization="energy"
            ),
        ),
    ],
)
def test_custom_shapes_match_built_in_models(custom_function, built_in):
    custom = ElectronDistribution(
        function=custom_function,
        gamma_min=1e2,
        gamma_max=1e8,
        normalization="energy",
    )
    custom_details, custom_spectrum = _spectrum(custom)
    _, built_in_spectrum = _spectrum(built_in)
    B = float(custom_details.fwd.B_comv[0, 0, 2])
    frequency = np.logspace(np.log10(4e6 * B * 1e4), np.log10(4e6 * B * 1e16), 120)
    assert np.allclose(custom_spectrum(frequency), built_in_spectrum(frequency), rtol=2e-13, atol=0)


def test_custom_flat_is_exactly_equivalent_to_built_in_flat():
    custom = ElectronDistribution(
        function=lambda gamma: np.ones_like(gamma),
        gamma_min=1e2,
        gamma_max=1e7,
        normalization="number",
    )
    built_in = FlatElectrons(gamma_min=1e2, gamma_max=1e7, normalization="number")
    time = np.logspace(3, 6, 10)
    frequency = np.full_like(time, 1e14)
    custom_flux = _model(custom).flux_density(time, frequency).total
    built_in_flux = _model(built_in).flux_density(time, frequency).total
    assert np.array_equal(custom_flux, built_in_flux)


def test_callable_is_vectorized_called_once_and_not_retained():
    class Shape:
        def __init__(self):
            self.calls = 0
            self.last_shape = None

        def __call__(self, gamma):
            self.calls += 1
            self.last_shape = gamma.shape
            return np.exp(-0.5 * (np.log(gamma / 1e5) / 0.4) ** 2)

    function = Shape()
    function_ref = weakref.ref(function)
    electrons = ElectronDistribution(
        function=function,
        gamma_min=1e2,
        gamma_max=1e8,
        normalization="energy",
    )
    assert function.calls == 1
    assert function.last_shape == electrons.gamma.shape

    del function
    gc.collect()
    assert function_ref() is None

    model = _model(electrons)
    time = np.logspace(3, 5, 5)
    assert np.all(model.flux_density(time, np.full_like(time, 1e14)).total > 0)
    assert model.details(1e3, 1e5).fwd.sync_spectrum is not None
    assert np.any(model.sky_image(np.array([1e4]), 1e14, fov=1e-8, npixel=16).image > 0)


def test_sampled_shape_properties_are_copies():
    electrons = ElectronDistribution(
        function=lambda gamma: gamma**-2,
        gamma_min=10,
        gamma_max=1e6,
        normalization="number",
        samples_per_decade=48,
    )
    assert electrons.electron_model == "custom"
    assert electrons.gamma_min == 10
    assert electrons.gamma_max == 1e6
    assert electrons.normalization == "number"
    assert electrons.samples_per_decade == 48
    gamma = electrons.gamma
    shape = electrons.shape
    gamma[:] = -1
    shape[:] = -1
    assert np.all(electrons.gamma >= 1)
    assert np.all(electrons.shape >= 0)


def test_pileup_emission_peaks_near_its_characteristic_frequency():
    gamma0 = 1e5
    sigma = 0.18
    electrons = ElectronDistribution(
        function=lambda gamma: np.exp(-0.5 * (np.log(gamma / gamma0) / sigma) ** 2),
        gamma_min=1e2,
        gamma_max=1e8,
        normalization="energy",
        samples_per_decade=64,
    )
    details, spectrum = _spectrum(electrons)
    B = float(details.fwd.B_comv[0, 0, 2])
    expected = 4e6 * B * gamma0**2
    frequency = np.logspace(np.log10(expected) - 2, np.log10(expected) + 2, 300)
    peak = frequency[np.argmax(spectrum(frequency))]
    assert abs(np.log10(peak / expected)) < 0.7


def test_broken_powerlaw_photon_spectrum_tracks_the_custom_break():
    gamma_break = 1e5
    p_low = 2.0
    p_high = 3.4

    def broken_powerlaw(gamma):
        ratio = gamma / gamma_break
        return np.where(gamma <= gamma_break, ratio**-p_low, ratio**-p_high)

    electrons = ElectronDistribution(
        function=broken_powerlaw,
        gamma_min=1e2,
        gamma_max=1e8,
        normalization="energy",
        samples_per_decade=64,
    )
    details, spectrum = _spectrum(electrons)
    B = float(details.fwd.B_comv[0, 0, 2])

    low_frequency = np.logspace(np.log10(4e6 * B * 1e6), np.log10(4e6 * B * 1e8), 50)
    high_frequency = np.logspace(np.log10(4e6 * B * 1e12), np.log10(4e6 * B * 1e14), 50)
    low_slope = np.polyfit(np.log(low_frequency), np.log(spectrum(low_frequency)), 1)[0]
    high_slope = np.polyfit(np.log(high_frequency), np.log(spectrum(high_frequency)), 1)[0]
    assert low_slope == pytest.approx(-(p_low - 1) / 2, abs=0.08)
    assert high_slope == pytest.approx(-(p_high - 1) / 2, abs=0.08)
    assert high_slope < low_slope - 0.45


def test_custom_shape_controls_spectrum_and_radiation_p_does_not_reshape_it():
    bounds = dict(gamma_min=1e2, gamma_max=1e8, normalization="energy")
    flat = ElectronDistribution(function=lambda gamma: np.ones_like(gamma), **bounds)
    steep = ElectronDistribution(function=lambda gamma: gamma**-3, **bounds)
    flat_details, flat_spectrum = _spectrum(flat, p=2.1)
    _, steep_spectrum = _spectrum(steep, p=2.1)
    _, same_flat_other_p = _spectrum(flat, p=2.9)
    B = float(flat_details.fwd.B_comv[0, 0, 2])
    frequency = np.logspace(np.log10(4e6 * B * 1e6), np.log10(4e6 * B * 1e12), 80)

    flat_values = flat_spectrum(frequency)
    steep_values = steep_spectrum(frequency)
    assert not np.allclose(flat_values / flat_values.max(), steep_values / steep_values.max(), rtol=1e-3)
    assert np.allclose(flat_values, same_flat_other_p(frequency), rtol=2e-13, atol=0)


def test_custom_distribution_is_not_reshaped_by_cooling_age():
    electrons = ElectronDistribution(
        function=lambda gamma: np.where(gamma < 1e5, gamma**-1.5, 1e2 * gamma**-2.5),
        gamma_min=1e2,
        gamma_max=1e8,
        normalization="number",
    )
    details = _model(electrons).details(1e2, 1e7)
    gamma = np.logspace(2.2, 7.8, 100)
    early = details.fwd.electron_column_distribution[0, 0, 2](gamma)
    late_index = details.fwd.Gamma.shape[2] - 3
    late = details.fwd.electron_column_distribution[0, 0, late_index](gamma)
    early /= early[20]
    late /= late[20]
    assert np.allclose(early, late, rtol=2e-13, atol=0)
    assert np.all(np.isnan(details.fwd.gamma_c))


@pytest.mark.parametrize(
    ("function", "match"),
    [
        (lambda gamma: 1.0, "one-dimensional array"),
        (lambda gamma: np.ones(gamma.size - 1), "matching the gamma sampling grid"),
        (lambda gamma: np.ones((gamma.size, 1)), "one-dimensional array"),
        (lambda gamma: np.full_like(gamma, np.nan), "finite and non-negative"),
        (lambda gamma: np.full_like(gamma, np.inf), "finite and non-negative"),
        (lambda gamma: -np.ones_like(gamma), "finite and non-negative"),
        (lambda gamma: np.ones_like(gamma, dtype=complex), "real values"),
        (lambda gamma: np.zeros_like(gamma), "positive number integral"),
    ],
)
def test_invalid_callable_outputs(function, match):
    with pytest.raises((TypeError, ValueError), match=match):
        ElectronDistribution(
            function=function,
            gamma_min=1,
            gamma_max=1e4,
            normalization="energy",
        )


def test_callable_exception_preserves_cause_and_context():
    def fails(_gamma):
        raise LookupError("deliberate shape failure")

    with pytest.raises(RuntimeError, match="failed while sampling electron distribution") as error:
        ElectronDistribution(function=fails, gamma_min=1, gamma_max=10, normalization="energy")
    assert isinstance(error.value.__cause__, LookupError)
    assert "deliberate shape failure" in str(error.value.__cause__)


@pytest.mark.parametrize("samples_per_decade", [0, 3, 513, 10000])
def test_invalid_sampling_resolution(samples_per_decade):
    with pytest.raises(ValueError, match="samples_per_decade must be between 4 and 512"):
        ElectronDistribution(
            function=lambda gamma: np.ones_like(gamma),
            gamma_min=1,
            gamma_max=10,
            normalization="energy",
            samples_per_decade=samples_per_decade,
        )


@pytest.mark.parametrize(
    ("gamma_min", "gamma_max", "normalization"),
    [
        (0.9, 10, "energy"),
        (10, 10, "energy"),
        (np.nan, 10, "number"),
        (1, np.inf, "number"),
        (1, 10, "automatic"),
    ],
)
def test_invalid_custom_bounds_and_normalization(gamma_min, gamma_max, normalization):
    with pytest.raises((ValueError, RuntimeError)):
        ElectronDistribution(
            function=lambda gamma: np.ones_like(gamma),
            gamma_min=gamma_min,
            gamma_max=gamma_max,
            normalization=normalization,
        )


def test_noncallable_and_unsupported_ssc_are_rejected():
    with pytest.raises(TypeError, match="function must be callable"):
        ElectronDistribution(function=3.0, gamma_min=1, gamma_max=10, normalization="energy")

    electrons = ElectronDistribution(
        function=lambda gamma: np.ones_like(gamma), gamma_min=1, gamma_max=10, normalization="energy"
    )
    with pytest.raises(ValueError, match="ssc=True is not currently supported"):
        Radiation(eps_e=0.1, eps_B=1e-3, p=2.3, ssc=True, electrons=electrons)
    with pytest.raises(ValueError, match="kn=True is not currently supported"):
        Radiation(eps_e=0.1, eps_B=1e-3, p=2.3, kn=True, electrons=electrons)
