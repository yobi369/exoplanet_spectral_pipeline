import numpy as np
import pytest

from exospec import (
    SpectrumInputError,
    analyze_spectrum,
    flag_outliers,
    illustrative_templates,
    load_spectrum_csv,
    make_synthetic_spectrum,
)


def test_synthetic_data_is_deterministic():
    first = make_synthetic_spectrum(seed=12)
    second = make_synthetic_spectrum(seed=12)
    for left, right in zip(first, second, strict=True):
        np.testing.assert_array_equal(left, right)


def test_rejects_bad_wavelength_order():
    with pytest.raises(SpectrumInputError, match="strictly increasing"):
        analyze_spectrum(
            np.array([1, 2, 3, 5, 4], dtype=float),
            np.ones(5),
            np.ones(5),
            templates={},
        )


@pytest.mark.parametrize(
    "wave,flux,sigma",
    [
        ([1, 2, 3, 4, 5], [1, 1, np.nan, 1, 1], [1, 1, 1, 1, 1]),
        ([1, 2, 3, 4, 5], [1, 1, 1, 1, 1], [1, 1, 0, 1, 1]),
        ([1, 2, 3, 4, 5], [1, 1, 1, 1, 1], [1, 1, 1, 1]),
    ],
)
def test_rejects_nonfinite_nonpositive_or_mismatched_inputs(wave, flux, sigma):
    with pytest.raises(SpectrumInputError):
        analyze_spectrum(np.asarray(wave), np.asarray(flux), np.asarray(sigma), templates={})


def test_constant_flux_is_valid_and_fit_is_finite():
    wave = np.linspace(500, 2500, 200)
    flux = np.ones_like(wave)
    sigma = np.full_like(wave, 0.01)
    result = analyze_spectrum(wave, flux, sigma)
    assert result.n_points == len(wave)
    assert result.fits
    for fit in result.fits.values():
        assert np.isfinite(fit["amplitude"])
        assert np.isfinite(fit["amplitude_uncertainty"])
        assert fit["amplitude_uncertainty"] > 0


def test_template_injection_recovers_positive_coefficient():
    wave = np.linspace(500, 2500, 1000)
    template = illustrative_templates(wave)["H2O_toy"]
    sigma = np.full_like(wave, 0.001)
    flux = 1.0 - 0.8 * template
    result = analyze_spectrum(wave, flux, sigma, templates={"injected_toy": template})
    assert result.fits["injected_toy"]["amplitude"] == pytest.approx(-0.8, abs=0.02)


def test_outlier_flagging_flags_single_spike_without_mutating_input():
    values = np.ones(101)
    values[50] = 100.0
    before = values.copy()
    mask = flag_outliers(values, threshold=5)
    assert mask[50]
    np.testing.assert_array_equal(values, before)


def test_csv_loader(tmp_path):
    path = tmp_path / "spectrum.csv"
    path.write_text(
        "wavelength,flux,uncertainty\n1,1.0,0.1\n2,1.1,0.1\n"
        "3,0.9,0.1\n4,1.0,0.1\n5,1.0,0.1\n",
        encoding="utf-8",
    )
    wave, flux, sigma = load_spectrum_csv(path)
    assert len(wave) == len(flux) == len(sigma) == 5


def test_csv_requires_uncertainties(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("wavelength,flux\n1,1\n2,1\n", encoding="utf-8")
    with pytest.raises(SpectrumInputError, match="uncertainty"):
        load_spectrum_csv(path)
