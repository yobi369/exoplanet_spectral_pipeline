import sys
import types

import numpy as np
import pytest

from exospec import SpectrumInputError
from exospec.fits_io import load_spectrum_fits


class FakeHDU:
    def __init__(self, data):
        self.data = data


class FakeHDUList(list):
    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False


def install_fake_astropy(monkeypatch, hdus):
    fits = types.ModuleType("astropy.io.fits")
    fits.open = lambda *_args, **_kwargs: FakeHDUList(hdus)
    astropy = types.ModuleType("astropy")
    astropy_io = types.ModuleType("astropy.io")
    astropy.io = astropy_io
    astropy_io.fits = fits
    monkeypatch.setitem(sys.modules, "astropy", astropy)
    monkeypatch.setitem(sys.modules, "astropy.io", astropy_io)
    monkeypatch.setitem(sys.modules, "astropy.io.fits", fits)


def sample_table():
    return np.array(
        [
            (1.0, 1.01, 0.1),
            (2.0, 0.99, 0.1),
            (3.0, 1.00, 0.1),
            (4.0, 1.02, 0.1),
            (5.0, 0.98, 0.1),
        ],
        dtype=[
            ("wavelength", "f8"),
            ("flux", "f8"),
            ("flux_error", "f8"),
        ],
    )


def test_fits_loader_selects_columns_case_insensitively(monkeypatch, tmp_path):
    install_fake_astropy(monkeypatch, [FakeHDU(None), FakeHDU(sample_table())])
    wave, flux, uncertainty = load_spectrum_fits(
        tmp_path / "spectrum.fits",
        wavelength_column="WAVELENGTH",
        flux_column="FLUX",
        uncertainty_column="FLUX_ERROR",
    )
    np.testing.assert_array_equal(wave, [1, 2, 3, 4, 5])
    np.testing.assert_allclose(flux, [1.01, 0.99, 1.00, 1.02, 0.98])
    np.testing.assert_array_equal(uncertainty, np.full(5, 0.1))


def test_fits_loader_rejects_missing_columns(monkeypatch, tmp_path):
    install_fake_astropy(monkeypatch, [FakeHDU(None), FakeHDU(sample_table())])
    with pytest.raises(SpectrumInputError, match="Missing FITS uncertainty column"):
        load_spectrum_fits(
            tmp_path / "spectrum.fits", uncertainty_column="NO_SUCH_COLUMN"
        )


def test_fits_loader_rejects_nonpositive_uncertainty(monkeypatch, tmp_path):
    table = sample_table()
    table["flux_error"][2] = 0.0
    install_fake_astropy(monkeypatch, [FakeHDU(None), FakeHDU(table)])
    with pytest.raises(SpectrumInputError, match="strictly positive"):
        load_spectrum_fits(tmp_path / "spectrum.fits")


def test_fits_loader_rejects_multidimensional_columns(monkeypatch, tmp_path):
    table = np.zeros(
        5,
        dtype=[("wavelength", "f8"), ("flux", "f8", (2,)), ("flux_error", "f8")],
    )
    install_fake_astropy(monkeypatch, [FakeHDU(None), FakeHDU(table)])
    with pytest.raises(SpectrumInputError, match="one-dimensional"):
        load_spectrum_fits(tmp_path / "spectrum.fits")


def test_fits_loader_rejects_missing_hdu(monkeypatch, tmp_path):
    install_fake_astropy(monkeypatch, [FakeHDU(None)])
    with pytest.raises(SpectrumInputError, match="does not exist"):
        load_spectrum_fits(tmp_path / "spectrum.fits", hdu=2)
