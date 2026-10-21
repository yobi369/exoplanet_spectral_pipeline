"""Optional FITS-table ingestion for calibrated one-dimensional spectra.

This module reads arrays from an explicitly selected FITS table. It does not
calibrate, extract, normalize, mask, or reinterpret the source product.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from .pipeline import SpectrumInputError, _validated_arrays

FloatArray = NDArray[np.float64]


def load_spectrum_fits(
    path: str | Path,
    *,
    hdu: int = 1,
    wavelength_column: str = "WAVELENGTH",
    flux_column: str = "FLUX",
    uncertainty_column: str = "FLUX_ERROR",
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Load wavelength, flux, and 1-sigma uncertainty columns from a FITS table.

    Requires the optional `astropy` dependency. Column names are matched
    case-insensitively. Ambiguous/missing columns, non-1-D columns, masked or
    invalid values, and invalid wavelength/uncertainty arrays are rejected.
    """
    try:
        from astropy.io import fits
    except ImportError as exc:
        raise ImportError(
            "FITS input requires Astropy. Install it with "
            "`python -m pip install 'exoplanet-spectral-pipeline[astronomy]'`."
        ) from exc

    if not isinstance(hdu, int) or hdu < 0:
        raise SpectrumInputError("hdu must be a non-negative integer.")
    requested = {
        "wavelength": wavelength_column,
        "flux": flux_column,
        "uncertainty": uncertainty_column,
    }
    if any(not isinstance(name, str) or not name.strip() for name in requested.values()):
        raise SpectrumInputError("FITS column names must be non-empty strings.")

    try:
        with fits.open(Path(path), memmap=False) as hdul:
            if hdu >= len(hdul):
                raise SpectrumInputError(f"FITS HDU {hdu} does not exist.")
            data = getattr(hdul[hdu], "data", None)
            names = getattr(getattr(data, "dtype", None), "names", None)
            if data is None or not names:
                raise SpectrumInputError(f"FITS HDU {hdu} is not a binary/table HDU.")
            lookup: dict[str, str] = {}
            for name in names:
                key = str(name).casefold()
                if key in lookup:
                    raise SpectrumInputError(
                        f"FITS HDU has ambiguous case-insensitive column name {name!r}."
                    )
                lookup[key] = str(name)
            selected: dict[str, FloatArray] = {}
            for role, column in requested.items():
                actual = lookup.get(column.casefold())
                if actual is None:
                    available = ", ".join(str(name) for name in names)
                    raise SpectrumInputError(
                        f"Missing FITS {role} column {column!r}; available columns: {available}."
                    )
                values = np.asarray(data[actual])
                if values.ndim != 1:
                    raise SpectrumInputError(
                        f"FITS column {actual!r} must be one-dimensional; got shape {values.shape}."
                    )
                selected[role] = np.asarray(values, dtype=float)
    except SpectrumInputError:
        raise
    except (OSError, ValueError, TypeError) as exc:
        raise SpectrumInputError(f"Could not read FITS spectrum {str(path)!r}: {exc}") from exc

    return _validated_arrays(
        selected["wavelength"], selected["flux"], selected["uncertainty"]
    )
