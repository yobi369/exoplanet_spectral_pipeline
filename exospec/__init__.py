"""Auditable exploratory spectral validation and template fitting."""

from .pipeline import (
    AnalysisResult,
    SpectrumInputError,
    analyze_spectrum,
    flag_outliers,
    illustrative_templates,
    load_spectrum_csv,
    make_synthetic_spectrum,
)

__all__ = [
    "AnalysisResult",
    "SpectrumInputError",
    "analyze_spectrum",
    "flag_outliers",
    "illustrative_templates",
    "load_spectrum_csv",
    "load_spectrum_fits",
    "make_synthetic_spectrum",
]

def __getattr__(name: str):
    """Load optional FITS support only when explicitly requested."""
    if name == "load_spectrum_fits":
        from .fits_io import load_spectrum_fits

        return load_spectrum_fits
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
