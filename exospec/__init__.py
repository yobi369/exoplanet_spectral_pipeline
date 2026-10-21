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
    "make_synthetic_spectrum",
]
