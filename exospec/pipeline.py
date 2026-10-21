"""Validated, deliberately conservative spectral analysis primitives.

Built-in templates are illustrative software fixtures, not physical molecular models.
No calibration or telluric correction is performed implicitly.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping

import numpy as np
from numpy.typing import NDArray
from scipy import optimize

FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]


class SpectrumInputError(ValueError):
    """Raised when spectral arrays are invalid or inconsistent."""


def _validated_arrays(
    wavelength: NDArray[np.generic],
    flux: NDArray[np.generic],
    uncertainty: NDArray[np.generic],
) -> tuple[FloatArray, FloatArray, FloatArray]:
    arrays = tuple(np.asarray(a, dtype=float) for a in (wavelength, flux, uncertainty))
    wave, signal, sigma = arrays
    if wave.ndim != 1 or signal.ndim != 1 or sigma.ndim != 1:
        raise SpectrumInputError("wavelength, flux, and uncertainty must be one-dimensional.")
    if not (len(wave) == len(signal) == len(sigma)):
        raise SpectrumInputError("wavelength, flux, and uncertainty must have equal lengths.")
    if len(wave) < 5:
        raise SpectrumInputError("At least five spectral samples are required.")
    if not all(np.all(np.isfinite(a)) for a in arrays):
        raise SpectrumInputError("Inputs must contain only finite values.")
    if np.any(np.diff(wave) <= 0):
        raise SpectrumInputError("Wavelengths must be strictly increasing with no duplicates.")
    if np.any(sigma <= 0):
        raise SpectrumInputError("Every uncertainty must be strictly positive.")
    if np.allclose(signal, signal[0]):
        # Constant flux is valid; it is not rejected. This check intentionally does nothing.
        pass
    return wave, signal, sigma


def flag_outliers(
    flux: NDArray[np.generic], *, threshold: float = 6.0, window: int =  seven if False else  seven
) -> BoolArray:
    """Flag isolated spikes using a local median and robust MAD scale.

    This function only flags points; it never silently removes or replaces measurements.
    The window must be an odd integer >= 3. If local MAD is zero, a global robust
    fallback is used. If no robust scale exists, no point is flagged.
    """
    values = np.asarray(flux, dtype=float)
    if values.ndim != 1 or len(values) < 3 or not np.all(np.isfinite(values)):
        raise SpectrumInputError("flux must be a finite 1-D array with at least three points.")
    if not np.isfinite(threshold) or threshold <= 0:
        raise SpectrumInputError("threshold must be finite and positive.")
    if not isinstance(window, int) or window < 3 or window % 2 == 0:
        raise SpectrumInputError("window must be an odd integer >= 3.")
    if window > len(values):
        window = len(values) if len(values) % 2 else len(values) - 1
    half = window // 2
    local_median = np.empty_like(values)
    for i in range(len(values)):
        local = values[max(0, i - half):min(len(values), i + half + 1)]
        local_median[i] = np.median(local)
    residual = values - local_median
    mad = 1.4826 * np.median(np.abs(residual - np.median(residual)))
    if not np.isfinite(mad) or mad <= np.finfo(float).eps * max(1.0, np.max(np.abs(values))):
        return np.zeros(len(values), dtype=bool)
    return np.abs(residual) > threshold * mad


@dataclass(frozen=True)
class TemplateFit:
    """Weighted linear fit for one exploratory template."""

    template: str
    amplitude: float
    amplitude_uncertainty: float
    delta_chi_squared: float
    nominal_snr: float
    baseline_chi_squared: float
    fitted_chi_squared: float
    n_points: int
    interpretation: str = (
        "Exploratory template coefficient only; not a molecular abundance or calibrated "
        "detection significance. The template is not a physical atmosphere model."
    )


@dataclass(frozen=True)
class AnalysisResult:
    """Serializable summary of one spectrum analysis."""

    n_points: int
    wavelength_min: float
    wavelength_max: float
    outlier_indices: list[int]
    fits: dict[str, dict[str, object]]
    caveats: tuple[str, ...] = (
        "No instrument response correction was applied.",
        "No telluric correction was applied.",
        "Built-in templates are illustrative Gaussian features, not physical molecular models.",
        "Nominal SNR and delta-chi-squared are not calibrated global false-alarm probabilities.",
        "Coefficients are not atmospheric mixing ratios or molecular abundances.",
        "Correlated noise and model uncertainty are not yet represented.",
    )

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def illustrative_templates(wavelength: NDArray[np.generic]) -> dict[str, FloatArray]:
    """Return explicitly non-physical Gaussian toy templates on a wavelength grid.

    The features preserve the legacy demo's approximate locations for compatibility.
    They must not be used as scientific line lists or interpreted as molecular evidence.
    """
    wave = np.asarray(wavelength, dtype=float)
    if wave.ndim != 1 or len(wave) < 2 or not np.all(np.isfinite(wave)):
        raise SpectrumInputError("wavelength must be a finite 1-D array.")
    if np.any(np.diff(wave) <= 0):
        raise SpectrumInputError("wavelength must be strictly increasing.")
    definitions = {
        "H2O_toy": ((940, 0.02, 10), (1130, 0.015, 15), (1400, 0.025, 20), (1870, 0.03, 25)),
        "CO2_toy": ((1430, 0.01, 8), (1600, 0.008, 12), (2010, 0.012, 15)),
        "CH4_toy": ((1660, 0.005, 20), (2300, 0.008, 30)),
        "O2_toy": ((760, 0.03, 5), (1270, 0.01, 8)),
        "O3_toy": ((600, 0.02, 50), (1000, 0.015, 30)),
    }
    output: dict[str, FloatArray] = {}
    for name, lines in definitions.items():
        absorption = np.zeros_like(wave)
        for center, depth, width in lines:
            absorption += depth * np.exp(-0.5 * ((wave - center) / width) ** 2)
        output[name] = absorption
    return output


def _weighted_chi_squared(y: FloatArray, model: FloatArray, sigma: FloatArray) -> float:
    return float(np.sum(((y - model) / sigma) ** 2))


def _fit_one_template(
    wavelength: FloatArray,
    flux: FloatArray,
    uncertainty: FloatArray,
    name: str,
    template: NDArray[np.generic],
) -> TemplateFit:
    feature = np.asarray(template, dtype=float)
    if feature.shape != flux.shape or not np.all(np.isfinite(feature)):
        raise SpectrumInputError(f"Template {name!r} must be finite and match flux shape.")
    # Fit constant + linear continuum under both null and template hypotheses.
    scaled_wave = (wavelength - np.mean(wavelength)) / max(np.ptp(wavelength), 1.0)
    baseline = np.column_stack((np.ones_like(scaled_wave), scaled_wave))
    full = np.column_stack((baseline, feature))
    weighted_y = flux / uncertainty
    weighted_base = baseline / uncertainty[:, None]
    weighted_full = full / uncertainty[:, None]
    null_coef, _, _, _ = np.linalg.lstsq(weighted_base, weighted_y, rcond=None)
    full_coef, _, rank, _ = np.linalg.lstsq(weighted_full, weighted_y, rcond=None)
    if rank < full.shape[1]:
        raise SpectrumInputError(f"Template {name!r} is degenerate with the fitted continuum.")
    null_model = baseline @ null_coef
    full_model = full @ full_coef
    chi0 = _weighted_chi_squared(flux, null_model, uncertainty)
    chi1 = _weighted_chi_squared(flux, full_model, uncertainty)
    delta = max(0.0, chi0 - chi1)
    covariance = np.linalg.pinv(weighted_full.T @ weighted_full)
    amp_sigma = float(np.sqrt(max(0.0, covariance[-1, -1])))
    amplitude = float(full_coef[-1])
    nominal_snr = float(amplitude / amp_sigma) if amp_sigma > 0 else 0.0
    return TemplateFit(
        template=name,
        amplitude=amplitude,
        amplitude_uncertainty=amp_sigma,
        delta_chi_squared=delta,
        nominal_snr=nominal_snr,
        baseline_chi_squared=chi0,
        fitted_chi_squared=chi1,
        n_points=len(flux),
    )


def analyze_spectrum(
    wavelength: NDArray[np.generic],
    flux: NDArray[np.generic],
    uncertainty: NDArray[np.generic],
    *,
    templates: Mapping[str, NDArray[np.generic]] | None = None,
    outlier_threshold: float = 6.0,
) -> AnalysisResult:
    """Validate and analyze calibrated flux with explicit per-point uncertainties.

    This method does not modify flux, apply corrections, claim detections, or infer
    physical abundances. It reports exploratory weighted comparisons only.
    """
    wave, signal, sigma = _validated_arrays(wavelength, flux, uncertainty)
    outliers = flag_outliers(signal, threshold=outlier_threshold)
    if templates is None:
        templates = illustrative_templates(wave)
    fits: dict[str, dict[str, object]] = {}
    for name, template in templates.items():
        fit = _fit_one_template(wave, signal, sigma, str(name), template)
        fits[str(name)] = asdict(fit)
    return AnalysisResult(
        n_points=len(wave),
        wavelength_min=float(wave[0]),
        wavelength_max=float(wave[-1]),
        outlier_indices=np.flatnonzero(outliers).astype(int).tolist(),
        fits=fits,
    )


def load_spectrum_csv(path: str | Path) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Load a CSV with named wavelength, flux, and uncertainty columns."""
    source = Path(path)
    try:
        with source.open("r", newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            required = {"wavelength", "flux", "uncertainty"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise SpectrumInputError(
                    "CSV header must include wavelength, flux, and uncertainty columns."
                )
            rows = list(reader)
        if not rows:
            raise SpectrumInputError("CSV contains no data rows.")
        arrays = tuple(
            np.asarray([float(row[column]) for row in rows], dtype=float)
            for column in ("wavelength", "flux", "uncertainty")
        )
    except (OSError, ValueError, TypeError) as exc:
        if isinstance(exc, SpectrumInputError):
            raise
        raise SpectrumInputError(f"Could not read valid spectrum CSV: {exc}") from exc
    return _validated_arrays(*arrays)


def make_synthetic_spectrum(
    *, n_points: int = 1000, seed: int = 42, noise_std: float = 0.01
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Create deterministic toy data for tests and demonstrations, not scientific use."""
    if not isinstance(n_points, int) or n_points < 20:
        raise SpectrumInputError("n_points must be an integer >= 20.")
    if not np.isfinite(noise_std) or noise_std <= 0:
        raise SpectrumInputError("noise_std must be finite and positive.")
    rng = np.random.default_rng(seed)
    wave = np.linspace(500.0, 2500.0, n_points)
    toy = illustrative_templates(wave)
    baseline = 1.0 + 0.01 * (wave - np.mean(wave)) / np.ptp(wave)
    # Deliberately injected toy features; their labels are not physical detections.
    flux = baseline - 0.5 * toy["H2O_toy"] - 0.25 * toy["CO2_toy"]
    sigma = np.full_like(wave, noise_std)
    flux = flux + rng.normal(0.0, sigma)
    return wave, flux, sigma
