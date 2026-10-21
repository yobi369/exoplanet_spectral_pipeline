# Exoplanet Spectral Pipeline

A small, auditable Python toolkit for validating one-dimensional spectra and running **exploratory, uncertainty-weighted template comparisons**. It is being developed as a research-software foundation, not as a validated atmospheric-retrieval package.

## Scientific status — read before interpreting results

The current built-in molecular templates are **illustrative Gaussian templates**, not laboratory line lists, radiative-transfer models, or instrument-convolved planetary atmosphere models. Their scores are exploratory matched-filter statistics, **not calibrated detection significances or evidence of molecular presence**. Fitted template amplitudes are not molecular abundances.

The pipeline deliberately does not invent instrument-response or telluric corrections. Supply calibrated data and any justified correction models yourself. For publication-grade work, compare against validated retrieval tools and instrument-specific pipelines, propagate correlated/systematic uncertainties, and validate false-positive rates on injection/recovery experiments.

## Features

- Input validation for finite, ordered wavelength arrays, flux, and positive uncertainties.
- Robust spike/outlier flags using a median/MAD statistic; flagged points are retained unless explicitly replaced.
- Optional uncertainty-weighted template fitting with a linear continuum nuisance term.
- Deterministic synthetic data generation for software tests and demonstrations.
- CSV input and machine-readable JSON output.
- Unit tests and GitHub Actions checks.

## Install

Python 3.10 or newer.

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
python -m pip install -e ".[dev]"
pytest
```

## Run on a CSV spectrum

CSV must contain columns `wavelength`, `flux`, and `uncertainty`. Wavelength units and flux units must be consistent with your data source and documented in your analysis.

```bash
exospec analyze spectrum.csv --output results.json
```

The JSON report includes input summary, flagged outliers, template-fit coefficients, and explicit limitations. It does not apply guessed calibration, denoising, or atmospheric corrections.

## Python API

```python
import numpy as np
from exospec import analyze_spectrum, illustrative_templates

wavelength = np.linspace(500, 2500, 1000)  # replace with calibrated observations
flux = np.ones_like(wavelength)            # replace with measured flux
uncertainty = np.full_like(wavelength, 0.01)

result = analyze_spectrum(
    wavelength, flux, uncertainty,
    templates=illustrative_templates(wavelength),
)
print(result.to_dict())
```

This example demonstrates the interface only; constant flux is not a scientific observation.

## Data sources and provenance

For real projects, use the calibrated products and documentation from the relevant archive. JWST products are available through [MAST](https://archive.stsci.edu/missions-and-data/jwst) and can be retrieved programmatically with [astroquery.mast](https://jwst-docs.stsci.edu/accessing-jwst-data/mast-api-access). Published exoplanet spectra are also indexed by the [NASA Exoplanet Archive atmospheric spectroscopy service](https://exoplanetarchive.ipac.caltech.edu/docs/atmospheres/atmospheres_home.html).

Record at least the archive identifier/DOI, target, instrument/mode, data-product filename, calibration pipeline/version, extraction/reduction reference, units, masks, and every preprocessing choice. Respect the data product's citation and acknowledgement requirements.

## Scope and roadmap

This version establishes safer inputs, transparent statistics, tests, and a reproducible CLI. It does **not** yet provide validated JWST/HST ingestion, physical radiative transfer, Bayesian retrieval, covariance-aware noise models, or an independently benchmarked scientific result. Those require instrument-specific implementation and validation datasets; they must not be simulated by hard-coded correction curves.

## Development

```bash
python -m pip install -e ".[dev]"
pytest
ruff check .
```

## License

MIT. See [LICENSE](LICENSE).
