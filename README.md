# Exoplanet Spectral Pipeline

An auditable Python toolkit for validating one-dimensional spectra and running **exploratory, uncertainty-weighted template comparisons** from CSV or selected FITS binary-table products. It is a research-software foundation, not a validated atmospheric-retrieval package.

## Scientific status — read before interpreting results

The built-in molecular templates are **illustrative Gaussian templates**, not laboratory line lists, radiative-transfer models, or instrument-convolved planetary atmosphere models. Their scores are exploratory matched-filter statistics, **not calibrated detection significances or evidence of molecular presence**. Fitted template amplitudes are not molecular abundances.

The pipeline deliberately does not invent instrument-response or telluric corrections. Supply calibrated data and justified correction models yourself. For publication-grade work, compare against validated retrieval tools and instrument-specific pipelines, propagate correlated/systematic uncertainties, and validate false-positive rates on injection/recovery experiments.

## Features

- Input validation for finite, ordered wavelength arrays, flux, and positive uncertainties.
- Robust spike/outlier flags; flagged points are retained rather than silently deleted.
- Optional uncertainty-weighted template fitting with a linear continuum nuisance term.
- CSV ingestion plus optional FITS binary-table ingestion with explicit HDU/column selection.
- JSON reports that preserve input format and FITS column mapping.
- Deterministic synthetic data, a one-command demo, tests, and GitHub Actions checks.

## Install

Python 3.10 or newer.

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
python -m pip install -e ".[dev]"
pytest
```

For FITS support, install the optional astronomy dependency:

```bash
python -m pip install -e ".[astronomy]"
```

## Run a complete local demo

```bash
exospec demo --output-dir exospec-demo --seed 42
```

This creates a seeded synthetic spectrum and JSON report. It is **synthetic software-test data**, not an observation or scientific evidence.

## Analyze a CSV spectrum

CSV must contain `wavelength`, `flux`, and `uncertainty` columns:

```bash
exospec analyze spectrum.csv --output results.json
```

## Analyze a FITS binary table

For a calibrated FITS product with one-dimensional table columns:

```bash
exospec analyze-fits spectrum.fits --hdu 1 \
  --wavelength-column WAVELENGTH --flux-column FLUX \
  --uncertainty-column FLUX_ERROR --output results.json
```

HDU and column names vary by instrument and product; inspect the source product's documentation and pass the correct names. The reader rejects unsupported layouts instead of silently reshaping them. It does not calibrate or extract data, apply masks, correct telluric absorption, convert units, or verify provenance. FITS support requires `astropy`.

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

The constant-flux example demonstrates the interface only; it is not a scientific observation.

## Data provenance

JWST calibrated products are available through [MAST](https://archive.stsci.edu/missions-and-data/jwst); see the [MAST API documentation](https://jwst-docs.stsci.edu/accessing-jwst-data/mast-api-access). Published spectra are also indexed by the [NASA Exoplanet Archive atmospheric spectroscopy service](https://exoplanetarchive.ipac.caltech.edu/docs/atmospheres/atmospheres_home.html).

Record the archive identifier/DOI, target, instrument/mode, product filename, calibration pipeline/version, reduction reference, units, masks, and preprocessing choices. Follow the product's citation requirements.

## Scope and roadmap

This release adds a conservative FITS table reader for calibrated one-dimensional inputs. It does **not** yet provide instrument-specific JWST/HST ingestion, physical radiative transfer, Bayesian retrieval, covariance-aware noise models, or an independently benchmarked scientific result. Those require instrument-specific implementation and validation datasets; they must not be simulated by hard-coded correction curves.

## Development

```bash
python -m pip install -e ".[dev]"
pytest
ruff check .
```

## License

MIT. See [LICENSE](LICENSE).
