# Scientific validation plan

## Current release boundary

The package provides input validation, outlier flagging, and uncertainty-weighted fits to explicitly illustrative templates. It is **not a validated exoplanet atmospheric retrieval system**. It has not yet been benchmarked against an established retrieval tool or used to publish a scientific result.

## Required gates before research-grade claims

### Gate 1 — software correctness
- Run unit tests and CI on supported Python versions.
- Test malformed inputs, non-finite values, non-positive uncertainties, duplicate/non-monotonic wavelengths, edge cases, and serialization.
- Keep all synthetic random generators seeded in tests.
- Keep input data immutable; all flags and preprocessing steps must be reported.

### Gate 2 — observational-data provenance
- Select a published, public spectrum with an archive identifier and paper citation.
- Preserve original data unchanged and record source URL/DOI, retrieval date, file checksum, target, instrument, observing mode, units, reduction pipeline/version, and any masks.
- Parse wavelength, measurement, and uncertainty columns explicitly; reject unknown or ambiguous units rather than guessing.
- Do not apply synthetic instrument or telluric corrections to archival data.

### Gate 3 — benchmark and calibration
- Compare output to the published spectrum and a documented reference implementation appropriate to the same data product.
- Run injection/recovery experiments across signal amplitude, wavelength coverage, resolution, uncertainty, and correlated-noise regimes.
- Measure false-positive rates on null/no-injection data and report sensitivity/false-positive trade-offs.
- Calibrate any detection statistic for the actual search procedure, including the number of templates and wavelength ranges tested. A nominal coefficient/error ratio is not a global detection significance.

### Gate 4 — physical inference
- Replace toy templates with versioned, sourced molecular opacity/line-list data and a documented radiative-transfer/forward model.
- Convolve the forward model with the instrument line-spread function and match the data's spectral resolution and binning.
- Model stellar contamination, clouds/hazes, instrument systematics, and correlated noise where relevant to the data.
- Use uncertainty-aware retrieval with parameter priors, convergence checks, posterior diagnostics, and posterior predictive checks.
- Validate retrieved parameters on known benchmark cases and report sensitivity to model assumptions.

### Gate 5 — independent review
- Publish the exact command, environment, configuration, dataset identifier, and expected output needed to reproduce each benchmark.
- Ask an astronomer or atmospheric-retrieval specialist to review the method before making scientific claims.
- Keep software quality claims separate from claims about atmospheric composition.

## Data provenance record template

For every real dataset, record:

- Target and observing program:
- Archive and persistent identifier/DOI:
- Data product filename(s) and SHA-256:
- Instrument, mode, wavelength units, flux/transit-depth units:
- Pipeline and calibration-reference versions:
- Original paper and reduction reference:
- Input columns and uncertainty interpretation:
- Masks, bins, and any transformations:
- Model/template source and version:
- Software commit, Python version, dependency versions:
- Known caveats and comparison result:

## Decision rule

Until all relevant gates are passed, describe outputs as exploratory fits or software-validation results—not molecule detections, atmospheric abundances, habitability assessments, or discoveries.
