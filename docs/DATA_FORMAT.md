# Input data format

The command line accepts CSV files with exactly named required columns:

```csv
wavelength,flux,uncertainty
1000.0,1.002,0.010
1001.0,0.998,0.011
1002.0,1.001,0.009
1003.0,0.999,0.010
1004.0,1.003,0.012
```

The five rows above are **format-only example values**, not observed data. Use at least five valid samples; real analyses should generally use the full appropriate dataset.

## Required semantics

- `wavelength`: finite, strictly increasing values, with a single documented unit.
- `flux`: finite measurement values, in the units defined by the source product.
- `uncertainty`: finite, strictly positive, one-sigma uncertainty per sample in the same units as flux.

This version assumes independent per-point uncertainties during fitting. It does not model covariance between spectral bins. If the source data has covariance information, preserve it separately and do not assume this tool accounts for it.

## Example

```bash
exospec analyze spectrum.csv --output results.json
```

The tool does not automatically normalize, denoise, calibrate, mask, or correct telluric absorption. Apply only scientifically justified preprocessing for the relevant instrument, and document each step and its parameters.
