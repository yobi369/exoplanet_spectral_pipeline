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

## FITS binary-table input

Install the optional astronomy dependency:

```bash
python -m pip install -e ".[astronomy]"
```

For a calibrated FITS product whose selected binary-table HDU contains one row per wavelength sample:

```bash
exospec analyze-fits jwst_x1d.fits --hdu 1 \
  --wavelength-column WAVELENGTH --flux-column FLUX \
  --uncertainty-column FLUX_ERROR --output results.json
```

Column names are matched case-insensitively. The HDU and columns are configurable because FITS products differ by instrument, pipeline, and product type. Only one-dimensional table columns are accepted; image HDUs, multi-order/multi-aperture arrays, missing columns, non-finite values, non-increasing wavelengths, and non-positive uncertainties are rejected rather than guessed or silently reshaped. Some products store uncertainty under another column name, such as `FLUX_ERR`; pass the actual name with `--uncertainty-column`.

This is a reader, not an instrument pipeline: it does not perform extraction, wavelength calibration, flux calibration, background subtraction, quality masking, telluric correction, or unit conversion. Confirm that the selected product is already calibrated and that its uncertainty column is a one-sigma uncertainty in the flux units. The JSON report records the input format, selected HDU, and column mapping but cannot verify scientific provenance.

## Example

```bash
exospec analyze spectrum.csv --output results.json
```

The tool does not automatically normalize, denoise, calibrate, mask, or correct telluric absorption. Apply only scientifically justified preprocessing for the relevant instrument, and document each step and its parameters.
