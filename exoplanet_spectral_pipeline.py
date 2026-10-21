"""Compatibility entry point for the exospec package.

The original prototype used hard-coded simulated calibration/telluric corrections and
uncalibrated correlation scores. Those routines are no longer used because they could
produce scientifically misleading output. Use the validated API or CLI documented in
README.md. Built-in templates are illustrative only.
"""

from exospec import analyze_spectrum, illustrative_templates, make_synthetic_spectrum


def main() -> None:
    wavelength, flux, uncertainty = make_synthetic_spectrum()
    result = analyze_spectrum(
        wavelength,
        flux,
        uncertainty,
        templates=illustrative_templates(wavelength),
    )
    print("Exploratory demonstration only; no scientific detections are claimed.")
    print(f"Samples: {result.n_points}")
    print(f"Wavelength range: {result.wavelength_min:.2f}–{result.wavelength_max:.2f}")
    for caveat in result.caveats:
        print(f"LIMITATION: {caveat}")


if __name__ == "__main__":
    main()
