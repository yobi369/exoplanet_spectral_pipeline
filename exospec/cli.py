"""Command-line entry point for reproducible CSV analysis."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .pipeline import analyze_spectrum, load_spectrum_csv


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run exploratory uncertainty-weighted template fits on a calibrated spectrum."
    )
    parser.add_argument("command", choices=["analyze"])
    parser.add_argument("csv", type=Path, help="CSV with wavelength,flux,uncertainty columns")
    parser.add_argument("--output", type=Path, help="write JSON report to this path")
    args = parser.parse_args()
    wavelength, flux, uncertainty = load_spectrum_csv(args.csv)
    report = analyze_spectrum(wavelength, flux, uncertainty).to_dict()
    encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False)
    if args.output:
        args.output.write_text(encoded + "\n", encoding="utf-8")
    else:
        print(encoded)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
