"""Command-line interface for reproducible exploratory spectrum analysis."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from .pipeline import analyze_spectrum, load_spectrum_csv, make_synthetic_spectrum


def _write_report(report: dict[str, object], output: Path | None) -> None:
    encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded, encoding="utf-8")
    else:
        print(encoded, end="")


def _demo(output_dir: Path, seed: int) -> None:
    """Write a fully reproducible synthetic demonstration and its analysis report."""
    output_dir.mkdir(parents=True, exist_ok=True)
    wavelength, flux, uncertainty = make_synthetic_spectrum(seed=seed)
    csv_path = output_dir / "synthetic_demo.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["wavelength", "flux", "uncertainty"])
        writer.writerows(zip(wavelength, flux, uncertainty, strict=True))
    report = analyze_spectrum(wavelength, flux, uncertainty).to_dict()
    report["input_file"] = csv_path.name
    report["data_status"] = (
        "SYNTHETIC DEMONSTRATION ONLY — not telescope observations and not scientific evidence."
    )
    _write_report(report, output_dir / "results.json")
    print(f"Wrote synthetic demo spectrum: {csv_path}")
    print(f"Wrote exploratory report: {output_dir / 'results.json'}")
    print("WARNING: synthetic data only; template fits are not calibrated molecular detections.")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Validate calibrated one-dimensional spectra and run exploratory "
            "uncertainty-weighted template fits."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    analyze = subparsers.add_parser("analyze", help="analyze a CSV spectrum")
    analyze.add_argument(
        "csv", type=Path, help="CSV with wavelength,flux,uncertainty columns"
    )
    analyze.add_argument("--output", type=Path, help="write JSON report to this path")

    demo = subparsers.add_parser(
        "demo", help="generate reproducible synthetic data and analyze it"
    )
    demo.add_argument("--output-dir", type=Path, default=Path("exospec-demo"))
    demo.add_argument("--seed", type=int, default=42)

    args = parser.parse_args()
    if args.command == "demo":
        _demo(args.output_dir, args.seed)
        return 0

    wavelength, flux, uncertainty = load_spectrum_csv(args.csv)
    report = analyze_spectrum(wavelength, flux, uncertainty).to_dict()
    report["input_file"] = args.csv.name
    report["data_status"] = (
        "USER-SUPPLIED DATA; provenance and calibration were not independently verified."
    )
    _write_report(report, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
