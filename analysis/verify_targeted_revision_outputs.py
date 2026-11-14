#!/usr/bin/env python3
"""Compare an independent gate/preparation rerun with the publication CSVs."""
import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--rerun-dir", type=Path, required=True)
args = parser.parse_args()
names = ["drift_gate_calibration_replicates.csv", "null_gate_verification_replicates.csv",
         "null_gate_verification_summary.csv", "drift_gate_calibration_summary.csv",
         "batch_gate_calibration_replicates.csv", "preparation_biology_replicates.csv",
         "batch_gate_calibration_summary.csv", "periodic_qc_counterexample.csv",
         "preparation_biology_summary.csv"]
count = 0
for name in names:
    def read(path):
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle))
    expected = read(ROOT / "manuscript/generated/review_validation" / name)
    observed = read(args.rerun_dir / name)
    if expected != observed:
        raise AssertionError(f"Independent rerun differs: {name}")
    count += len(expected)
print(f"PASS: {len(names)} targeted validation CSVs, {count} rows match exactly")
