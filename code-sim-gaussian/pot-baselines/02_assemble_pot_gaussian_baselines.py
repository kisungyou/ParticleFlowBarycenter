#!/usr/bin/env python3
"""Assemble P1 POT baseline metrics into one CSV."""
from __future__ import annotations
import csv
from pathlib import Path
import os

def run_suffix() -> str:
    label = os.environ.get("REV_RUN_LABEL", "")
    return f"-{label}" if label else ""

script_dir = Path(__file__).resolve().parent
run_dir = script_dir / f"runs{run_suffix()}"
files = sorted(run_dir.glob("result_*.csv"))
if not files:
    raise SystemExit(f"No result CSV files found in {run_dir}")
rows = []
for fp in files:
    with fp.open(newline="") as f:
        rows.extend(csv.DictReader(f))
out = script_dir / f"assembled_pot_gaussian_baselines{run_suffix()}.csv"
with out.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)
print(f"Wrote {out} with {len(rows)} rows")
