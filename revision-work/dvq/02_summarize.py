"""Summarize or verify the three-partition paired budget check using saved metrics."""
from pathlib import Path
import csv
import json
import statistics as st
import argparse
import math

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input-dir", type=Path, default=Path(__file__).resolve().parent)
parser.add_argument("--check", action="store_true", help="Verify retained CSV summaries without writing files")
args = parser.parse_args()
root = args.input_dir.resolve()
with (root / "paired_metrics.csv").open() as f:
    rows = list(csv.DictReader(f))
numeric = ["S", "k", "repeat_id", "distortion", "ARI", "NMI", "aggregation_sec", "construction_sec"]
for row in rows:
    for key in numeric:
        row[key] = float(row[key]) if row[key] not in ("NA", "") else math.nan
assert len(rows) == 24
assert all(row["converged"] == "TRUE" for row in rows)

def export(name, items):
    if args.check:
        with (root / name).open() as f:
            saved = list(csv.DictReader(f))
        assert len(saved) == len(items), name
        for expected, actual in zip(items, saved):
            for key, value in expected.items():
                if isinstance(value, (int, float)):
                    assert math.isclose(value, float(actual[key]), rel_tol=1e-12, abs_tol=1e-12), (name, key)
                else:
                    assert value == actual[key], (name, key)
    else:
        with (root / name).open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(items[0]))
            writer.writeheader()
            writer.writerows(items)

summaries = []
deltas = []
for data in ("pbmc", "fashion"):
    for S in (5, 10):
        group = [r for r in rows if r["data"] == data and r["S"] == S]
        for method in ("DVQ", "Summary k-means"):
            g = [r for r in group if r["method"] == method]
            assert len(g) == 3
            summary = {"data": data, "S": S, "k": 10, "method": method, "n_repetitions": len(g)}
            for metric in ("distortion", "ARI", "NMI", "aggregation_sec", "construction_sec"):
                values = [r[metric] for r in g]
                summary[metric + "_mean"] = st.mean(values) if all(map(math.isfinite, values)) else math.nan
                summary[metric + "_sd"] = st.stdev(values) if all(map(math.isfinite, values)) else math.nan
            summaries.append(summary)
        pairs = []
        for rep in (1, 2, 3):
            a = next(r for r in group if r["repeat_id"] == rep and r["method"] == "DVQ")
            b = next(r for r in group if r["repeat_id"] == rep and r["method"] == "Summary k-means")
            assert a["seed"] == b["seed"]
            pairs.append({"data": data, "S": S, "repeat_id": rep,
                          "distortion_difference": a["distortion"] - b["distortion"],
                          "distortion_relative_difference_pct": 100 * (a["distortion"] / b["distortion"] - 1),
                          "ARI_difference": a["ARI"] - b["ARI"],
                          "NMI_difference": a["NMI"] - b["NMI"]})
        deltas.extend(pairs)
export("summary_metrics.csv", summaries)
export("paired_differences.csv", deltas)
delta_summary = []
for data in ("pbmc", "fashion"):
    for S in (5, 10):
        g = [r for r in deltas if r["data"] == data and r["S"] == S]
        r = {"data": data, "S": S, "n_repetitions": len(g)}
        for key in ("distortion_difference", "distortion_relative_difference_pct", "ARI_difference", "NMI_difference"):
            r[key + "_mean"] = st.mean(t[key] for t in g)
            r[key + "_sd"] = st.stdev(t[key] for t in g)
        delta_summary.append(r)
export("paired_difference_summary.csv", delta_summary)

print(json.dumps(delta_summary, indent=2))
print("Saved paired summaries agree." if args.check else "Wrote summary and paired-difference CSV files.")
