"""Rebuild or check medoid summaries from saved predictions; no OT solves.

Uses only the Python standard library. --check leaves saved evidence unchanged.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import statistics

ROOT = Path(os.environ.get("MNIST_WORK_DIR", Path(__file__).resolve().parent)).resolve()


def assemble():
    records = []
    for rep in (1, 2, 3):
        folder = ROOT / f"rep{rep}"
        record = json.loads((folder / "metrics.json").read_text())
        with (folder / "predictions.csv").open() as fp:
            predictions = list(csv.DictReader(fp))
        assert len(predictions) == 2500
        confusion = [[0] * 10 for _ in range(10)]
        for row in predictions:
            confusion[int(row["actual"])][int(row["predicted"])] += 1
        precision, recall, f1 = [], [], []
        for lab in range(10):
            tp = confusion[lab][lab]
            predicted = sum(row[lab] for row in confusion)
            actual = sum(confusion[lab])
            p, r = (tp / predicted if predicted else 0), tp / actual
            precision.append(p)
            recall.append(r)
            f1.append(2 * p * r / (p + r) if p + r else 0)
        recomputed = {
            "Accuracy": sum(confusion[i][i] for i in range(10)) / len(predictions),
            "MacroPrecision": statistics.mean(precision),
            "MacroRecall": statistics.mean(recall),
            "MacroF1": statistics.mean(f1),
        }
        for key, value in recomputed.items():
            assert math.isclose(value, record[key], rel_tol=0, abs_tol=1e-12), (rep, key)
        records.append(record)
    summary = {"completed_splits": 3, "method": records[0]["method"]}
    for key in records[0]:
        if key not in ("rep", "method"):
            values = [row[key] for row in records]
            summary[key + "_mean"] = statistics.mean(values)
            summary[key + "_sd"] = statistics.stdev(values)
    return records, summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    records, summary = assemble()
    if args.check:
        saved = json.loads((ROOT / "summary.json").read_text())
        for key, value in summary.items():
            if isinstance(value, (int, float)):
                assert math.isclose(value, saved[key], rel_tol=1e-12, abs_tol=1e-12), key
            else:
                assert value == saved[key], key
        with (ROOT / "split_metrics.csv").open() as fp:
            rows = list(csv.DictReader(fp))
        assert len(rows) == 3
        for record, row in zip(records, rows):
            for key, value in record.items():
                if isinstance(value, (int, float)):
                    assert math.isclose(value, float(row[key]), rel_tol=1e-12, abs_tol=1e-12), key
                else:
                    assert value == row[key], key
        print("Saved predictions, split metrics, means, and sample SDs agree.")
    else:
        with (ROOT / "split_metrics.csv").open("w", newline="") as fp:
            writer = csv.DictWriter(fp, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
        (ROOT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
