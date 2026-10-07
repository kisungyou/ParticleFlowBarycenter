"""Exact-W2 candidate-restricted digit medoid; run after 01_prepare.R.

Dependencies: numpy, scipy, POT. Run with --rep 1,2,3 separately or all.
Every class and test-class distance row is checkpointed. Inputs and split
indices are serialized from the original R preprocessing and split procedure.
"""
import os
for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[name] = "1"
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time
import warnings
import numpy as np
import scipy
from scipy.spatial.distance import cdist
import ot

ROOT = Path(os.environ.get("MNIST_WORK_DIR", Path(__file__).resolve().parent)).resolve()

def write_json(path, value):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2))
    tmp.replace(path)

def read_inputs():
    cloud_path = ROOT / "point_clouds.csv.gz"
    if not cloud_path.exists():
        raise FileNotFoundError(f"Missing {cloud_path}; run 01_prepare.R first")
    metadata_path = ROOT / "python_metadata.json"
    if metadata_path.exists():
        expected = json.loads(metadata_path.read_text())["input_sha256"]
        assert hashlib.sha256(gzip.decompress(cloud_path.read_bytes())).hexdigest() == expected["point_clouds.csv"], "Point-cloud input checksum differs from saved results"
        assert hashlib.sha256((ROOT / "splits.json").read_bytes()).hexdigest() == expected["splits.json"], "Split checksum differs from saved results"
    pts = np.loadtxt(cloud_path, delimiter=",", skiprows=1)
    cuts = np.flatnonzero(np.r_[True, np.diff(pts[:, 0]) != 0, True])
    clouds = {int(pts[a, 0]): np.ascontiguousarray(pts[a:b, 1:])
              for a, b in zip(cuts[:-1], cuts[1:])}
    images = np.loadtxt(ROOT / "images.csv", delimiter=",", skiprows=1, dtype=int)
    labels = dict(zip(images[:, 0].tolist(), images[:, 1].tolist()))
    return clouds, labels, json.loads((ROOT / "splits.json").read_text())

def w2sq(x, y):
    return float(ot.emd2(np.full(len(x), 1.0 / len(x)),
                        np.full(len(y), 1.0 / len(y)),
                        cdist(x, y, metric="sqeuclidean"),
                        numItermax=1000000, numThreads=1,
                        check_marginals=True))

def validate(clouds, save=False):
    refs = np.loadtxt(ROOT / "r_reference_distances.csv", delimiter=",", skiprows=1)
    t0 = time.perf_counter()
    computed = np.array([w2sq(clouds[int(i)], clouds[int(j)]) for i, j, _ in refs])
    err = np.abs(computed - refs[:, 2])
    assert np.max(err) < 1e-9, f"R/Python discrepancy: {np.max(err)}"
    result = {"pairs": len(refs), "max_abs_squared_W2_difference": float(err.max()),
              "elapsed_sec": time.perf_counter() - t0,
              "threshold": 1e-9, "passed": True}
    if save:
        write_json(ROOT / "distance_validation.json", result)
    return result

def score(pred, actual):
    cm = np.zeros((10, 10), dtype=int)
    np.add.at(cm, (actual, pred), 1)
    tp = np.diag(cm)
    precision = np.divide(tp, cm.sum(axis=0), where=cm.sum(axis=0) != 0,
                          out=np.zeros(10, dtype=float))
    recall = np.divide(tp, cm.sum(axis=1), where=cm.sum(axis=1) != 0,
                       out=np.zeros(10, dtype=float))
    f1 = np.divide(2 * precision * recall, precision + recall,
                   where=(precision + recall) != 0, out=np.zeros(10, dtype=float))
    return {"Accuracy": float(tp.sum() / cm.sum()),
            "MacroPrecision": float(precision.mean()),
            "MacroRecall": float(recall.mean()), "MacroF1": float(f1.mean())}, cm

def run_rep(split, clouds, labels):
    rep = split["rep"]
    out = ROOT / f"rep{rep}"
    out.mkdir(exist_ok=True)
    train_ids = np.array(split["train_idx"], dtype=int)
    test_ids = np.array(split["test_idx"], dtype=int)
    assert not set(train_ids) & set(test_ids)
    prototypes = []
    class_records = []
    for lab in range(10):
        fp = out / f"prototype_{lab}.json"
        if fp.exists():
            record = json.loads(fp.read_text())
        else:
            ids = train_ids[np.array([labels[int(i)] == lab for i in train_ids])]
            candidates = split["candidates"][str(lab)]
            assert len(ids) == 100 and len(candidates) == 10
            t0 = time.perf_counter()
            costs = np.array([[0.0 if int(i) == int(j) else
                               w2sq(clouds[int(i)], clouds[int(j)])
                               for j in ids] for i in candidates])
            obj = costs.mean(axis=1)
            best = int(np.argmin(obj))
            selected = int(candidates[best])
            elapsed = time.perf_counter() - t0
            record = {"rep": rep, "label": lab, "selected_image": selected,
                      "selected_atoms": len(clouds[selected]),
                      "candidates": candidates, "training_ids": ids.tolist(),
                      "objectives": obj.tolist(), "selected_objective": float(obj[best]),
                      "construction_sec": elapsed, "exact_ot_solves": 990}
            np.save(out / f"training_costs_{lab}.npy", costs)
            write_json(fp, record)
            print(f"rep {rep}, digit {lab}: selected image {selected}, "
                  f"{len(clouds[selected])} atoms; construction {elapsed:.2f}s", flush=True)
        prototypes.append(record["selected_image"])
        class_records.append(record)
    rows = []
    classification_sec = 0.0
    for lab, proto in enumerate(prototypes):
        fp = out / f"test_distances_{lab}.npz"
        if fp.exists():
            saved = np.load(fp)
            distances = saved["distances"]
            elapsed = float(saved["elapsed_sec"])
        else:
            t0 = time.perf_counter()
            distances = np.array([w2sq(clouds[proto], clouds[int(j)]) for j in test_ids])
            elapsed = time.perf_counter() - t0
            np.savez(fp, distances=distances, elapsed_sec=elapsed)
            print(f"rep {rep}, test distances to digit {lab}: {elapsed:.2f}s", flush=True)
        rows.append(distances)
        classification_sec += elapsed
    dmat = np.array(rows)
    pred = dmat.argmin(axis=0)
    actual = np.array([labels[int(i)] for i in test_ids])
    scores, cm = score(pred, actual)
    metrics = {"rep": rep, "method": "10-candidate_point-cloud_medoid",
               "n_train": len(train_ids), "n_test": len(test_ids),
               "construction_sec": sum(r["construction_sec"] for r in class_records),
               "classification_sec": classification_sec,
               "mean_atoms": float(np.mean([r["selected_atoms"] for r in class_records])),
               "min_atoms": min(r["selected_atoms"] for r in class_records),
               "max_atoms": max(r["selected_atoms"] for r in class_records), **scores}
    write_json(out / "metrics.json", metrics)
    np.savetxt(out / "predictions.csv", np.column_stack([test_ids, actual, pred]),
               delimiter=",", header="image_id,actual,predicted", comments="", fmt="%d")
    np.savetxt(out / "confusion.csv", cm, delimiter=",", fmt="%d")
    print(json.dumps(metrics), flush=True)
    return metrics

def summarize():
    paths = sorted(ROOT.glob("rep[123]/metrics.json"))
    metrics = [json.loads(p.read_text()) for p in paths]
    if not metrics:
        return
    with (ROOT / "split_metrics.csv").open("w") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(metrics[0]))
        writer.writeheader(); writer.writerows(metrics)
    numeric = [k for k in metrics[0] if k not in ("rep", "method")]
    summary = {"completed_splits": len(metrics), "method": metrics[0]["method"]}
    for k in numeric:
        values = np.array([r[k] for r in metrics])
        summary[k + "_mean"] = float(values.mean())
        summary[k + "_sd"] = float(values.std(ddof=1)) if len(values) > 1 else None
    write_json(ROOT / "summary.json", summary)
    return summary

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rep", type=int, choices=[1, 2, 3])
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    warnings.simplefilter("error", UserWarning)
    clouds, labels, splits = read_inputs()
    validation = validate(clouds, save=not args.validate_only)
    print("Exact distance validation:", validation, flush=True)
    cpu = platform.processor() or platform.machine()
    if platform.system() == "Darwin":
        try:
            cpu_read = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], text=True, capture_output=True)
            if cpu_read.returncode == 0:
                cpu = cpu_read.stdout.strip()
        except OSError:
            pass
    metadata = {"python": platform.python_version(), "numpy": np.__version__,
                "scipy": scipy.__version__, "POT": ot.__version__,
                "platform": platform.platform(), "CPU": cpu,
                "OT": "POT emd2, network simplex, float64, numThreads=1, numItermax=1000000",
                "input_sha256": {
                    "point_clouds.csv": hashlib.sha256(gzip.decompress((ROOT / "point_clouds.csv.gz").read_bytes())).hexdigest(),
                    "splits.json": hashlib.sha256((ROOT / "splits.json").read_bytes()).hexdigest()},
                "timing": "Single-process sequential class construction and classification; preprocessing, file I/O, and distance validation excluded; stored checkpoints preserve measured computation timings.",
                "thread_environment": {k: os.environ[k] for k in
                    ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                     "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"]}}
    if not args.validate_only:
        # Keep the original machine/timing provenance when reusing saved checkpoints.
        write_json(ROOT / "execution_metadata.json", metadata)
        if not (ROOT / "python_metadata.json").exists():
            write_json(ROOT / "python_metadata.json", metadata)
        for split in splits:
            if args.rep is None or split["rep"] == args.rep:
                run_rep(split, clouds, labels)
                summarize()
