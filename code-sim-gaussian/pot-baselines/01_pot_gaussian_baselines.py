#!/usr/bin/env python3
"""Revision Experiment P1: POT free-support and Sinkhorn baselines.

Purpose
-------
Benchmark representative existing OT barycenter solvers on the same Gaussian
oracle setting used in the paper. This directly addresses reviewer requests for
POT/entropic/free-support baselines. Results are saved per task; visualization is
left for later.

Usage
-----
python 01_pot_gaussian_baselines.py <task_id>
"""

from __future__ import annotations

import csv
import math
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

try:
    import ot
except ImportError as exc:  # pragma: no cover
    raise SystemExit("Missing Python package 'POT' (import name: ot). Install with: pip install POT") from exc

try:
    from scipy.linalg import sqrtm
except ImportError as exc:  # pragma: no cover
    raise SystemExit("Missing scipy. Install with: pip install scipy") from exc

try:
    from sklearn.cluster import KMeans
except Exception:  # pragma: no cover
    KMeans = None


@dataclass(frozen=True)
class Config:
    rep: int
    support: int
    method: str
    reg: float


def is_smoke() -> bool:
    return os.environ.get("REV_SMOKE", "0").lower() in {"1", "true", "yes", "y", "smoke"}

def run_suffix() -> str:
    label = os.environ.get("REV_RUN_LABEL", "")
    return f"-{label}" if label else ""


def experiment_scale() -> str:
    return os.environ.get("REV_EXPERIMENT_SCALE", "core").lower()


def make_grid() -> list[Config]:
    scale = experiment_scale()
    if scale == "heavy":
        reps = range(1, 31)
        supports = [50, 100, 200, 500]
        sinkhorn_regs = [0.1, 1.0, 5.0, 20.0]
    else:
        # Core reviewer-response baseline grid: representative exact POT and
        # entropic/Sinkhorn comparisons without making this a full solver study.
        reps = range(1, 4)
        supports = [50, 100, 200]
        sinkhorn_regs = [1.0, 5.0]

    configs: list[Config] = []
    for rep in reps:
        for support in supports:
            configs.append(Config(rep=rep, support=support, method="pot_exact", reg=math.nan))
            for reg in sinkhorn_regs:
                configs.append(Config(rep=rep, support=support, method="pot_sinkhorn", reg=reg))
    return configs


def spd_sqrt(A: np.ndarray) -> np.ndarray:
    out = sqrtm((A + A.T) / 2.0).real
    return (out + out.T) / 2.0


def spd_invsqrt(A: np.ndarray) -> np.ndarray:
    vals, vecs = np.linalg.eigh((A + A.T) / 2.0)
    vals = np.maximum(vals, 1e-14)
    return (vecs * (1.0 / np.sqrt(vals))) @ vecs.T


def gaussian_barycenter_cov(covs: np.ndarray, weights: np.ndarray | None = None,
                            maxiter: int = 500, tol: float = 1e-12) -> np.ndarray:
    # covs shape: (K, d, d)
    K, d, _ = covs.shape
    if weights is None:
        weights = np.ones(K) / K
    S = np.mean(covs, axis=0)
    S = (S + S.T) / 2.0
    for _ in range(maxiter):
        S_sqrt = spd_sqrt(S)
        A = np.zeros((d, d))
        for k in range(K):
            A += weights[k] * spd_sqrt(S_sqrt @ covs[k] @ S_sqrt)
        S_next = spd_invsqrt(S) @ (A @ A) @ spd_invsqrt(S)
        S_next = (S_next + S_next.T) / 2.0
        if np.linalg.norm(S_next - S, ord="fro") <= tol * (1.0 + np.linalg.norm(S, ord="fro")):
            S = S_next
            break
        S = S_next
    return S


def bures(A: np.ndarray, B: np.ndarray) -> float:
    A = (A + A.T) / 2.0
    B = (B + B.T) / 2.0
    As = spd_sqrt(A)
    middle = spd_sqrt(As @ B @ As)
    val = np.trace(A) + np.trace(B) - 2.0 * np.trace(middle)
    return float(np.sqrt(max(0.0, val)))


def pairwise_sqeuclidean(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    X2 = np.sum(X * X, axis=1)[:, None]
    Y2 = np.sum(Y * Y, axis=1)[None, :]
    M = X2 + Y2 - 2.0 * X @ Y.T
    return np.maximum(M, 0.0)


def make_gaussian_problem(rep: int, n_per_measure: int = 500, d: int = 2,
                          num_measures: int = 4, loc_scale: float = 10.0):
    rng = np.random.default_rng(700000 + rep)
    signs = np.array([[1, 1], [1, -1], [-1, 1], [-1, -1]], dtype=float)
    means = np.zeros((num_measures, d))
    covs = np.zeros((num_measures, d, d))
    measures: list[np.ndarray] = []
    for k in range(num_measures):
        sign = np.ones(d)
        sign[:2] = signs[k % 4, :2]
        means[k] = loc_scale * sign + rng.normal(scale=1.0, size=d)
        A = rng.normal(size=(max(d + 2, 4), d))
        cov = A.T @ A
        covs[k] = cov
        measures.append(rng.multivariate_normal(means[k], cov, size=n_per_measure))
    weights = np.ones(num_measures) / num_measures
    bary_mean = weights @ means
    bary_cov = gaussian_barycenter_cov(covs, weights)
    return measures, means, covs, bary_mean, bary_cov


def initial_support(measures: list[np.ndarray], support: int, seed: int) -> np.ndarray:
    pooled = np.vstack(measures)
    rng = np.random.default_rng(seed)
    if KMeans is not None and support <= pooled.shape[0]:
        km = KMeans(n_clusters=support, n_init=5, random_state=seed, max_iter=100)
        return km.fit(pooled).cluster_centers_.astype(float)
    idx = rng.choice(pooled.shape[0], size=support, replace=(support > pooled.shape[0]))
    return pooled[idx].copy()


def objective(X: np.ndarray, measures: list[np.ndarray]) -> float:
    b = ot.unif(X.shape[0])
    vals = []
    for Y in measures:
        a = ot.unif(Y.shape[0])
        vals.append(ot.emd2(b, a, pairwise_sqeuclidean(X, Y)))
    return float(np.mean(vals))


def semidiscrete_w2(X: np.ndarray, bary_mean: np.ndarray, bary_cov: np.ndarray,
                    seed: int, n_truth: int = 1000) -> float:
    rng = np.random.default_rng(seed)
    truth = rng.multivariate_normal(bary_mean, bary_cov, size=n_truth)
    return float(np.sqrt(ot.emd2(ot.unif(X.shape[0]), ot.unif(n_truth), pairwise_sqeuclidean(X, truth))))


def run_config(task_id: int, cfg: Config, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    result_npz = out_dir / f"result_{task_id:05d}.npz"
    result_csv = out_dir / f"result_{task_id:05d}.csv"
    if result_npz.exists() and result_csv.exists():
        return

    scale = experiment_scale()
    n_per_measure = 80 if is_smoke() else (500 if scale == "heavy" else 300)
    max_iter = 5 if is_smoke() else (100 if scale == "heavy" else 50)
    n_truth = 100 if is_smoke() else (1000 if scale == "heavy" else 500)

    measures, means, covs, bary_mean, bary_cov = make_gaussian_problem(cfg.rep, n_per_measure=n_per_measure)
    measure_weights = [ot.unif(Y.shape[0]) for Y in measures]
    b = ot.unif(cfg.support)
    X_init = initial_support(measures, cfg.support, seed=800000 + 1000 * cfg.rep + cfg.support)

    t0 = time.perf_counter()
    status = "ok"
    try:
        if cfg.method == "pot_exact":
            X = ot.lp.free_support_barycenter(
                measures, measure_weights, X_init, b,
                weights=None, numItermax=max_iter, stopThr=1e-7, verbose=False, log=False
            )
        elif cfg.method == "pot_sinkhorn":
            X = ot.bregman.free_support_sinkhorn_barycenter(
                measures, measure_weights, X_init, cfg.reg, b,
                weights=None, numItermax=max_iter, stopThr=1e-7, verbose=False, log=False
            )
        else:  # pragma: no cover
            raise ValueError(f"Unknown method {cfg.method}")
    except Exception as exc:  # save failures explicitly for honest reporting
        status = "failed: " + repr(exc)
        X = np.full_like(X_init, np.nan)
    runtime = time.perf_counter() - t0

    if np.isfinite(X).all():
        obj = objective(X, measures)
        sw2 = semidiscrete_w2(X, bary_mean, bary_cov, seed=900000 + cfg.rep, n_truth=n_truth)
        mean_error = float(np.linalg.norm(X.mean(axis=0) - bary_mean))
        cov_error = bures(np.cov(X.T), bary_cov)
    else:
        obj = sw2 = mean_error = cov_error = math.nan

    np.savez_compressed(
        result_npz,
        support=X,
        means=means,
        covs=covs,
        bary_mean=bary_mean,
        bary_cov=bary_cov,
        rep=cfg.rep,
        support_size=cfg.support,
        method=cfg.method,
        reg=cfg.reg,
        runtime_sec=runtime,
        status=status,
        experiment_scale=scale,
        n_per_measure=n_per_measure,
        max_iter=max_iter,
        n_truth=n_truth,
    )
    with result_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "task_id", "rep", "support", "method", "reg", "experiment_scale",
            "n_per_measure", "max_iter", "n_truth", "runtime_sec",
            "objective", "semi_w2", "mean_error", "cov_error", "status"
        ])
        writer.writeheader()
        writer.writerow({
            "task_id": task_id,
            "rep": cfg.rep,
            "support": cfg.support,
            "method": cfg.method,
            "reg": cfg.reg,
            "experiment_scale": scale,
            "n_per_measure": n_per_measure,
            "max_iter": max_iter,
            "n_truth": n_truth,
            "runtime_sec": runtime,
            "objective": obj,
            "semi_w2": sw2,
            "mean_error": mean_error,
            "cov_error": cov_error,
            "status": status,
        })


def main() -> None:
    task_id = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    grid = make_grid()
    if task_id > len(grid):
        return
    script_dir = Path(__file__).resolve().parent
    run_config(task_id, grid[task_id - 1], script_dir / f"runs{run_suffix()}")


if __name__ == "__main__":
    main()
