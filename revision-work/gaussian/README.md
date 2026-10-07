# Paired Gaussian comparison

These are the inputs and numerical records used by the current manuscript's
paired Gaussian comparison (Figure 3) and its supplementary multistart and
oracle-sample checks. Publication plots are built by `../../manuscript-figures/`.
The separate resolution/damping experiment is in `../refresh/`.

From the repository root:

```sh
python revision-work/gaussian/paired_benchmark.py --check
python revision-work/gaussian/summarize.py
```

The first command checks all 150 paired records, all 15 multistart records,
input checksums, marginal feasibility, exact-method agreement, and three
accuracy calculations from saved supports. The second rebuilds the numerical
supplement tables, summary JSON, and the oracle-sample sensitivity CSV. Neither
command optimizes a barycenter. The scripts locate their inputs relative to
their own files, so they also work from another working directory.

For more extensive reproduction:

```sh
python revision-work/gaussian/paired_benchmark.py --assemble-only
python revision-work/gaussian/multistart.py
python revision-work/gaussian/paired_benchmark.py --run
python revision-work/gaussian/multistart.py --run
```

The first two commands rescore saved supports. `--run` computes missing fits
and reuses complete records. To recompute a particular fit, move its final CSV,
support CSV, and history CSV out of `runs/` before using `--run`; keep the
canonical input NPZ. `REPS` and `SUPPORTS` control an exploratory run, with
defaults `10` and `50,100,200`. Changing them also changes the assembled summary.
Do not mix an exploratory summary with the published figure.

## Dependencies and numerical protocol

Python needs NumPy, SciPy, pandas, POT, scikit-learn, and threadpoolctl. An R
rerun needs `Rscript` and T4transport. The recorded environment is in
`metadata.json` and `R-session.txt`; these retain software and hardware details
without personal filesystem paths. The portable R implementation uses the
public `T4transport::wasserstein` plan. The original transparent R loop called
the exact-plan routine directly; a representative public-interface rerun
reproduces its saved support within 3.2e-13 per coordinate. Fresh wall-clock times need not match
the archived measurements.

There are ten independent Gaussian problems, with four inputs containing 300
uniformly weighted observations each. Barycenter sizes are 50, 100, and 200;
outer weights are 1/4. Each case saves explicit pooled-k-means initialization,
all coordinates and masses, and oracle draws before any solver runs. R and
Python therefore receive the same problem. Forty-five `problem.npz` files
include the 30 main cases and 15 multistart cases. Redundant R input CSVs are
regenerated from NPZ only when an R computation is requested.

The exact implementations apply the same full barycentric-projection update.
They are an implementation comparison, not competing update algorithms.
Stabilized Sinkhorn uses epsilon/s² = 0.1, 0.3, and 1, where s² is the average
within-input squared distance from each input's own mean. Inner stopping uses
tolerance 1e-9 and a cap of 20,000 iterations. Fixed row/column cost offsets
improve numerical scaling without changing the coupling; objectives use the
original costs. The outer limit is 300 updates. Termination requires both a
relative objective change no greater than 1e-8, with denominator
`max(1, abs(previous), abs(current))`, and
`sum_i b_i ||M_i-Z_i||^2 / s² <= 1e-7`, checked at the returned support.
Regularized fits use their own objective and projection for stopping. Their
exact-OT residual is also recorded but need not vanish.

Every output is scored using the same unregularized objective and 1,000 common
oracle draws per problem. Atomic covariance uses the prescribed masses with
no sample-size correction. Figure summaries are means and sample SDs over ten
problems. Capped fits remain in the summaries: 0, 3, and 4 cap exits among 30
fits at the three regularization levels. All retained inner solves converged.
The exact methods agree to numerical precision. Runtime includes optimization
and convergence checks, excluding initialization and external accuracy scoring.

The multistart study uses five starts on each of three fixed problems at m=50.
Its reference is the best objective attained among those starts, not a known
global optimum. The oracle check compares nested sample sizes 500, 1,000, and
2,000 for three fixed exact fits at m=100.

## Retained files

- `paired_results.csv`, `multistart_results.csv`, `oracle_sensitivity.csv`:
  publication data for plots and supplementary checks.
- `runs/`: fitted supports, complete histories, and final records.
- `inputs/`, `input-checksums.sha256`: canonical serialized cases.
- `paired_diagnostics_table.tex`, `multistart_table.tex`, `summary.json`:
  generated numerical summaries.
- `validation.json`: results of the lightweight saved-output checks.

Superseded pilot outputs, temporary partial histories, font caches, and draft
manuscript prose are intentionally excluded. The exploratory regularization
ratio 0.03 was rejected after inaccurate inner solves. An initial inner cap of
4,000 was increased after one initialization reached it; that fit was rerun,
and all retained fits satisfy the final inner-solver checks.
