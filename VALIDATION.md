# Repository alignment checks

Checked against the current manuscript and supplement on 2026-10-07.

## Numerical evidence

- All four main numerical tables agree with the repository summaries at their printed precision, including the added medoid row.
- The 36 Gaussian resolution fits, 54 WASP fits, and 12 posterior sample-size fits retain their reported objective and residual checks.
- All 150 paired Gaussian outputs and 15 multistart outputs pass the archived-record checks. Exact R/POT objective differences are at most 5.69e-13. A representative public-R-interface replay agrees with its stored support within 3.20e-13. Seven regularized iteration-limit outputs remain included and labelled.
- Replaying the small Gaussian resolution case reproduces its stored support and 12 iterations exactly. The posterior sensitivity summaries are rescored from the bundled inputs and supports.
- MNIST summary verification reconstructs metrics from all 7,500 saved predictions. Independent R checks confirm split balance, candidate membership, and metrics. Twenty R/POT squared-distance comparisons agree within 1.53e-16.
- The 12 saved DVQ runs pass final OT, feasibility, stopping, assignment, and metric checks. The eight available historical matches agree. A fresh one-configuration replay reproduces the archived DVQ and k-means centers and scores.
- All ten figure builds match the current manuscript: regenerated PGFs and their 31 raster dependencies are byte-identical, and rendered PDFs match at 120 dpi. All 14 optional R plotting exports match the bundled plotting data.

## Repeat the checks

```sh
python3 verify_replication.py
python3 manuscript-figures/build_figures.py --check-published
python3 revision-work/gaussian/paired_benchmark.py --check
Rscript revision-work/refresh/refresh.R --check
Rscript revision-work/refresh/sensitivity.R --check
python3 revision-work/mnist/03_report.py --check
Rscript revision-work/mnist/04_verify.R
python3 revision-work/mnist/02_medoid.py --validate-only
python3 revision-work/dvq/02_summarize.py --check
Rscript revision-work/dvq/03_verify.R
```

The repository-wide check uses only Python's standard library. Experiment checks have the dependencies stated in their directory guides. Full figure generation is separate from numerical recomputation:

```sh
python3 manuscript-figures/build_figures.py
```

Python sources and R scripts pass syntax checks. The current workflows resolve their inputs from the repository, including compact posterior sensitivity inputs and saved initial supports. Published assets and archived numeric inputs are covered by `archive-sha256.json`.

These checks do not rerun all expensive experiments, verify MCMC mixing, recover missing historical inputs, or establish convergence for the old capped prototype runs. The original notebooks retain their historical role. Their descriptions now use the manuscript's conditional-mean interpretation and identify the sampled-fit k-means baseline.
