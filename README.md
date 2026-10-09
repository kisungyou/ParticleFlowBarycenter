# Free-support Wasserstein barycenter replication

Replication materials for **A Particle-Flow Algorithm for Free-Support Wasserstein Barycenters**. The current experiments use a classical full-step free-support update, its damped form, and an MM interpretation. The repository separates experiments checked with objective and residual criteria from earlier experiments using iteration caps and native stopping rules. It does not treat a capped output as verified convergence.

## Start with the current results

The saved inputs, numerical outputs, and final figure assets are included. Expensive experiments do not need to be rerun to inspect the reported results or rebuild the figures.

```bash
# Standard-library checks of archived outputs, pairings, and file integrity
python3 verify_replication.py

# Install the tested Python dependencies, then rebuild all ten figures
python3 -m pip install -r requirements.txt
python3 manuscript-figures/build_figures.py
```

Figure generation also requires `pdflatex` and the LaTeX packages listed in the [figure guide](manuscript-figures/README.md). It produces PDF and PGF figures with manuscript typography and line/marker distinctions suitable for black-and-white printing. Original manuscript exports are in `manuscript-figures/published/`. Rebuilds go to `manuscript-figures/build/`.

The [manuscript map](MANUSCRIPT.md) identifies the inputs and scripts for every figure and table, and records the limits of the available archive.

## Current experiment workflows

| Directory | Contents |
| --- | --- |
| [revision-work/refresh](revision-work/refresh/README.md) | Gaussian resolution, WASP aggregation, and sensitivity to posterior sample size. Saved supports, input records, and joint stopping diagnostics. |
| [revision-work/gaussian](revision-work/gaussian/README.md) | Paired exact R/POT and stabilized Sinkhorn comparison on shared inputs, weights, initialization, and oracle draws, plus initialization sensitivity. |
| [revision-work/mnist](revision-work/mnist/README.md) | Candidate-restricted point-cloud medoid on the same three balanced MNIST splits, with distance validation and classification results. |
| [revision-work/dvq](revision-work/dvq/README.md) | DVQ and pooled-center k-means on identical PBMC/FASHION summaries and equal summary budgets. |
| [manuscript-figures](manuscript-figures/README.md) | Current PDF/PGF assets, plotting data, typography settings, and rebuilding scripts. |

Each experiment guide gives commands for checking saved results, rebuilding summaries, and recomputing experiments where supported. Recomputed timings depend on the machine and are not expected to equal archived timings.

## Software and reproducibility

The archived controlled computations used R 4.5.1, T4transport 0.1.9, Python 3.11.4, NumPy 2.4.6, SciPy 1.17.1, and POT 0.9.7.post1. The figure dependencies are also pinned in `requirements.txt`. Individual experiment folders preserve additional software and hardware records.

Install the R dependencies with:

```bash
Rscript install-packages.R
```

The recorded T4transport installation has no source-commit identifier. Its native damped update and initialization behavior therefore cannot be identified by version alone. Guides state any required API features, and controlled replay uses saved inputs and initialization where available. Installing the latest package does not establish exact reproduction of the archived native solver.

## Interpretation of the reported results

- The Gaussian resolution and WASP runs use objective and residual checks. The paired Gaussian solver comparison retains iteration-limit outputs and reports their statuses.
- WASP intervals concern the conditional mean or coefficients. They omit future-observation noise, and their inclusion rates do not establish predictive or repeated-data calibration. The three posterior realizations share one simulated dataset.
- MNIST barycentric prototypes use a native 50-iteration cap. Their joint-rule convergence was not verified. The restricted medoid searches ten candidates per class and does not match the barycenter's atom budget.
- Original DVQ comparisons use method-specific seeds. The added PBMC/FASHION comparison pairs identical local-center summaries. Close near-zero scores for PBMC and 20NEWS do not establish effective label recovery.
- Complete numeric plot data were not recovered for several older displays. Those figures retain the original plotted layers and disclose this in the figure guide. They are not presented as regenerated observations.
- The processed 20NEWS matrix is absent. Full recomputation of the original three-dataset study requires that missing input. See [data/README.md](data/README.md).

## Original notebooks

The original Quarto workflows remain available for running the earlier native-stopping experiments:

- [Gaussian](code-sim-gaussian/sim-gaussian.md)
- [WASP](code-sim-wasp/sim-wasp.md)
- [MNIST](code-real-digits/real-digits.md)
- [Clustering](code-real-clustering/real-clustering.md)

They default to small smoke runs. Their `full` mode runs the earlier grids and does not substitute for the current workflows above. The saved smoke results are examples, not the full evidence for the manuscript.

```bash
quarto render code-sim-gaussian/sim-gaussian.qmd
quarto render code-sim-wasp/sim-wasp.qmd
quarto render code-real-digits/real-digits.qmd
quarto render code-real-clustering/real-clustering.qmd
```
