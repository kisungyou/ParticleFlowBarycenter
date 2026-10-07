# Map to the current manuscript

This map refers to the October 2026 second JCGS revision. It connects the reported displays to the archived computations. The supplement's accompanying `_revision_workwork` records are packaged here as `revision-work/`, with repository-relative paths.

## Main figures

All final PDF/PGF files are in `manuscript-figures/published/`. The figure guide explains the saved-data build and distinguishes numeric data from preserved image layers.

| Figure | File stem | Numerical source or provenance |
| --- | --- | --- |
| 1 | `fig-sim-gauss-2` | `revision-work/refresh/gaussian_resolution.csv` and the 36 saved fits. Both step sizes are present. Near-overlap does not imply identical support configurations. |
| 2 | `fig-sim-gauss-3` | Original finite-budget runtime plot layers. The complete numeric sweep outputs were not recovered. |
| 3 | `fig-sim-gauss-4` | `revision-work/gaussian/paired_results.csv`, shared serialized inputs, and common unregularized scoring. Panels are (A), (B), (C). |
| 4 | `fig-sim-normal-1` | Saved posterior display data, the analytic posterior, and the checked WASP fit. |
| 5 | `fig-sim-normal-2` | `revision-work/refresh/wasp_metrics.csv` and `wasp_full_reference.csv`. |
| 6 | `fig-sim-normal-3` | The same WASP records, with conditional-mean and coefficient-inclusion interpretation. |
| 7 | `fig-real-digits-1` | Saved MNIST image and its Otsu threshold, foreground mask, histogram, and point cloud. |
| 8 | `fig-real-digits-2` | Original prototype and density layers. Complete saved prototype objects were not recovered. |
| 9 | `fig-real-clustering-1` | Numeric PCA data for PBMC/FASHION and the preserved 20NEWS panel. |
| 10 | `fig-real-clustering-2` | Original clustering-summary plot layers. Complete numeric panel data were not recovered. |

## Main tables

The CSVs below are in `manuscript-figures/data/manuscript-tables/`. They preserve reported summaries even when the complete earlier per-run output is unavailable. Such summaries support inspection of the published table, not recovery of missing run-level data.

| Table | Records |
| --- | --- |
| 1, digit classification | `table-real-digits-classification-summary.csv`, `table-real-digits-baseline-summary.csv`, and `revision-work/mnist/split_metrics.csv` for the restricted medoid row. |
| 2, digit diagnostics | `table-real-digits-centroid-summary.csv` and classification timing from `table-real-digits-classification-summary.csv`. |
| 3, central clustering comparisons | `table-real-clustering-central-summary.csv`. These are the original descriptive comparisons. |
| 4, DVQ diagnostics | `table-real-clustering-dvq-diagnostics.csv`. Native iteration histories do not certify the joint stopping rule. |

## Supplement tables

| Table | Records |
| --- | --- |
| 1, Gaussian resolution diagnostics | `revision-work/refresh/gaussian_resolution.csv` and saved Gaussian fit records. |
| 2, paired Gaussian solver diagnostics | `revision-work/gaussian/paired_results.csv` and saved supports/histories. |
| 3, initialization sensitivity | `revision-work/gaussian/multistart_results.csv`. |
| 4, posterior draw-count sensitivity | `revision-work/refresh/draw_sensitivity.csv` and the compact saved sensitivity inputs. |
| 5, WASP diagnostics | `revision-work/refresh/wasp_metrics.csv` and saved WASP fits. |
| 6, restricted point-cloud medoid | `revision-work/mnist/split_metrics.csv`, splits, selected prototypes, and distance validation. |
| 7, paired DVQ differences | `revision-work/dvq/paired_differences.csv` and `paired_difference_summary.csv`, backed by identical saved local-center summaries. |

## Changes from the initial repository

The current workflow adds the shared-input Gaussian comparison, explicit stopping diagnostics, posterior sample-size sensitivity, restricted MNIST medoid, and paired DVQ comparison. It also includes every current figure with PDF/PGF typography, non-color line/marker cues, and clearly documented original plot layers where data are incomplete.

The original notebooks remain usable as historical workflows. Their documentation now separates those computations from the current results, corrects the WASP interpretation, and records the unavailable 20NEWS data without claiming a verified reason for its absence.
