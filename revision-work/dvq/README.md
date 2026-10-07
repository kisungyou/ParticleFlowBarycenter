# Paired DVQ comparison on identical local summaries

This folder contains the bounded PBMC/FASHION comparison in the second-revision supplement. Twelve configurations combine two datasets, `S=5,10` partitions, and three repetitions, with `k=10` final prototypes and `d=20`. Both aggregators receive the identical saved local-center summaries: `Skd` coordinates, no local cluster-mass counts, and uniform local and outer masses. The comparator fits k-means to the pooled `Sk` centers. Every original observation is assigned to its nearest final center before distortion, ARI, and NMI are evaluated.

The paired distortion changes (DVQ relative to summary k-means) average +1.685% and +0.684% for PBMC at S=5 and S=10, and −3.481% and +0.609% for FASHION. Means and sample SDs over three partitions are descriptive; they do not establish a general ranking.

## Verify retained evidence

From the repository root:

```sh
python3 revision-work/dvq/02_summarize.py --check
Rscript revision-work/dvq/03_verify.R
```

The first command requires only Python's standard library. It rebuilds all means, SDs, and within-partition differences from `paired_metrics.csv` and checks the retained CSV summaries without changing them.

The R command requires `T4transport`, `mclustcomp`, and `jsonlite`. It verifies the data checksums, saved initial supports, final exact OT objectives and feasibility, objective histories, stopping tests, all-observation assignments, distortion, ARI/NMI, and summary budgets. It also checks the historical cap audit and, when the matching legacy files are present, the eight recorded historical matches. It does not recompute local k-means or the full experiment.

All twelve retained DVQ runs satisfy the joint rule within three to five updates. The native-to-controlled support discrepancy is at most 3.638e-12. The retained audit covers 64 available unpaired DVQ outputs, none at their 50-iteration cap; their maximum count is 22. These are different evidence sets and should not be pooled.

## Replay or rebuild

```sh
Rscript revision-work/dvq/01_paired_budget.R
python3 revision-work/dvq/02_summarize.py --input-dir revision-work/dvq/recomputed
```

Default replay loads the bundled local summaries and exact initial support coordinates, runs the explicit full-step MM loop using public `T4transport::wasserstein`, and reruns server-side k-means. It writes to `recomputed/`, preserving all retained evidence. The explicit loop requires both normalized objective change at most 1e-8 and normalized root fixed-point residual at most 1e-6, with at most 100 updates. The scale is the pooled-summary variance. Replay does not depend on the package's private Gaussian initializer. Local-construction time is marked unavailable because replay starts from already prepared summaries; aggregation and evaluation are timed afresh.

A one-configuration replay is available for a quick check:

```sh
Rscript revision-work/dvq/01_paired_budget.R --key pbmc_S5_rep1 --output-dir /tmp/dvq-one-run
```

To rebuild partitions, local centers, and the Gaussian initializer, rather than replay those inputs:

```sh
Rscript revision-work/dvq/01_paired_budget.R --rebuild-inputs --output-dir /tmp/dvq-rebuilt
```

Rebuilding requires the retained `T4transport:::` initializer `aux_ginit`; the script checks for it and stops clearly when it is absent. The recorded package reports version 0.1.9, but version alone does not identify the retained implementation. No installed source commit was available to pin. Default replay avoids this private-API dependency and is the recommended route for checking the controlled algorithm. Package versions and original hardware/threading information are in `metadata.json`.

Partition seeds are `1000000 + 10000*r + 600 + j`, where j=3 for S=5 and j=4 for S=10. Local fit i uses seed+100*i, ten starts, 100 iterations, and at most 5,000 observations. The saved initial atoms record the retained Gaussian initializer at seed+999. Server-side k-means also uses seed+999 and ten starts; the methods share inputs, not a common optimization initializer.

## Inputs, outputs, and missing legacy data

- `inputs/*.rds`: local centers, partition labels, seed, data checksum, and exact saved initial atoms. Paths are repository-relative.
- `runs/*.rds`: controlled/native barycenters, k-means fit, optimization histories, and all-observation assignments and metrics.
- `paired_metrics.csv`, `summary_metrics.csv`, `paired_differences.csv`, `paired_difference_summary.csv`: complete controlled numerical evidence.
- `available_historical_metrics.csv`, `historical_cap_audit.csv`, `reproduction_checks.csv`: the bounded audit of the available legacy records.
- `provenance.json`: original evidence hashes; portable input records add the corresponding saved initial atoms and replace personal paths, without changing numerical arrays.

The existing [`data/processed_pbmc.RData`](../../data/processed_pbmc.RData) and [`data/processed_fmnist.RData`](../../data/processed_fmnist.RData) are reused; no duplicate processed matrices are stored here. The processed `20NEWS` input is unavailable, so no paired NEWS result is claimed or regenerated. The larger historical collection behind the unpaired manuscript tables is incomplete locally; this folder does not promise to reconstruct that full study. Existing scripts and available outputs remain under `code-real-clustering/`.

Original local-construction and aggregation times are retained separately. Other computations were active on the original host, so those timings are descriptive. Replay timings belong to the new machine/run and must not replace the reported historical timings silently.
