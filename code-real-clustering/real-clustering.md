# Distributed vector quantization with barycentric prototypes
Kisung You

- [Overview](#overview)
- [Setup](#setup)
- [Computation](#computation)
- [Figures and tables](#figures-and-tables)
- [Notes](#notes)

## Overview

> **Original workflow.** This notebook uses the original iteration caps and native
> stopping criteria. The [current replication guide](../README.md) points to the
> checked experiments, saved results, and final PDF/PGF figures. Full mode here
> runs the original workflow and does not reproduce the added paired comparisons
> or certify convergence under the joint objective/residual rule.

This notebook runs the original distributed vector quantization (DVQ) experiment. The goal is to evaluate whether free-support Wasserstein
barycenters can be used as compressed clustering prototypes. The
notebook compares DVQ with sampled-fit `k`-means and spherical `k`-means.

The manuscript includes PBMC, 20NEWS, and FASHION. The processed 20NEWS matrix is not available in this archive. This notebook therefore runs PBMC and FASHION by default. The current paired comparison likewise includes only these two datasets.

For the default two-dataset run, place the following processed files in
`../data/`:

- `processed_pbmc.RData`
- `processed_fmnist.RData`

Each file should contain a processed feature matrix `X_reduced` and a
label vector `Y`.

## Setup

Edit the first cell to choose a smoke test or the full original grid.

``` r
# User settings ------------------------------------------------------------
run_mode <- "smoke"       # "smoke" or "full"
run_compute <- TRUE       # FALSE: skip computation and only make figures from existing results
install_packages <- FALSE # TRUE: install missing R packages from CRAN when possible
jobs <- 1L                # use 1 for smoke; increase for full runs on a workstation
experiment_scale <- "core"
run_label <- ""           # automatic: "smoke-nonews" or "core-nonews"

# Data sets to include. The processed 20NEWS matrix is not bundled.
# To run the original three-dataset workflow after obtaining that matrix, use
# c("pbmc", "news", "fashion") and put processed_news.RData in ../data/.
cluster_datasets <- c("pbmc", "fashion")
# -------------------------------------------------------------------------
```

``` r
source("auxiliary.R")
load_or_install(
  c("T4transport", "skmeans", "cluster", "clusterCrit", "mclustcomp", "ggplot2", "patchwork", "dplyr", "readr", "tidyr"),
  install = install_packages
)
cluster_datasets <- match.arg(cluster_datasets, choices = c("pbmc", "news", "fashion"), several.ok = TRUE)
cluster_datasets <- unique(cluster_datasets)

# Use a distinct label for the default two-dataset run so that outputs are not
# mixed with older three-dataset results.
if (!nzchar(run_label) && !all(c("pbmc", "news", "fashion") %in% cluster_datasets)) {
  run_label <- if (run_mode == "smoke") "smoke-nonews" else "core-nonews"
}
run_label <- set_run_environment(run_mode, run_label, experiment_scale)

data_dir <- normalizePath(file.path("..", "data"), mustWork = TRUE)
required_data_files <- c(
  pbmc = "processed_pbmc.RData",
  news = "processed_news.RData",
  fashion = "processed_fmnist.RData"
)
for (dataset_name in cluster_datasets) {
  require_file(
    file.path(data_dir, required_data_files[[dataset_name]]),
    paste0(toupper(dataset_name), " processed data")
  )
}
Sys.setenv(
  CLUSTER_DATA_DIR = data_dir,
  CLUSTER_DATASETS = paste(cluster_datasets, collapse = ",")
)
message("Run mode: ", run_mode)
message("Run label: ", run_label)
message("Data sets: ", paste(cluster_datasets, collapse = ", "))
message("Using clustering data directory: ", data_dir)
```

## Computation

The full grid varies the selected data sets, cluster count, DVQ subset
count, and step size. The smoke grid runs a single task to test data
loading, clustering metrics, saving, and assembly.

``` r
if (run_compute) {
  if (run_mode == "smoke") {
    c1_tasks <- 1L
  } else {
    n_methods <- if (tolower(experiment_scale) == "heavy") 8L else 6L
    n_cluster_indices <- if (tolower(experiment_scale) == "heavy") 11L else 3L
    n_repetitions <- if (tolower(experiment_scale) == "heavy") 10L else 3L
    n_tasks <- length(cluster_datasets) * n_methods * n_cluster_indices * n_repetitions
    c1_tasks <- seq_len(n_tasks)
  }
  run_rscript_grid(
    file.path("real-2-clustering-extended", "01_clustering_alpha_runtime.R"),
    c1_tasks,
    jobs = jobs,
    label = "DVQ clustering"
  )
  run_rscript_once(file.path("real-2-clustering-extended", "02_assemble_clustering_alpha_runtime.R"))
}
```

## Figures and tables

In full mode, the visualization script produces the PCA visualization,
clustering quality figure, and summary tables for the selected data
sets. In smoke mode we print the assembled result.

``` r
if (run_mode == "full") {
  source("vis-real-clustering.R")
} else {
  cl_file <- file.path("real-2-clustering-extended", paste0("assembled_cluster_extended-", run_label, ".csv"))
  if (file.exists(cl_file)) print(readr::read_csv(cl_file, show_col_types = FALSE))
}
```

    # A tibble: 1 × 21
      task_id data  algorithm base_algorithm split alpha cluster_index nclust
        <dbl> <chr> <chr>     <chr>          <lgl> <lgl>         <dbl>  <dbl>
    1       1 pbmc  kmeans    kmeans         NA    NA                4      8
    # ℹ 13 more variables: repeat_id <dbl>, runtime_sec <dbl>, niter <lgl>,
    #   history_first <lgl>, history_last <lgl>, history_min <lgl>,
    #   history_length <dbl>, monotone_violations <lgl>,
    #   max_history_increase <lgl>, ARI <dbl>, NMI <dbl>, Silhouette <dbl>,
    #   CH <dbl>

## Notes

The reference labels are used only for evaluation through ARI and NMI.
The Silhouette and Calinski–Harabasz indices are internal geometric
criteria. The original method-specific seeds make the comparisons descriptive rather than paired. The current matched-summary comparison is in `../revision-work/dvq`. Running 20NEWS also requires its missing processed matrix.
