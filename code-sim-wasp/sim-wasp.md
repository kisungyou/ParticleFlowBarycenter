# Bayesian posterior aggregation with Wasserstein barycenters
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

This notebook runs the original Bayesian posterior aggregation experiment.
The full-data posterior is available analytically for the simulated
conjugate linear model, while precomputed MCMC/subset posterior draws
are supplied in the data file. The notebook then runs the free-support
barycenter aggregation step and visualizes Wasserstein, moment, conditional-mean, interval-inclusion, and runtime diagnostics. The intervals concern the conditional mean and omit future-observation noise. Their inclusion rates do not establish predictive or repeated-data calibration.

Before rendering this notebook, place the following file in `../data/`:

- `normal-simplified-reduced.RData`

## Setup

Edit only the first cell for a quick smoke run or the full original grid.

``` r
# User settings ------------------------------------------------------------
run_mode <- "smoke"       # "smoke" or "full"
run_compute <- TRUE       # FALSE: skip computation and only make figures from existing results
install_packages <- FALSE # TRUE: install missing R packages from CRAN when possible
jobs <- 1L                # use 1 for smoke; increase for full runs on a workstation
experiment_scale <- "core"
run_label <- ""           # automatic: "smoke" for smoke, "core" for full
# -------------------------------------------------------------------------
```

``` r
source("auxiliary.R")
load_or_install(
  c("T4transport", "mvtnorm", "MASS", "ggplot2", "patchwork", "dplyr", "readr", "tidyr"),
  install = install_packages
)
run_label <- set_run_environment(run_mode, run_label, experiment_scale)
wasp_data <- require_file(file.path("..", "data", "normal-simplified-reduced.RData"), "WASP data file")
Sys.setenv(WASP_DATA_FILE = wasp_data)
message("Run mode: ", run_mode)
message("Run label: ", run_label)
message("Using WASP data: ", wasp_data)
```

## Computation

The full grid varies the number of subsets, barycenter support size, and
step size. The smoke grid runs only a few tasks so that package loading,
data access, T4transport calls, and assembly can be checked quickly.

``` r
if (run_compute) {
  w1_tasks <- if (run_mode == "smoke") 1:2 else 1:54
  run_rscript_grid(
    file.path("simulation-normal-extended", "01_wasp_alpha_barycenters.R"),
    w1_tasks,
    jobs = jobs,
    label = "WASP barycenters"
  )
  run_rscript_once(file.path("simulation-normal-extended", "02_assemble_wasp_alpha_metrics.R"))
}
```

## Figures and tables

The original figures are generated in full mode. Use `../manuscript-figures/build_figures.py` for the current manuscript figures. Smoke mode prints the
small assembled result so the reader can confirm that the pipeline is
working.

``` r
if (run_mode == "full") {
  source("vis-simulation-normal.R")
} else {
  wasp_file <- file.path("simulation-normal-extended", paste0("assembled_wasp_alpha_metrics-", run_label, ".csv"))
  if (file.exists(wasp_file)) print(readr::read_csv(wasp_file, show_col_types = FALSE))
}
```

    # A tibble: 2 × 22
      task_id original_index   rep nsplit support alpha runtime_sec niter
        <dbl>          <dbl> <dbl>  <dbl>   <dbl> <dbl>       <dbl> <dbl>
    1       1            151     1      5      50     1       0.14      5
    2       2            152     2      5      50     1       0.138     5
    # ℹ 14 more variables: history_first <dbl>, history_last <dbl>,
    #   history_min <dbl>, history_length <dbl>, monotone_violations <dbl>,
    #   max_history_increase <dbl>, semi_w2 <dbl>, mean_error <dbl>,
    #   cov_error <dbl>, marginal_coverage <dbl>, marginal_width <dbl>,
    #   pred_rmse_mean <dbl>, pred_interval_coverage <dbl>,
    #   pred_interval_width <dbl>

## Notes

The original CSV field names `pred_rmse_mean`, `pred_interval_coverage`, and `marginal_coverage` are preserved for compatibility. They denote conditional-mean RMSE, conditional-mean inclusion, and coefficient inclusion, respectively.

The precomputed MCMC file is used so that the notebook focuses on the
barycenter aggregation step. Regenerating all MCMC draws is much more
expensive and is not needed to check the barycenter computation.
