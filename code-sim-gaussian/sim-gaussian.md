# Gaussian barycenters: resolution, step size, and baselines
Kisung You

- [Overview](#overview)
- [Setup](#setup)
- [Computation](#computation)
- [Figures and tables](#figures-and-tables)
- [Notes](#notes)

## Overview

This notebook reproduces the Gaussian barycenter experiment from the
paper. The purpose is to study finite-resolution behavior in a setting
where the population barycenter is known, to compare the full update
with a damped update, and to give representative computational baselines
using POT.

The notebook is self-contained except for software dependencies. The
Gaussian inputs are generated synthetically, so no external data file is
required.

## Setup

The first code cell is the only place where most users should edit the
run. Use `run_mode <- "smoke"` for a quick check and
`run_mode <- "full"` for the manuscript-scale reviewer-response grid.
Full mode can be computationally expensive because exact OT subproblems
dominate runtime.

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
  c("T4transport", "mvtnorm", "abind", "ggplot2", "patchwork", "dplyr", "readr", "tidyr"),
  install = install_packages
)
run_label <- set_run_environment(run_mode, run_label, experiment_scale)
message("Run mode: ", run_mode)
message("Run label: ", run_label)
```

## Computation

The finite-resolution experiment uses support sizes and step sizes from
the paper. The runtime experiment uses one-factor sweeps in support
size, target support size, number of measures, dimension, and step size.
The POT baseline requires Python with the packages listed in
`../requirements.txt`.

``` r
if (run_compute) {
  if (run_mode == "smoke") {
    # Enough to test data generation, T4transport calls, saving, and assembly.
    g1_tasks <- 1:2
    g2_tasks <- 1:2
    pot_tasks <- 1
  } else {
    # Core reviewer-response grid after reducing repetitions to keep exact OT feasible.
    g1_tasks <- 1:36
    g2_tasks <- 1:54
    pot_tasks <- 1:27
  }

  run_rscript_grid(
    file.path("simulation-gaussian-extended", "01_gaussian_resolution_alpha.R"),
    g1_tasks,
    jobs = jobs,
    label = "Gaussian resolution / alpha"
  )
  run_rscript_once(file.path("simulation-gaussian-extended", "02_assemble_gaussian_resolution_alpha.R"))

  run_rscript_grid(
    file.path("simulation-gaussian-extended", "03_runtime_scaling.R"),
    g2_tasks,
    jobs = jobs,
    label = "Gaussian runtime scaling"
  )
  run_rscript_once(file.path("simulation-gaussian-extended", "04_assemble_runtime_scaling.R"))

  # POT is optional in smoke mode if Python/POT is unavailable. In full mode,
  # it is part of the baseline comparison used in the paper.
  try({
    run_python_grid(
      file.path("pot-baselines", "01_pot_gaussian_baselines.py"),
      pot_tasks,
      jobs = jobs,
      label = "POT baselines"
    )
    run_python_grid(
      file.path("pot-baselines", "02_assemble_pot_gaussian_baselines.py"),
      1,
      jobs = 1,
      label = "POT assembly"
    )
  }, silent = (run_mode == "smoke"))
}
```


    Python task failed: /Users/kyou/Desktop/ParticleFlowBarycenter/code-sim-gaussian/pot-baselines/01_pot_gaussian_baselines.py id=1
    Log file: /Users/kyou/Desktop/ParticleFlowBarycenter/code-sim-gaussian/pot-baselines/logs/01_pot_gaussian_baselines_00001.log
    Missing Python package 'POT' (import name: ot). Install with: pip install POT 

## Figures and tables

In full mode, the following visualization script reproduces the
manuscript figures and summary tables. In smoke mode, the same script
may have too little data for all panels, so we print the assembled smoke
summaries instead.

``` r
if (run_mode == "full") {
  source("vis-simulation-gaussian.R")
} else {
  g1_file <- file.path("simulation-gaussian-extended", paste0("assembled_gaussian_resolution_alpha-", run_label, ".csv"))
  g2_file <- file.path("simulation-gaussian-extended", paste0("assembled_runtime_scaling-", run_label, ".csv"))
  if (file.exists(g1_file)) print(readr::read_csv(g1_file, show_col_types = FALSE))
  if (file.exists(g2_file)) print(readr::read_csv(g2_file, show_col_types = FALSE))
}
```

    # A tibble: 2 × 16
      task_id   rep support alpha runtime_sec niter history_first history_last
        <dbl> <dbl>   <dbl> <dbl>       <dbl> <dbl>         <dbl>        <dbl>
    1       1     1      10     1      0.0350     5          215.         213.
    2       2     2      10     1      0.034      5          210.         209.
    # ℹ 8 more variables: history_min <dbl>, history_length <dbl>,
    #   monotone_violations <dbl>, max_history_increase <dbl>,
    #   final_objective <dbl>, semi_w2 <dbl>, mean_error <dbl>, cov_error <dbl>
    # A tibble: 2 × 17
      task_id sweep   rep     N target_n     d support alpha runtime_sec niter
        <dbl> <chr> <dbl> <dbl>    <dbl> <dbl>   <dbl> <dbl>       <dbl> <dbl>
    1       1 m         1     4       80     2      10     1       0.032     5
    2       2 m         1     4       80     2      10     1       0.031     5
    # ℹ 7 more variables: history_first <dbl>, history_last <dbl>,
    #   history_min <dbl>, history_length <dbl>, monotone_violations <dbl>,
    #   max_history_increase <dbl>, final_objective <dbl>

## Notes

The full update `alpha = 1` is the default because it exactly minimizes
the frozen-plan quadratic surrogate when OT subproblems are solved
exactly. The damped update `alpha = 0.5` is included as a sensitivity
analysis.
