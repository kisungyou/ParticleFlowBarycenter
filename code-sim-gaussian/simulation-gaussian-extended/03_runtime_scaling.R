# Revision Experiment G2: Synthetic runtime scaling.
# Purpose: answer reviewer requests for runtime/memory scaling versus support size m,
# number of input measures N, target support size n, ambient dimension d, and alpha.
#
# Usage: Rscript 03_runtime_scaling.R <task_id>

script_dir <- getwd()
source(file.path(dirname(script_dir), "common", "rev_utils.R"))
rev_require(c("T4transport", "mvtnorm", "abind"))

myid <- rev_arg_id(default = 1L)
out_dir <- rev_output_dir(script_dir, "runs-scaling")
rev_ensure_dir(out_dir)

# One-factor sweeps. The default core grid is much smaller than the original
# heavy grid, but still covers scaling in support size, target support size,
# number of input measures, ambient dimension, and step size.
rev_scale <- tolower(Sys.getenv("REV_EXPERIMENT_SCALE", unset = "core"))

if (rev_scale == "heavy") {
  vec_rep <- seq_len(10L)
  config_m <- data.frame(sweep = "m", rep = rep(vec_rep, each = 7),
                         N = 4L, n = 1000L, d = 2L,
                         m = rep(c(10L, 25L, 50L, 100L, 200L, 500L, 1000L), times = length(vec_rep)),
                         alpha = 1.0)
  config_n <- data.frame(sweep = "target_n", rep = rep(vec_rep, each = 5),
                         N = 4L, n = rep(c(50L, 100L, 200L, 500L, 1000L), times = length(vec_rep)),
                         d = 2L, m = 100L, alpha = 1.0)
  config_N <- data.frame(sweep = "N", rep = rep(vec_rep, each = 5),
                         N = rep(c(2L, 4L, 8L, 16L, 32L), times = length(vec_rep)),
                         n = 200L, d = 2L, m = 100L, alpha = 1.0)
  config_d <- data.frame(sweep = "d", rep = rep(vec_rep, each = 6),
                         N = 4L, n = 200L,
                         d = rep(c(2L, 5L, 10L, 25L, 50L, 100L), times = length(vec_rep)),
                         m = 100L, alpha = 1.0)
  config_alpha <- data.frame(sweep = "alpha", rep = rep(vec_rep, each = 4),
                             N = 4L, n = 500L, d = 2L, m = 200L,
                             alpha = rep(c(1.0, 0.75, 0.50, 0.25), times = length(vec_rep)))
} else {
  vec_rep <- seq_len(3L)
  config_m <- data.frame(sweep = "m", rep = rep(vec_rep, each = 5),
                         N = 4L, n = 300L, d = 2L,
                         m = rep(c(10L, 50L, 100L, 200L, 500L), times = length(vec_rep)),
                         alpha = 1.0)
  config_n <- data.frame(sweep = "target_n", rep = rep(vec_rep, each = 4),
                         N = 4L, n = rep(c(50L, 100L, 200L, 500L), times = length(vec_rep)),
                         d = 2L, m = 100L, alpha = 1.0)
  config_N <- data.frame(sweep = "N", rep = rep(vec_rep, each = 4),
                         N = rep(c(2L, 4L, 8L, 16L), times = length(vec_rep)),
                         n = 200L, d = 2L, m = 100L, alpha = 1.0)
  config_d <- data.frame(sweep = "d", rep = rep(vec_rep, each = 3),
                         N = 4L, n = 200L,
                         d = rep(c(2L, 10L, 50L), times = length(vec_rep)),
                         m = 100L, alpha = 1.0)
  config_alpha <- data.frame(sweep = "alpha", rep = rep(vec_rep, each = 2),
                             N = 4L, n = 300L, d = 2L, m = 200L,
                             alpha = rep(c(1.0, 0.50), times = length(vec_rep)))
}
config_grid <- rbind(config_m, config_n, config_N, config_d, config_alpha)
config_grid$task_id <- seq_len(nrow(config_grid))

if (myid > nrow(config_grid)) quit(save = "no", status = 0, runLast = FALSE)
now <- config_grid[myid, ]

save_file <- file.path(out_dir, sprintf("result_%05d.RData", myid))
if (file.exists(save_file)) quit(save = "no", status = 0, runLast = FALSE)

maxiter <- if (rev_scale == "heavy") rev_smoke_int(50L, 5L) else rev_smoke_int(20L, 5L)
abstol <- 1e-8

now_N <- if (rev_is_smoke()) min(as.integer(now$N), 4L) else as.integer(now$N)
now_n <- if (rev_is_smoke()) min(as.integer(now$n), 80L) else as.integer(now$n)
now_d <- if (rev_is_smoke()) min(as.integer(now$d), 5L) else as.integer(now$d)
now_m <- if (rev_is_smoke()) min(as.integer(now$m), 10L) else as.integer(now$m)

# Generate random Gaussian measures. The seed is tied to the full configuration
# to make reruns reproducible while keeping each task independent. Smoke mode
# uses a smaller problem with the same code path and separate output directory.
data_seed <- 400000L + 10000L * now$rep + myid
dat <- rev_make_gaussian_measures(seed = data_seed,
                                  num_measures = now_N,
                                  n_per_measure = now_n,
                                  d = now_d,
                                  loc_scale = 5,
                                  wishart_df = max(now_d + 2L, 4L),
                                  loc_sd = 0.5)

init_seed <- 500000L + 10000L * now$rep + now_m + as.integer(100 * now$alpha)
set.seed(init_seed)
fit <- rev_elapsed(
  T4transport::rbaryGD(dat$measures,
                       num_support = now_m,
                       alpha = as.numeric(now$alpha),
                       maxiter = maxiter,
                       abstol = abstol)
)
est_bary <- fit$value
hsum <- rev_history_summary(est_bary$history)
final_objective <- rev_empirical_objective(est_bary$support, dat$measures)

metrics <- data.frame(
  task_id = myid,
  sweep = as.character(now$sweep),
  rep = as.integer(now$rep),
  N = now_N,
  target_n = now_n,
  d = now_d,
  support = now_m,
  alpha = as.numeric(now$alpha),
  runtime_sec = fit$elapsed,
  niter = est_bary$niter,
  history_first = hsum$first,
  history_last = hsum$last,
  history_min = hsum$min,
  history_length = hsum$n,
  monotone_violations = hsum$monotone_violations,
  max_history_increase = hsum$max_increase,
  final_objective = final_objective
)
metadata <- list(experiment = "G2_runtime_scaling",
                 maxiter = maxiter, abstol = abstol,
                 experiment_scale = rev_scale,
                 data_seed = data_seed, init_seed = init_seed)

save(metadata, metrics, est_bary, file = save_file)
