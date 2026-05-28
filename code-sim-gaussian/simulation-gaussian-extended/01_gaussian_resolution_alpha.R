# Revision Experiment G1: Gaussian accuracy, resolution, step-size, and runtime.
# Purpose: answer reviewer concerns about finite resolution (m), alpha sensitivity,
# monotonicity diagnostics, and runtime in a setting with an analytic Gaussian barycenter.
#
# Usage: Rscript 01_gaussian_resolution_alpha.R <task_id>

script_dir <- getwd()
source(file.path(dirname(script_dir), "common", "rev_utils.R"))

rev_require(c("T4transport", "mvtnorm", "abind"))

myid <- rev_arg_id(default = 1L)
out_dir <- rev_output_dir(script_dir, "runs")
rev_ensure_dir(out_dir)

# Two grid sizes are supported. The default "core" grid is the reviewer-response
# grid: it keeps large-m evidence and alpha sensitivity while avoiding a
# prohibitively expensive full factorial exact-OT study. Set
# REV_EXPERIMENT_SCALE=heavy to recover the original large grid.
rev_scale <- tolower(Sys.getenv("REV_EXPERIMENT_SCALE", unset = "core"))

if (rev_scale == "heavy") {
  vec_supp <- c(10L, 25L, 50L, 100L, 200L, 500L, 1000L)
  vec_alpha <- c(1.00, 0.75, 0.50, 0.25)
  vec_rep <- seq_len(50L)
} else {
  vec_supp <- c(10L, 25L, 50L, 100L, 200L, 500L)
  vec_alpha <- c(1.00, 0.50)
  vec_rep <- seq_len(3L)
}

config_grid <- expand.grid(rep = vec_rep, supp = vec_supp, alpha = vec_alpha)
if (myid > nrow(config_grid)) quit(save = "no", status = 0, runLast = FALSE)

now_rep <- as.integer(config_grid$rep[myid])
now_supp <- as.integer(config_grid$supp[myid])
now_alpha <- as.numeric(config_grid$alpha[myid])

save_file <- file.path(out_dir, sprintf("result_%05d.RData", myid))
if (file.exists(save_file)) quit(save = "no", status = 0, runLast = FALSE)

# Constants used for all tasks. The core grid is intentionally lighter than the
# original heavy grid but still uses the same code path.
num_measures <- 4L
if (rev_scale == "heavy") {
  n_per_measure <- rev_smoke_int(1000L, 80L)
  oracle_sample_n <- rev_smoke_int(2000L, 100L)
  maxiter <- rev_smoke_int(100L, 5L)
} else {
  n_per_measure <- rev_smoke_int(300L, 80L)
  oracle_sample_n <- rev_smoke_int(500L, 100L)
  maxiter <- rev_smoke_int(50L, 5L)
}
abstol <- 1e-8

# Generate the same data for all alpha/support settings within a repetition.
data_seed <- 100000L + now_rep
dat <- rev_make_gaussian_measures(seed = data_seed,
                                  num_measures = num_measures,
                                  n_per_measure = n_per_measure,
                                  d = 2L,
                                  loc_scale = 10,
                                  wishart_df = 4L,
                                  loc_sd = 1)

est_gauss <- T4transport::gaussbarypd(means = dat$means, vars = dat$covs)

# Force the same randomized k-means initialization across alpha values for the
# same repetition/support size, if the initializer uses R's RNG.
init_seed <- 200000L + 1000L * now_rep + now_supp
set.seed(init_seed)
fit <- rev_elapsed(
  T4transport::rbaryGD(dat$measures,
                       num_support = now_supp,
                       alpha = now_alpha,
                       maxiter = maxiter,
                       abstol = abstol)
)
est_bary <- fit$value
runtime_sec <- fit$elapsed

# Accuracy relative to the analytic Gaussian barycenter.
set.seed(300000L + now_rep)
oracle_draws <- mvtnorm::rmvnorm(n = oracle_sample_n,
                                 mean = est_gauss$mean,
                                 sigma = est_gauss$var)
semi_w2 <- T4transport::wasserstein(est_bary$support, oracle_draws, p = 2)$distance
mean_error <- sqrt(sum((colMeans(est_bary$support) - est_gauss$mean)^2))
cov_error <- rev_bures(stats::cov(est_bary$support), est_gauss$var)
final_objective <- rev_empirical_objective(est_bary$support, dat$measures)
hsum <- rev_history_summary(est_bary$history)

metadata <- list(
  experiment = "G1_gaussian_resolution_alpha",
  task_id = myid,
  rep = now_rep,
  support = now_supp,
  alpha = now_alpha,
  n_per_measure = n_per_measure,
  oracle_sample_n = oracle_sample_n,
  maxiter = maxiter,
  abstol = abstol,
  experiment_scale = rev_scale,
  data_seed = data_seed,
  init_seed = init_seed
)

metrics <- data.frame(
  task_id = myid,
  rep = now_rep,
  support = now_supp,
  alpha = now_alpha,
  runtime_sec = runtime_sec,
  niter = est_bary$niter,
  history_first = hsum$first,
  history_last = hsum$last,
  history_min = hsum$min,
  history_length = hsum$n,
  monotone_violations = hsum$monotone_violations,
  max_history_increase = hsum$max_increase,
  final_objective = final_objective,
  semi_w2 = semi_w2,
  mean_error = mean_error,
  cov_error = cov_error
)

save(metadata, metrics, dat, est_gauss, est_bary, file = save_file)
