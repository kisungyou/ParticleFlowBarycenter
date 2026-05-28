# Revision Experiment W1: WASP alpha/support sensitivity with runtime histories.
# Purpose: reuse existing MCMC subset posterior draws and recompute barycenters
# with several fixed step sizes alpha. This separates barycenter-solver behavior
# from MCMC sampling and subset generation.
#
# Usage: Rscript 01_wasp_alpha_barycenters.R <task_id>

script_dir <- getwd()
source(file.path(dirname(script_dir), "common", "rev_utils.R"))
rev_require(c("T4transport"))

find_data_file <- function() {
  env <- Sys.getenv("WASP_DATA_FILE", unset = "")
  candidates <- c(
    env,
    file.path(script_dir, "normal-simplified-reduced.RData"),
    file.path(dirname(script_dir), "simulation-normal", "normal-simplified-reduced.RData"),
    file.path(dirname(dirname(script_dir)), "simulation-normal", "normal-simplified-reduced.RData"),
    file.path(getwd(), "normal-simplified-reduced.RData")
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    stop("Could not find normal-simplified-reduced.RData. Set WASP_DATA_FILE=/path/to/normal-simplified-reduced.RData.",
         call. = FALSE)
  }
  hit[1L]
}

myid <- rev_arg_id(default = 1L)
out_dir <- rev_output_dir(script_dir, "runs")
rev_ensure_dir(out_dir)

# Use all 50 original data/MCMC repetitions but selected subset counts. This is
# enough to address split robustness without recomputing the entire MCMC suite.
rev_scale <- tolower(Sys.getenv("REV_EXPERIMENT_SCALE", unset = "core"))
if (rev_scale == "heavy") {
  vec_rep <- seq_len(50L)
  vec_nsplit <- c(2L, 5L, 10L, 20L)
  vec_supp <- c(10L, 25L, 50L, 100L, 200L)
  vec_alpha <- c(1.00, 0.75, 0.50, 0.25)
} else {
  vec_rep <- seq_len(3L)
  vec_nsplit <- c(5L, 10L, 20L)
  vec_supp <- c(50L, 100L, 200L)
  vec_alpha <- c(1.00, 0.50)
}
config_grid <- expand.grid(rep = vec_rep, nsplit = vec_nsplit,
                           support = vec_supp, alpha = vec_alpha)
if (myid > nrow(config_grid)) quit(save = "no", status = 0, runLast = FALSE)
now_rep <- as.integer(config_grid$rep[myid])
now_nsplit <- as.integer(config_grid$nsplit[myid])
now_supp <- as.integer(config_grid$support[myid])
now_alpha <- as.numeric(config_grid$alpha[myid])

save_file <- file.path(out_dir, sprintf("result_%05d.RData", myid))
if (file.exists(save_file)) quit(save = "no", status = 0, runLast = FALSE)

data_file <- find_data_file()
load(data_file)  # common_posterior_mean, common_posterior_cov, common_full_mcmc, individual_nsplits, individual_subfits

original_grid <- expand.grid(rep = seq_len(50L), nsplits = seq(from = 2L, to = 20L, by = 1L))
target_idx <- which(original_grid$rep == now_rep & original_grid$nsplits == now_nsplit)
if (length(target_idx) != 1L) stop("Could not map task to original MCMC object.", call. = FALSE)
if (individual_nsplits[target_idx] != now_nsplit) {
  stop("Stored nsplit does not match requested nsplit.", call. = FALSE)
}

target_subposteriors <- individual_subfits[[target_idx]]
maxiter <- if (rev_scale == "heavy") rev_smoke_int(100L, 5L) else rev_smoke_int(50L, 5L)
abstol <- 1e-8
init_seed <- 600000L + 10000L * now_rep + 100L * now_nsplit + now_supp
set.seed(init_seed)
fit <- rev_elapsed(
  T4transport::rbaryGD(target_subposteriors,
                       num_support = now_supp,
                       alpha = now_alpha,
                       maxiter = maxiter,
                       abstol = abstol)
)
est_bary <- fit$value
hsum <- rev_history_summary(est_bary$history)

metadata <- list(experiment = "W1_wasp_alpha_barycenters",
                 task_id = myid,
                 data_file = data_file,
                 original_index = target_idx,
                 rep = now_rep,
                 nsplit = now_nsplit,
                 support = now_supp,
                 alpha = now_alpha,
                 maxiter = maxiter,
                 abstol = abstol,
  experiment_scale = rev_scale,
                 init_seed = init_seed)
metrics <- data.frame(
  task_id = myid,
  original_index = target_idx,
  rep = now_rep,
  nsplit = now_nsplit,
  support = now_supp,
  alpha = now_alpha,
  runtime_sec = fit$elapsed,
  niter = est_bary$niter,
  history_first = hsum$first,
  history_last = hsum$last,
  history_min = hsum$min,
  history_length = hsum$n,
  monotone_violations = hsum$monotone_violations,
  max_history_increase = hsum$max_increase
)

save(metadata, metrics, est_bary, file = save_file)
