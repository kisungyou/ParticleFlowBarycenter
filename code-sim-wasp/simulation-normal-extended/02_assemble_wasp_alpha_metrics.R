# Assemble W1 results and compute posterior aggregation diagnostics.
# Adds mean/covariance error, semi-discrete W2, marginal coverage/width, and
# posterior predictive mean error/credible interval diagnostics.

script_dir <- getwd()
source(file.path(dirname(script_dir), "common", "rev_utils.R"))
rev_require(c("T4transport", "mvtnorm"))

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
  if (length(hit) == 0L) stop("Could not find normal-simplified-reduced.RData. Set WASP_DATA_FILE.", call. = FALSE)
  hit[1L]
}

data_file <- find_data_file()
load(data_file)
analytic_mean <- common_posterior_mean
analytic_cov <- common_posterior_cov
beta_true <- c(1, 1)

set.seed(700000L)
oracle_draws <- mvtnorm::rmvnorm(n = rev_smoke_int(1000L, 100L), mean = analytic_mean, sigma = analytic_cov)
X_test <- matrix(rnorm(rev_smoke_int(2000L, 100L) * 2L), ncol = 2L)
eta_true <- as.vector(X_test %*% beta_true)

summarize_draws <- function(draws) {
  draws <- as.matrix(draws)
  semi_w2 <- T4transport::wasserstein(draws, oracle_draws, p = 2)$distance
  mean_error <- sqrt(sum((colMeans(draws) - analytic_mean)^2))
  cov_error <- rev_bures(stats::cov(draws), analytic_cov)
  q025 <- apply(draws, 2, stats::quantile, probs = 0.025, names = FALSE)
  q975 <- apply(draws, 2, stats::quantile, probs = 0.975, names = FALSE)
  marginal_coverage <- mean((beta_true >= q025) & (beta_true <= q975))
  marginal_width <- mean(q975 - q025)
  eta_draws <- X_test %*% t(draws)  # n_test by n_draws
  eta_mean <- rowMeans(eta_draws)
  eta_q025 <- apply(eta_draws, 1, stats::quantile, probs = 0.025, names = FALSE)
  eta_q975 <- apply(eta_draws, 1, stats::quantile, probs = 0.975, names = FALSE)
  pred_rmse_mean <- sqrt(mean((eta_mean - eta_true)^2))
  pred_interval_coverage <- mean((eta_true >= eta_q025) & (eta_true <= eta_q975))
  pred_interval_width <- mean(eta_q975 - eta_q025)
  c(semi_w2 = semi_w2,
    mean_error = mean_error,
    cov_error = cov_error,
    marginal_coverage = marginal_coverage,
    marginal_width = marginal_width,
    pred_rmse_mean = pred_rmse_mean,
    pred_interval_coverage = pred_interval_coverage,
    pred_interval_width = pred_interval_width)
}

run_dir <- rev_output_dir(script_dir, "runs")
files <- sort(list.files(run_dir, pattern = "^result_.*\\.RData$", full.names = TRUE))
if (length(files) == 0L) stop("No result files found in ", run_dir, call. = FALSE)

rows <- vector("list", length(files))
for (i in seq_along(files)) {
  load(files[i])
  extra <- summarize_draws(est_bary$support)
  rows[[i]] <- cbind(metrics, as.data.frame(as.list(extra)))
  print(paste0("assembled ", i, "/", length(files)))
}
df_wasp_alpha <- do.call(rbind, rows)
df_wasp_alpha <- df_wasp_alpha[order(df_wasp_alpha$rep, df_wasp_alpha$nsplit,
                                     df_wasp_alpha$support, df_wasp_alpha$alpha), ]

# Full-data MCMC reference computed once from stored draws.
full_reference <- as.data.frame(as.list(summarize_draws(common_full_mcmc)))
full_reference$method <- "full_mcmc"

save(df_wasp_alpha, full_reference, analytic_mean, analytic_cov,
     file = rev_output_file(script_dir, "assembled_wasp_alpha_metrics"))
write.csv(df_wasp_alpha, file = rev_output_file(script_dir, "assembled_wasp_alpha_metrics", ext = ".csv"), row.names = FALSE)
write.csv(full_reference, file = rev_output_file(script_dir, "assembled_wasp_full_reference", ext = ".csv"), row.names = FALSE)
