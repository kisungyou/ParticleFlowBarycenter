# Revision Experiment D1: Digit barycentric prototypes with alpha/support/split sensitivity.
# Purpose: quantify support-size and step-size effects with reproducible Otsu
# preprocessing and multiple train/test splits.
#
# Usage: Rscript 01_digit_centroids_alpha.R <task_id>

script_dir <- getwd()
source(file.path(dirname(script_dir), "common", "rev_utils.R"))
rev_require(c("T4transport"))

make_digit_split <- function(labels, rep_id, n_train_per_digit = 100L, n_test_per_digit = 250L) {
  set.seed(800000L + rep_id)
  labs <- sort(unique(labels))
  train_idx <- integer(0)
  test_idx <- integer(0)
  for (lab in labs) {
    pool <- which(labels == lab)
    pool <- sample(pool, length(pool), replace = FALSE)
    ntr <- min(n_train_per_digit, length(pool))
    train_lab <- pool[seq_len(ntr)]
    rem <- setdiff(pool, train_lab)
    nts <- min(n_test_per_digit, length(rem))
    test_lab <- rem[seq_len(nts)]
    train_idx <- c(train_idx, train_lab)
    test_idx <- c(test_idx, test_lab)
  }
  list(train_idx = sample(train_idx, length(train_idx), replace = FALSE),
       test_idx = sample(test_idx, length(test_idx), replace = FALSE))
}

myid <- rev_arg_id(default = 1L)
out_dir <- rev_output_dir(script_dir, "centroids")
rev_ensure_dir(out_dir)

rev_scale <- tolower(Sys.getenv("REV_EXPERIMENT_SCALE", unset = "core"))
if (rev_scale == "heavy") {
  vec_rep <- seq_len(5L)
  vec_supp <- c(20L, 40L, 80L, 160L, 200L)
  vec_alpha <- c(1.00, 0.75, 0.50, 0.25)
} else {
  vec_rep <- seq_len(3L)
  vec_supp <- c(40L, 80L, 160L)
  vec_alpha <- c(1.00, 0.50)
}
vec_digit <- 0:9
config_grid <- expand.grid(rep = vec_rep, digit = vec_digit,
                           support = vec_supp, alpha = vec_alpha)
if (myid > nrow(config_grid)) quit(save = "no", status = 0, runLast = FALSE)
now_rep <- as.integer(config_grid$rep[myid])
now_digit <- as.integer(config_grid$digit[myid])
now_supp <- as.integer(config_grid$support[myid])
now_alpha <- as.numeric(config_grid$alpha[myid])

save_file <- file.path(out_dir, sprintf("result_%05d.RData", myid))
if (file.exists(save_file)) quit(save = "no", status = 0, runLast = FALSE)

data(digits, package = "T4transport")
labels <- as.integer(digits$label)
split <- make_digit_split(labels, now_rep, n_train_per_digit = rev_smoke_int(100L, 10L), n_test_per_digit = rev_smoke_int(250L, 5L))
idx_digit <- split$train_idx[labels[split$train_idx] == now_digit]
sub_data <- digits$image[idx_digit]
sub_coords <- vector("list", length(sub_data))
num_active_pixels <- integer(length(sub_data))
otsu_thresholds <- numeric(length(sub_data))
for (i in seq_along(sub_data)) {
  otsu_thresholds[i] <- rev_otsu(sub_data[[i]], num_bins = 256L)
  sub_coords[[i]] <- rev_img2coords(sub_data[[i]], num_bins = 256L)
  num_active_pixels[i] <- nrow(sub_coords[[i]])
}

maxiter <- if (rev_scale == "heavy") rev_smoke_int(100L, 5L) else rev_smoke_int(50L, 5L)
abstol <- 1e-8
init_seed <- 900000L + 10000L * now_rep + 100L * now_digit + now_supp
set.seed(init_seed)
fit <- rev_elapsed(
  T4transport::rbaryGD(sub_coords,
                       num_support = now_supp,
                       alpha = now_alpha,
                       maxiter = maxiter,
                       abstol = abstol)
)
est_bary <- fit$value
hsum <- rev_history_summary(est_bary$history)

# Frechet variation around the prototype.
est_variation <- mean(vapply(sub_coords, function(X) {
  T4transport::wasserstein(est_bary$support, X, p = 2)$distance^2
}, numeric(1)))

metrics <- data.frame(
  task_id = myid,
  rep = now_rep,
  digit = now_digit,
  support = now_supp,
  alpha = now_alpha,
  n_train = length(sub_coords),
  active_pixels_mean = mean(num_active_pixels),
  active_pixels_sd = stats::sd(num_active_pixels),
  otsu_threshold_mean = mean(otsu_thresholds),
  runtime_sec = fit$elapsed,
  niter = est_bary$niter,
  history_first = hsum$first,
  history_last = hsum$last,
  history_min = hsum$min,
  history_length = hsum$n,
  monotone_violations = hsum$monotone_violations,
  max_history_increase = hsum$max_increase,
  frechet_variation = est_variation
)
metadata <- list(experiment = "D1_digit_centroids_alpha",
                 split = split,
                 maxiter = maxiter,
                 abstol = abstol,
  experiment_scale = rev_scale,
                 init_seed = init_seed,
                 preprocessing = "Otsu thresholding with 256 bins; foreground pixels converted to centered coordinates in [-1,1]^2 with uniform weights")

save(metadata, metrics, est_bary, sub_coords, sub_data, file = save_file)
