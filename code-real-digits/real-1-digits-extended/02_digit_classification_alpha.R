# Revision Experiment D2: nearest barycentric prototype classification.
# Loads D1 centroids and evaluates balanced test splits. Visualization is deferred.
#
# Usage: Rscript 02_digit_classification_alpha.R <task_id>

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
out_dir <- rev_output_dir(script_dir, "classification")
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
config_grid <- expand.grid(rep = vec_rep, support = vec_supp, alpha = vec_alpha)
if (myid > nrow(config_grid)) quit(save = "no", status = 0, runLast = FALSE)
now_rep <- as.integer(config_grid$rep[myid])
now_supp <- as.integer(config_grid$support[myid])
now_alpha <- as.numeric(config_grid$alpha[myid])

save_file <- file.path(out_dir, sprintf("result_%05d.RData", myid))
if (file.exists(save_file)) quit(save = "no", status = 0, runLast = FALSE)

data(digits, package = "T4transport")
labels <- as.integer(digits$label)
split <- make_digit_split(labels, now_rep, n_train_per_digit = rev_smoke_int(100L, 10L), n_test_per_digit = rev_smoke_int(250L, 5L))
test_label <- labels[split$test_idx]
test_image <- vector("list", length(split$test_idx))
for (i in seq_along(split$test_idx)) {
  test_image[[i]] <- rev_img2coords(digits$image[[split$test_idx[i]]], num_bins = 256L)
}

# Load all 10 class prototypes for this rep/support/alpha.
vec_digit <- 0:9
centroids <- vector("list", length(vec_digit))
centroid_metrics <- vector("list", length(vec_digit))
centroid_grid <- expand.grid(rep = vec_rep, digit = vec_digit,
                             support = vec_supp,
                             alpha = vec_alpha)
for (i in seq_along(vec_digit)) {
  idx <- which(centroid_grid$rep == now_rep & centroid_grid$digit == vec_digit[i] &
                 centroid_grid$support == now_supp & abs(centroid_grid$alpha - now_alpha) < 1e-12)
  if (length(idx) != 1L) stop("Could not find centroid task index.", call. = FALSE)
  fp <- file.path(rev_output_dir(script_dir, "centroids"), sprintf("result_%05d.RData", idx))
  if (!file.exists(fp)) stop("Missing centroid file: ", fp, call. = FALSE)
  load(fp)
  centroids[[i]] <- est_bary$support
  centroid_metrics[[i]] <- metrics
}

fit <- rev_elapsed({
  distmat <- matrix(0, nrow = length(vec_digit), ncol = length(test_image))
  for (i in seq_along(vec_digit)) {
    for (j in seq_along(test_image)) {
      distmat[i, j] <- T4transport::wasserstein(centroids[[i]], test_image[[j]], p = 2)$distance
    }
  }
  pred <- vec_digit[max.col(-t(distmat), ties.method = "first")]
  list(pred = pred, distmat = distmat)
})

scores <- rev_classification_metrics(fit$value$pred, test_label)
metrics <- data.frame(
  task_id = myid,
  rep = now_rep,
  support = now_supp,
  alpha = now_alpha,
  method = "wasserstein_barycentric_prototype",
  n_test = length(test_label),
  runtime_sec = fit$elapsed,
  Accuracy = scores["Accuracy"],
  MacroPrecision = scores["MacroPrecision"],
  MacroRecall = scores["MacroRecall"],
  MacroF1 = scores["MacroF1"]
)
pred_label <- fit$value$pred
distmat <- fit$value$distmat
confusion <- table(actual = test_label, predicted = pred_label)
metadata <- list(experiment = "D2_digit_classification_alpha",
                 split = split,
                 preprocessing = "Otsu thresholding with 256 bins; foreground pixels converted to centered coordinates in [-1,1]^2 with uniform weights",
  experiment_scale = rev_scale)

save(metadata, metrics, confusion, pred_label, test_label, distmat, centroid_metrics,
     file = save_file)
