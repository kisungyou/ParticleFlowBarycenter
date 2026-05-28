# Revision Experiment D3: simple Euclidean raster baselines for digit classification.
# Baselines: class mean in normalized pixel space and 1-nearest-neighbor in pixel space.
#
# Usage: Rscript 03_digit_euclidean_baselines.R <task_id>

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
out_dir <- rev_output_dir(script_dir, "baselines")
rev_ensure_dir(out_dir)
rev_scale <- tolower(Sys.getenv("REV_EXPERIMENT_SCALE", unset = "core"))
vec_rep <- if (rev_scale == "heavy") seq_len(5L) else seq_len(3L)
if (myid > length(vec_rep)) quit(save = "no", status = 0, runLast = FALSE)
now_rep <- vec_rep[myid]
save_file <- file.path(out_dir, sprintf("result_%05d.RData", myid))
if (file.exists(save_file)) quit(save = "no", status = 0, runLast = FALSE)

data(digits, package = "T4transport")
labels <- as.integer(digits$label)
split <- make_digit_split(labels, now_rep, n_train_per_digit = rev_smoke_int(100L, 10L), n_test_per_digit = rev_smoke_int(250L, 5L))
train_label <- labels[split$train_idx]
test_label <- labels[split$test_idx]
X_train <- t(vapply(digits$image[split$train_idx], rev_normalized_image_vector, numeric(length(rev_normalized_image_vector(digits$image[[split$train_idx[1]]])))))
X_test <- t(vapply(digits$image[split$test_idx], rev_normalized_image_vector, numeric(ncol(X_train))))
vec_digit <- sort(unique(labels))

# Class mean prototype in normalized raster space.
fit_mean <- rev_elapsed({
  class_means <- matrix(0, nrow = length(vec_digit), ncol = ncol(X_train))
  for (i in seq_along(vec_digit)) {
    class_means[i, ] <- colMeans(X_train[train_label == vec_digit[i], , drop = FALSE])
  }
  dmat <- as.matrix(stats::dist(rbind(X_test, class_means)))[seq_len(nrow(X_test)), nrow(X_test) + seq_len(nrow(class_means)), drop = FALSE]
  pred <- vec_digit[max.col(-dmat, ties.method = "first")]
  pred
})
score_mean <- rev_classification_metrics(fit_mean$value, test_label)

# Euclidean 1NN in normalized raster space.
fit_1nn <- rev_elapsed({
  dmat <- as.matrix(stats::dist(rbind(X_test, X_train)))[seq_len(nrow(X_test)), nrow(X_test) + seq_len(nrow(X_train)), drop = FALSE]
  pred <- train_label[max.col(-dmat, ties.method = "first")]
  pred
})
score_1nn <- rev_classification_metrics(fit_1nn$value, test_label)

metrics <- rbind(
  data.frame(task_id = myid, rep = now_rep, method = "euclidean_class_mean", n_test = length(test_label),
             runtime_sec = fit_mean$elapsed,
             Accuracy = score_mean["Accuracy"], MacroPrecision = score_mean["MacroPrecision"],
             MacroRecall = score_mean["MacroRecall"], MacroF1 = score_mean["MacroF1"]),
  data.frame(task_id = myid, rep = now_rep, method = "euclidean_1nn", n_test = length(test_label),
             runtime_sec = fit_1nn$elapsed,
             Accuracy = score_1nn["Accuracy"], MacroPrecision = score_1nn["MacroPrecision"],
             MacroRecall = score_1nn["MacroRecall"], MacroF1 = score_1nn["MacroF1"])
)
pred_mean <- fit_mean$value
pred_1nn <- fit_1nn$value
confusion_mean <- table(actual = test_label, predicted = pred_mean)
confusion_1nn <- table(actual = test_label, predicted = pred_1nn)
metadata <- list(experiment = "D3_digit_euclidean_baselines", split = split,
                 experiment_scale = rev_scale)
save(metadata, metrics, confusion_mean, confusion_1nn,
     pred_mean, pred_1nn, test_label,
     file = save_file)
