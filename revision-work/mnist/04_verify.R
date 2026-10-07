# Independently check prediction metrics with the original R scoring function.
# Locate inputs/helpers relative to this script, from any working directory.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_arg <- gsub("~+~", " ", script_arg, fixed=TRUE) # Decode spaces encoded by Rscript.
stopifnot(length(script_arg) == 1L)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
repo_root <- normalizePath(file.path(script_dir, "../.."))
work_dir <- Sys.getenv("MNIST_WORK_DIR", unset = script_dir)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)
source(file.path(repo_root, "code-real-digits/common/rev_utils.R"))

stopifnot(requireNamespace("jsonlite", quietly = TRUE))
labels <- read.csv("images.csv")$label
splits <- jsonlite::fromJSON("splits.json", simplifyVector = FALSE)
results <- lapply(1:3, function(rep_id) {
  pred <- read.csv(sprintf("rep%d/predictions.csv", rep_id))
  saved <- jsonlite::fromJSON(sprintf("rep%d/metrics.json", rep_id))
  split <- splits[[rep_id]]
  train <- as.integer(unlist(split$train_idx))
  test <- as.integer(unlist(split$test_idx))
  stopifnot(length(intersect(train, test)) == 0L,
            all(table(labels[train]) == 100L), all(table(labels[test]) == 250L),
            identical(pred$image_id, test), identical(pred$actual, labels[test]))
  score <- rev_classification_metrics(pred$predicted, pred$actual)
  max_difference <- max(abs(score - unlist(saved[names(score)])))
  stopifnot(max_difference < 1e-12)
  for (lab in 0:9) {
    record <- jsonlite::fromJSON(sprintf("rep%d/prototype_%d.json", rep_id, lab))
    stopifnot(record$selected_image %in% record$candidates,
              all(record$candidates %in% train),
              all(labels[record$candidates] == lab),
              record$selected_image == record$candidates[which.min(record$objectives)])
  }
  list(rep = rep_id, n_train = length(train), n_test = length(test),
       independent_R_metric_max_difference = max_difference,
       balanced_disjoint_split = TRUE, training_only_candidates = TRUE)
})
jsonlite::write_json(list(passed = TRUE, splits = results), "verification.json",
                     pretty = TRUE, auto_unbox = TRUE)
cat("All three split metrics, split integrity, and candidate selections verified.\n")
