# Revision Experiment C1: clustering baselines, DVQ alpha sensitivity, runtime.
# Purpose: rerun the clustering experiment with repeated random partitions/seeds,
# explicit alpha settings for rbaryGD, and saved runtime diagnostics.
#
# Requires local clustering data directories under this folder:
#   real-2-clustering-extended/data-processed/processed_pbmc.RData
#   real-2-clustering-extended/data-processed/processed_news.RData
#   real-2-clustering-extended/data-processed/processed_fmnist.RData
# The data-raw directory may also be placed here for provenance/reproducibility,
# although this script only needs data-processed.
#
# Usage: Rscript 01_clustering_alpha_runtime.R <task_id>

script_dir <- getwd()
source(file.path(dirname(script_dir), "common", "rev_utils.R"))
rev_require(c("T4transport", "skmeans"))

myid <- rev_arg_id(default = 1L)
out_dir <- rev_output_dir(script_dir, "outcome")
rev_ensure_dir(out_dir)
rev_scale <- tolower(Sys.getenv("REV_EXPERIMENT_SCALE", unset = "core"))

if (rev_scale == "heavy") {
  method_grid <- rbind(
  data.frame(algorithm = "kmeans", split = NA_integer_, alpha = NA_real_),
  data.frame(algorithm = "skmeans", split = NA_integer_, alpha = NA_real_),
  expand.grid(algorithm = "DVQ", split = c(2L, 5L, 10L), alpha = c(1.00, 0.50))
  )
  vec_cluster_index <- seq_len(11L)
  vec_repeat_id <- seq_len(10L)
} else {
  method_grid <- rbind(
    data.frame(algorithm = "kmeans", split = NA_integer_, alpha = NA_real_),
    data.frame(algorithm = "skmeans", split = NA_integer_, alpha = NA_real_),
    expand.grid(algorithm = "DVQ", split = c(5L, 10L), alpha = c(1.00, 0.50))
  )
  vec_cluster_index <- c(4L, 6L, 8L)
  vec_repeat_id <- seq_len(3L)
}
parse_cluster_datasets <- function() {
  raw <- Sys.getenv("CLUSTER_DATASETS", unset = "pbmc,news,fashion")
  out <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  out <- out[nzchar(out)]
  allowed <- c("pbmc", "news", "fashion")
  bad <- setdiff(out, allowed)
  if (length(bad) > 0L) {
    stop("Unknown CLUSTER_DATASETS value(s): ", paste(bad, collapse = ", "),
         ". Allowed values are pbmc, news, fashion.", call. = FALSE)
  }
  unique(out)
}
vec_data <- parse_cluster_datasets()
config_grid <- expand.grid(method_id = seq_len(nrow(method_grid)),
                           data = vec_data,
                           cluster_index = vec_cluster_index,
                           repeat_id = vec_repeat_id)
if (myid > nrow(config_grid)) quit(save = "no", status = 0, runLast = FALSE)
now_cfg <- config_grid[myid, ]
now_method <- method_grid[as.integer(now_cfg$method_id), ]
now_data <- as.character(now_cfg$data)
now_cluster_index <- as.integer(now_cfg$cluster_index)
now_repeat_id <- as.integer(now_cfg$repeat_id)
now_alg <- as.character(now_method$algorithm)
now_split <- as.integer(now_method$split)
now_alpha <- as.numeric(now_method$alpha)

save_file <- file.path(out_dir, sprintf("result_%05d.RData", myid))
if (file.exists(save_file)) quit(save = "no", status = 0, runLast = FALSE)

load_dataset <- function(name) {
  local_processed <- file.path(script_dir, "data-processed")
  local_raw <- file.path(script_dir, "data-raw")
  data_dir <- Sys.getenv("CLUSTER_DATA_DIR", unset = local_processed)
  raw_dir <- Sys.getenv("CLUSTER_RAW_DIR", unset = local_raw)
  if (!dir.exists(data_dir)) {
    stop("Missing clustering processed-data directory: ", data_dir,
         ". Put data-processed under real-2-clustering-extended or export CLUSTER_DATA_DIR.",
         call. = FALSE)
  }
  file <- switch(name,
                 pbmc = file.path(data_dir, "processed_pbmc.RData"),
                 news = file.path(data_dir, "processed_news.RData"),
                 fashion = file.path(data_dir, "processed_fmnist.RData"))
  if (!file.exists(file)) stop("Missing processed data file: ", file,
                               ". Expected the processed file for dataset '", name, "' inside ",
                               data_dir, ".",
                               call. = FALSE)
  load(file)
  range_nclust <- switch(name,
                         pbmc = as.integer(seq(5, 15, by = 1)),
                         news = as.integer(seq(15, 25, by = 1)),
                         fashion = as.integer(seq(5, 15, by = 1)))
  list(X = as.matrix(X_reduced), Y = as.integer(Y), range_nclust = range_nclust, file = file)
}

aux_kmeans <- function(X, k, seed, max_n = 5000L) {
  centers <- rev_kmeans_centers(X, k, seed = seed, nstart = 10L, max_n = max_n)
  rev_assign_to_centers(X, centers)
}

aux_kmeans_centroids <- function(X, k, seed, max_n = 5000L) {
  rev_kmeans_centers(X, k, seed = seed, nstart = 10L, max_n = max_n)
}

aux_skmeans <- function(X, k, seed) {
  set.seed(seed)
  XX <- scale(X)
  norms <- sqrt(rowSums(XX^2))
  norms[norms == 0] <- 1
  XX <- XX / norms
  fit <- skmeans::skmeans(XX, k = k)
  as.integer(fit$cluster)
}

aux_DVQ <- function(X, k, num_splits, alpha, seed) {
  set.seed(seed)
  idx_perturb <- sample(rep(seq_len(num_splits), length.out = nrow(X)))
  subset_centers <- vector("list", length = num_splits)
  for (i in seq_len(num_splits)) {
    Xi <- X[idx_perturb == i, , drop = FALSE]
    subset_centers[[i]] <- aux_kmeans_centroids(Xi, k, seed = seed + 100L * i)
  }
  set.seed(seed + 999L)
  fit <- T4transport::rbaryGD(subset_centers,
                              num_support = k,
                              alpha = alpha,
                              maxiter = if (rev_scale == "heavy") 100L else 50L,
                              abstol = 1e-8)
  centers <- fit$support
  pred <- rev_assign_to_centers(X, centers)
  list(pred = pred, centers = centers, history = fit$history,
       niter = fit$niter, alpha = fit$alpha)
}

dat <- load_dataset(now_data)
X <- dat$X
Y <- dat$Y
now_nclust <- dat$range_nclust[now_cluster_index]
seed <- 1000000L + 10000L * now_repeat_id + 100L * now_cluster_index + as.integer(now_cfg$method_id)

fit <- rev_elapsed({
  if (now_alg == "kmeans") {
    list(pred = aux_kmeans(X, now_nclust, seed = seed), centers = NULL,
         history = numeric(0), niter = NA_integer_, alpha = NA_real_)
  } else if (now_alg == "skmeans") {
    list(pred = aux_skmeans(X, now_nclust, seed = seed), centers = NULL,
         history = numeric(0), niter = NA_integer_, alpha = NA_real_)
  } else if (now_alg == "DVQ") {
    aux_DVQ(X, now_nclust, now_split, now_alpha, seed = seed)
  } else {
    stop("Unknown algorithm", call. = FALSE)
  }
})

hsum <- rev_history_summary(fit$value$history)
method_label <- if (now_alg == "DVQ") sprintf("DVQ-%d-alpha-%.2f", now_split, now_alpha) else now_alg
metrics <- data.frame(
  task_id = myid,
  data = now_data,
  algorithm = method_label,
  base_algorithm = now_alg,
  split = ifelse(is.na(now_split), NA_integer_, now_split),
  alpha = ifelse(is.na(now_alpha), NA_real_, now_alpha),
  cluster_index = now_cluster_index,
  nclust = now_nclust,
  repeat_id = now_repeat_id,
  runtime_sec = fit$elapsed,
  niter = fit$value$niter,
  history_first = hsum$first,
  history_last = hsum$last,
  history_min = hsum$min,
  history_length = hsum$n,
  monotone_violations = hsum$monotone_violations,
  max_history_increase = hsum$max_increase
)
metadata <- list(experiment = "C1_clustering_alpha_runtime",
                 data_file = dat$file,
                 data_dir = Sys.getenv("CLUSTER_DATA_DIR", unset = file.path(script_dir, "data-processed")),
                 raw_dir = Sys.getenv("CLUSTER_RAW_DIR", unset = file.path(script_dir, "data-raw")),
                 seed = seed,
                 experiment_scale = rev_scale)
pred_label <- fit$value$pred
centers <- fit$value$centers
history <- fit$value$history
save(metadata, metrics, pred_label, centers, history, file = save_file)
