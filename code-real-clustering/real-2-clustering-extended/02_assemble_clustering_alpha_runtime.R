# Assemble C1 clustering results and compute ARI/NMI/Silhouette/CH.
script_dir <- getwd()
source(file.path(dirname(script_dir), "common", "rev_utils.R"))
rev_require(c("cluster", "clusterCrit", "mclustcomp"))

load_dataset <- function(name) {
  local_processed <- file.path(script_dir, "data-processed")
  data_dir <- Sys.getenv("CLUSTER_DATA_DIR", unset = local_processed)
  if (!dir.exists(data_dir)) {
    stop("Missing clustering processed-data directory: ", data_dir,
         ". Put data-processed under real-2-clustering-extended or export CLUSTER_DATA_DIR.",
         call. = FALSE)
  }
  file <- switch(name,
                 pbmc = file.path(data_dir, "processed_pbmc.RData"),
                 news = file.path(data_dir, "processed_news.RData"),
                 fashion = file.path(data_dir, "processed_fmnist.RData"))
  if (!file.exists(file)) {
    stop("Missing processed data file: ", file,
         ". Expected the processed file for dataset '", name, "' inside ",
         data_dir, ".",
         call. = FALSE)
  }
  load(file)
  list(X = as.matrix(X_reduced), Y = as.integer(Y), D = stats::dist(X_reduced), file = file)
}

aux_clust_scores <- function(lab_true, lab_pred, D, X) {
  score_comp <- mclustcomp::mclustcomp(lab_true, lab_pred, types = c("adjrand", "nmi1"))
  score_ARI <- score_comp$scores[1]
  score_NMI <- score_comp$scores[2]
  score_SIL <- mean(cluster::silhouette(lab_pred, D)[, "sil_width"])
  score_CH <- clusterCrit::intCriteria(X, as.integer(lab_pred), "Calinski_Harabasz")$calinski_harabasz
  c(ARI = score_ARI, NMI = score_NMI, Silhouette = score_SIL, CH = score_CH)
}

run_dir <- rev_output_dir(script_dir, "outcome")
files <- sort(list.files(run_dir, pattern = "^result_.*\\.RData$", full.names = TRUE))
if (length(files) == 0L) stop("No result files found in ", run_dir, call. = FALSE)

# First read the lightweight metrics from each file so smoke tests only load the
# datasets that were actually touched. This avoids requiring all large datasets
# for a one-task path check.
metric_cache <- vector("list", length(files))
needed_names <- character(0)
for (i in seq_along(files)) {
  load(files[i])
  metric_cache[[i]] <- list(metrics = metrics, pred_label = pred_label)
  needed_names <- unique(c(needed_names, as.character(metrics$data[1])))
}
datasets <- setNames(lapply(needed_names, load_dataset), needed_names)

rows <- vector("list", length(files))
for (i in seq_along(files)) {
  metrics <- metric_cache[[i]]$metrics
  pred_label <- metric_cache[[i]]$pred_label
  ds <- datasets[[as.character(metrics$data[1])]]
  scores <- aux_clust_scores(ds$Y, pred_label, ds$D, ds$X)
  rows[[i]] <- cbind(metrics, as.data.frame(as.list(scores)))
  print(paste0("assembled ", i, "/", length(files)))
}
df_cluster_extended <- do.call(rbind, rows)
df_cluster_extended <- df_cluster_extended[order(df_cluster_extended$data,
                                                 df_cluster_extended$nclust,
                                                 df_cluster_extended$algorithm,
                                                 df_cluster_extended$repeat_id), ]
save(df_cluster_extended, file = rev_output_file(script_dir, "assembled_cluster_extended"))
write.csv(df_cluster_extended, file = rev_output_file(script_dir, "assembled_cluster_extended", ext = ".csv"), row.names = FALSE)
