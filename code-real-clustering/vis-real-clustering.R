# final visualization for revised distributed vector quantization experiment
# This script follows the style of the original vis-real-clustering.R:
# ggplot2 + patchwork, colorblind-friendly colors/linetypes, theme_bw(),
# panel labels, and figures saved under ./figures.

# init --------------------------------------------------------------------
rm(list=ls())
graphics.off()
pacman::p_load(rstudioapi,
               ggplot2,
               patchwork,
               dplyr,
               readr,
               tidyr)

# Set working directory to the script location when run in RStudio; otherwise
# keep the current working directory, which should be the project root.
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  active_path <- rstudioapi::getActiveDocumentContext()$path
  if (nzchar(active_path)) setwd(dirname(active_path))
}

# helpers -----------------------------------------------------------------
find_file <- function(...) {
  candidates <- unlist(list(...))
  candidates <- candidates[nzchar(candidates)]
  for (p in candidates) {
    if (file.exists(p)) return(p)
  }
  stop("None of the candidate files exists: ", paste(candidates, collapse = ", "), call. = FALSE)
}

find_optional_file <- function(...) {
  candidates <- unlist(list(...))
  candidates <- candidates[nzchar(candidates)]
  for (p in candidates) {
    if (file.exists(p)) return(p)
  }
  NA_character_
}

aux_ci <- function(data, group_vars, value_var) {
  data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      mean = mean(.data[[value_var]], na.rm = TRUE),
      sd = stats::sd(.data[[value_var]], na.rm = TRUE),
      n = sum(!is.na(.data[[value_var]])),
      se = dplyr::if_else(n > 1 & is.finite(sd), sd/sqrt(n), 0),
      lower = mean - 1.96*se,
      upper = mean + 1.96*se,
      .groups = "drop"
    ) %>%
    dplyr::filter(n > 0, is.finite(mean))
}

aux_theme <- function(angle_x = 0) {
  theme_bw(base_size = 12) +
    theme(
      aspect.ratio = 1,
      legend.position = "bottom",
      axis.text.x = element_text(angle = angle_x, hjust = ifelse(angle_x == 0, 0.5, 1)),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(margin = margin(b = 7))
    )
}

copy_to_manuscript <- function(path) {
  manuscript_fig_dir <- file.path(dirname(getwd()), "manuscript", "figures")
  if (dir.exists(manuscript_fig_dir)) {
    file.copy(path, manuscript_fig_dir, overwrite = TRUE, recursive = TRUE)
  }
}

# paths and data -----------------------------------------------------------
run_label <- Sys.getenv("REV_RUN_LABEL", unset = "core")
if (!nzchar(run_label)) run_label <- "core"

path_cluster <- find_file(
  file.path(getwd(), "real-2-clustering-extended", paste0("assembled_cluster_extended-", run_label, ".csv")),
  file.path(getwd(), "real-2-clustering-extended", "assembled_cluster_extended-core.csv"),
  file.path(getwd(), "real-2-clustering-extended", "assembled_cluster_extended.csv"),
  file.path(getwd(), paste0("assembled_cluster_extended-", run_label, ".csv")),
  file.path(getwd(), "assembled_cluster_extended-core.csv"),
  file.path(getwd(), "assembled_cluster_extended.csv")
)

df_cluster <- readr::read_csv(path_cluster, show_col_types = FALSE)

fig_dir <- file.path(getwd(), "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# plotting parameters -----------------------------------------------------
cb_colors <- c(
  "#0072B2", # blue
  "#E69F00", # orange
  "#009E73", # bluish green
  "#CC79A7", # reddish purple
  "#56B4E9", # sky blue
  "#D55E00", # vermillion
  "#000000", # black
  "#999999"  # gray
)

# normalize column names and labels ---------------------------------------
if ("Sil" %in% names(df_cluster) && !("Silhouette" %in% names(df_cluster))) {
  df_cluster$Silhouette <- df_cluster$Sil
}
if ("Calinski_Harabasz" %in% names(df_cluster) && !("CH" %in% names(df_cluster))) {
  df_cluster$CH <- df_cluster$Calinski_Harabasz
}

# friendly dataset and algorithm labels
parse_cluster_datasets <- function() {
  raw <- Sys.getenv("CLUSTER_DATASETS", unset = "pbmc,news,fashion")
  out <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  out <- out[nzchar(out)]
  allowed <- c("pbmc", "news", "fashion")
  out <- unique(out[out %in% allowed])
  if (length(out) == 0L) out <- allowed
  out
}
selected_datasets <- parse_cluster_datasets()

pretty_data <- function(x) {
  out <- as.character(x)
  out[out == "pbmc"] <- "PBMC"
  out[out == "news"] <- "NEWS"
  out[out %in% c("fashion", "fmnist", "fashion-mnist")] <- "FASHION"
  out
}

pretty_algorithm <- function(base_algorithm, split, alpha, algorithm) {
  out <- as.character(algorithm)
  is_km <- base_algorithm == "kmeans" | algorithm == "kmeans"
  is_skm <- base_algorithm == "skmeans" | algorithm == "skmeans"
  is_dvq <- base_algorithm == "DVQ" | grepl("^DVQ", algorithm)
  out[is_km] <- "k-means"
  out[is_skm] <- "spherical k-means"
  if (any(is_dvq)) {
    out[is_dvq] <- sprintf("DVQ-%s, alpha = %s",
                           as.character(split[is_dvq]),
                           format(as.numeric(alpha[is_dvq]), trim = TRUE, nsmall = 1))
  }
  out
}

if (!("base_algorithm" %in% names(df_cluster))) df_cluster$base_algorithm <- NA_character_
if (!("split" %in% names(df_cluster))) df_cluster$split <- NA_integer_
if (!("alpha" %in% names(df_cluster))) df_cluster$alpha <- NA_real_

df_cluster <- df_cluster %>%
  dplyr::mutate(
    data = as.character(data),
    data_lab = pretty_data(data),
    algorithm_lab = pretty_algorithm(base_algorithm, split, alpha, algorithm),
    nclust = as.integer(nclust),
    repeat_id = if ("repeat_id" %in% names(.)) as.integer(repeat_id) else NA_integer_
  )

alg_levels <- c("k-means", "spherical k-means",
                "DVQ-5, alpha = 1.0", "DVQ-5, alpha = 0.5",
                "DVQ-10, alpha = 1.0", "DVQ-10, alpha = 0.5",
                "DVQ-2, alpha = 1.0", "DVQ-2, alpha = 0.5")
alg_levels <- alg_levels[alg_levels %in% unique(df_cluster$algorithm_lab)]
if (length(alg_levels) == 0L) alg_levels <- sort(unique(df_cluster$algorithm_lab))

df_cluster <- df_cluster %>%
  dplyr::mutate(
    data_lab = factor(data_lab, levels = pretty_data(selected_datasets)),
    algorithm_lab = factor(algorithm_lab, levels = alg_levels)
  )

alg_cols <- setNames(cb_colors[seq_along(levels(df_cluster$algorithm_lab))],
                     levels(df_cluster$algorithm_lab))
alg_ltys <- setNames(rep(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),
                         length.out = length(levels(df_cluster$algorithm_lab))),
                     levels(df_cluster$algorithm_lab))

# figure 1 : embeddings ---------------------------------------------------
# Prefer the original UMAP objects when available. If they are unavailable,
# fall back to a PCA visualization from the processed data files.
load_old_embedding <- function(name) {
  obj_name <- switch(name,
                     pbmc = "pbmc",
                     news = "news",
                     fashion = "fashion")
  candidates <- c(
    file.path(getwd(), "real-2-clustering-extended", "outcome", paste0("vis_", obj_name, ".RData")),
    file.path(getwd(), "real-2-clustering", "outcome", paste0("vis_", obj_name, ".RData")),
    file.path(getwd(), paste0("vis_", obj_name, ".RData"))
  )
  ff <- find_optional_file(candidates)
  if (is.na(ff)) return(NULL)
  env <- new.env(parent = emptyenv())
  load(ff, envir = env)
  x_obj <- paste0("X_", obj_name, "_umap")
  y_obj <- paste0("Y_", obj_name)
  if (!exists(x_obj, envir = env) || !exists(y_obj, envir = env)) return(NULL)
  list(X = get(x_obj, envir = env), Y = get(y_obj, envir = env), method = "UMAP")
}

load_processed_embedding <- function(name, max_plot = 5000L) {
  local_processed <- file.path(getwd(), "real-2-clustering-extended", "data-processed")
  data_dir <- Sys.getenv("CLUSTER_DATA_DIR", unset = local_processed)
  ff <- switch(name,
               pbmc = file.path(data_dir, "processed_pbmc.RData"),
               news = file.path(data_dir, "processed_news.RData"),
               fashion = file.path(data_dir, "processed_fmnist.RData"))
  if (!file.exists(ff)) return(NULL)
  env <- new.env(parent = emptyenv())
  load(ff, envir = env)
  if (!exists("X_reduced", envir = env) || !exists("Y", envir = env)) return(NULL)
  X <- as.matrix(get("X_reduced", envir = env))
  Y <- as.integer(get("Y", envir = env))
  set.seed(100 + match(name, c("pbmc", "news", "fashion")))
  idx <- seq_len(nrow(X))
  if (length(idx) > max_plot) idx <- sort(sample(idx, max_plot, replace = FALSE))
  Xp <- X[idx, , drop = FALSE]
  Yp <- Y[idx]
  if (ncol(Xp) == 1L) {
    emb <- cbind(Xp[,1], rep(0, nrow(Xp)))
    method <- "feature"
  } else if (ncol(Xp) == 2L) {
    emb <- Xp[,1:2, drop = FALSE]
    method <- "features"
  } else {
    pc <- stats::prcomp(Xp, center = TRUE, scale. = FALSE, rank. = 2)
    emb <- pc$x[,1:2, drop = FALSE]
    method <- "PCA"
  }
  list(X = emb, Y = Yp, method = method)
}

make_embedding_panel <- function(name, title) {
  emb <- load_old_embedding(name)
  if (is.null(emb)) emb <- load_processed_embedding(name)
  if (is.null(emb)) return(NULL)
  df <- data.frame(x = emb$X[,1], y = emb$X[,2], label = as.factor(emb$Y))
  ggplot(df, aes(x = x, y = y, color = label)) +
    geom_point(size = 0.45, alpha = 0.75) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(margin = margin(b = 7))) +
    coord_fixed(ratio = 1) +
    labs(x = paste0(emb$method, "1"), y = paste0(emb$method, "2")) +
    ggtitle(title)
}

panel_titles <- paste0("(", LETTERS[seq_along(selected_datasets)], ")")
panels_embed <- Map(make_embedding_panel, selected_datasets, panel_titles)
panels_embed <- panels_embed[!vapply(panels_embed, is.null, logical(1))]
if (length(panels_embed) > 0L) {
  fig1_final <- wrap_plots(panels_embed, nrow = 1)
  save_path <- file.path(fig_dir, "fig-real-clustering-1.png")
  ggsave(filename = save_path, plot = fig1_final, height = 3.5,
         width = 2.8 * length(panels_embed), unit = "in")
  copy_to_manuscript(save_path)
} else {
  warning("Skipping fig-real-clustering-1.png because embedding/processed data files were not found.")
}

# figure 2 : clustering metrics -------------------------------------------
metric_specs <- data.frame(
  metric = c("ARI", "NMI", "Silhouette", "CH"),
  label  = c("ARI", "NMI", "Silhouette", "CH"),
  stringsAsFactors = FALSE
)
metric_specs <- metric_specs[metric_specs$metric %in% names(df_cluster), , drop = FALSE]
if (nrow(metric_specs) < 4L) warning("Some clustering metric columns were not found; plotting available metrics only.")

metric_summary_list <- lapply(metric_specs$metric, function(mm) {
  aux_ci(df_cluster, c("data_lab", "nclust", "algorithm_lab"), mm) %>%
    dplyr::mutate(metric = mm)
})
metric_summary <- dplyr::bind_rows(metric_summary_list)

make_metric_panel <- function(data_label, metric_name, metric_label, row_title = NULL) {
  dat <- metric_summary %>%
    dplyr::filter(data_lab == data_label, metric == metric_name)
  p <- ggplot(dat, aes(x = nclust, y = mean,
                       color = algorithm_lab, linetype = algorithm_lab,
                       group = algorithm_lab)) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12,
                  linewidth = 0.35, alpha = 0.55) +
    scale_color_manual(values = alg_cols, name = NULL, drop = FALSE) +
    scale_linetype_manual(values = alg_ltys, name = NULL, drop = FALSE) +
    scale_x_continuous(breaks = sort(unique(dat$nclust))) +
    labs(x = "Number of clusters", y = metric_label) +
    aux_theme()
  if (!is.null(row_title)) {
    p <- p + ggtitle(row_title)
  } else {
    p <- p + ggtitle(" ") + theme(plot.title = element_text(color = "transparent"))
  }
  p
}

row_plot <- function(data_label, panel_title) {
  ps <- lapply(seq_len(nrow(metric_specs)), function(i) {
    make_metric_panel(data_label,
                      metric_specs$metric[i],
                      metric_specs$label[i],
                      if (i == 1L) panel_title else NULL)
  })
  wrap_plots(ps, nrow = 1)
}

available_data <- levels(droplevels(df_cluster$data_lab))
row_titles <- setNames(paste0("(", LETTERS[seq_along(available_data)], ") ", available_data), available_data)
rows <- lapply(available_data, function(dd) row_plot(dd, row_titles[[as.character(dd)]]))
fig2_final <- wrap_plots(rows, ncol = 1, guides = "collect") &
  theme(legend.position = "bottom",
        legend.key.width = grid::unit(1.35, "cm"),
        aspect.ratio = 1)

save_path <- file.path(fig_dir, "fig-real-clustering-2.png")
ggsave(filename = save_path, plot = fig2_final, height = 2.8 * length(rows), width = 9.2, unit = "in")
copy_to_manuscript(save_path)

# figure 3 : runtime and optimization diagnostics -------------------------
# Summarize elapsed time and MM diagnostics. Non-DVQ baselines have no
# barycenter history, so the iteration/monotonicity panels focus on DVQ.
central_index <- df_cluster %>%
  dplyr::group_by(data_lab) %>%
  dplyr::summarise(central_nclust = stats::median(sort(unique(nclust))), .groups = "drop")

df_central <- df_cluster %>%
  dplyr::inner_join(central_index, by = "data_lab") %>%
  dplyr::filter(nclust == central_nclust)

rt_central <- aux_ci(df_central, c("data_lab", "algorithm_lab"), "runtime_sec")
rt_curve <- aux_ci(df_cluster, c("data_lab", "nclust", "algorithm_lab"), "runtime_sec")

df_dvq <- df_cluster %>%
  dplyr::filter(grepl("^DVQ", as.character(algorithm_lab)))

niter_central <- aux_ci(df_central %>% dplyr::filter(grepl("^DVQ", as.character(algorithm_lab))),
                        c("data_lab", "algorithm_lab"), "niter")
viol_central <- aux_ci(df_central %>% dplyr::filter(grepl("^DVQ", as.character(algorithm_lab))),
                       c("data_lab", "algorithm_lab"), "monotone_violations")
last_central <- if ("history_last" %in% names(df_cluster)) {
  aux_ci(df_central %>% dplyr::filter(grepl("^DVQ", as.character(algorithm_lab))),
         c("data_lab", "algorithm_lab"), "history_last")
} else {
  data.frame()
}

fig3_a <- ggplot(rt_central, aes(x = algorithm_lab, y = mean, fill = algorithm_lab)) +
  geom_col(width = 0.75) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.35) +
  facet_wrap(~data_lab, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = alg_cols, guide = "none", drop = FALSE) +
  labs(x = NULL, y = "Runtime (sec)") +
  ggtitle("(A)") +
  aux_theme(angle_x = 45) +
  theme(aspect.ratio = 0.8)

fig3_b <- ggplot(rt_curve, aes(x = nclust, y = mean,
                               color = algorithm_lab, linetype = algorithm_lab,
                               group = algorithm_lab)) +
  geom_line(linewidth = 0.65) + geom_point(size = 1.4) +
  facet_wrap(~data_lab, scales = "free_x", nrow = 1) +
  scale_color_manual(values = alg_cols, name = NULL, drop = FALSE) +
  scale_linetype_manual(values = alg_ltys, name = NULL, drop = FALSE) +
  labs(x = "Number of clusters", y = "Runtime (sec)") +
  ggtitle("(B)") +
  aux_theme()

fig3_c <- ggplot(niter_central, aes(x = algorithm_lab, y = mean, fill = algorithm_lab)) +
  geom_col(width = 0.75) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.35) +
  facet_wrap(~data_lab, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = alg_cols, guide = "none", drop = FALSE) +
  labs(x = NULL, y = "Iterations") +
  ggtitle("(C)") +
  aux_theme(angle_x = 45) +
  theme(aspect.ratio = 0.8)

fig3_d <- ggplot(viol_central, aes(x = algorithm_lab, y = mean, fill = algorithm_lab)) +
  geom_col(width = 0.75) +
  facet_wrap(~data_lab, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = alg_cols, guide = "none", drop = FALSE) +
  labs(x = NULL, y = "Monotonicity violations") +
  ggtitle("(D)") +
  aux_theme(angle_x = 45) +
  theme(aspect.ratio = 0.8)

fig3_final <- (fig3_a / fig3_b / fig3_c / fig3_d) +
  plot_layout(guides = "collect", heights = c(1, 1, 1, 1)) &
  theme(legend.position = "bottom",
        legend.key.width = grid::unit(1.35, "cm"))

save_path <- file.path(fig_dir, "fig-real-clustering-3.png")
ggsave(filename = save_path, plot = fig3_final, height = 11.5, width = 9.5, unit = "in")
copy_to_manuscript(save_path)

# summary tables -----------------------------------------------------------
# Central cluster count summary: useful for a compact manuscript table.
central_summary <- df_central %>%
  dplyr::group_by(data_lab, nclust, algorithm_lab) %>%
  dplyr::summarise(
    ARI_mean = mean(ARI, na.rm = TRUE),
    NMI_mean = mean(NMI, na.rm = TRUE),
    Silhouette_mean = mean(Silhouette, na.rm = TRUE),
    CH_mean = mean(CH, na.rm = TRUE),
    runtime_mean = mean(runtime_sec, na.rm = TRUE),
    niter_mean = mean(niter, na.rm = TRUE),
    monotone_violations_total = sum(monotone_violations, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(central_summary, file.path(fig_dir, "table-real-clustering-central-summary.csv"))

# Best ARI by dataset and algorithm over the cluster-count sweep.
best_ari <- df_cluster %>%
  dplyr::group_by(data_lab, algorithm_lab, nclust) %>%
  dplyr::summarise(
    ARI_mean = mean(ARI, na.rm = TRUE),
    NMI_mean = mean(NMI, na.rm = TRUE),
    Silhouette_mean = mean(Silhouette, na.rm = TRUE),
    CH_mean = mean(CH, na.rm = TRUE),
    runtime_mean = mean(runtime_sec, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(data_lab, algorithm_lab) %>%
  dplyr::slice_max(order_by = ARI_mean, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()
readr::write_csv(best_ari, file.path(fig_dir, "table-real-clustering-best-ari-summary.csv"))

# DVQ-only optimization diagnostics across all cluster counts.
dvq_diag <- df_dvq %>%
  dplyr::group_by(data_lab, algorithm_lab) %>%
  dplyr::summarise(
    runtime_mean = mean(runtime_sec, na.rm = TRUE),
    niter_mean = mean(niter, na.rm = TRUE),
    monotone_violations_total = sum(monotone_violations, na.rm = TRUE),
    max_history_increase = max(max_history_increase, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(dvq_diag, file.path(fig_dir, "table-real-clustering-dvq-diagnostics.csv"))

message("Saved clustering revision figures and tables to: ", fig_dir)
