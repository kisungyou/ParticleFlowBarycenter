# final visualization for revised MNIST/digit prototype experiment
# Style follows the existing visualization scripts:
# ggplot2 + patchwork, colorblind-friendly colors, theme_classic(),
# panel labels (A), (B), etc., and figures saved under ./figures.

# init --------------------------------------------------------------------
rm(list=ls())
graphics.off()
pacman::p_load(rstudioapi,
               ggplot2,
               patchwork,
               dplyr,
               readr,
               tidyr,
               MASS)

# Set working directory to the script location when run in RStudio; otherwise
# keep the current working directory, which should be the project root.
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# helper for paths ---------------------------------------------------------
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
  return(NA_character_)
}

run_label <- Sys.getenv("REV_RUN_LABEL", unset = "core")
if (!nzchar(run_label)) run_label <- "core"

# assembled revision metrics
path_centroids <- find_file(
  file.path(getwd(), "real-1-digits-extended", paste0("assembled_digit_centroids-", run_label, ".csv")),
  file.path(getwd(), "real-1-digits-extended", "assembled_digit_centroids-core.csv"),
  file.path(getwd(), "real-1-digits-extended", "assembled_digit_centroids.csv"),
  file.path(getwd(), paste0("assembled_digit_centroids-", run_label, ".csv")),
  file.path(getwd(), "assembled_digit_centroids-core.csv"),
  file.path(getwd(), "assembled_digit_centroids.csv")
)

path_classification <- find_file(
  file.path(getwd(), "real-1-digits-extended", paste0("assembled_digit_classification-", run_label, ".csv")),
  file.path(getwd(), "real-1-digits-extended", "assembled_digit_classification-core.csv"),
  file.path(getwd(), "real-1-digits-extended", "assembled_digit_classification.csv"),
  file.path(getwd(), paste0("assembled_digit_classification-", run_label, ".csv")),
  file.path(getwd(), "assembled_digit_classification-core.csv"),
  file.path(getwd(), "assembled_digit_classification.csv")
)

path_baselines <- find_optional_file(
  file.path(getwd(), "real-1-digits-extended", paste0("assembled_digit_baselines-", run_label, ".csv")),
  file.path(getwd(), "real-1-digits-extended", "assembled_digit_baselines-core.csv"),
  file.path(getwd(), "real-1-digits-extended", "assembled_digit_baselines.csv"),
  file.path(getwd(), paste0("assembled_digit_baselines-", run_label, ".csv")),
  file.path(getwd(), "assembled_digit_baselines-core.csv"),
  file.path(getwd(), "assembled_digit_baselines.csv")
)

# result directory with prototype RData files, used for the qualitative prototype figure
centroid_dir_candidates <- c(
  file.path(getwd(), "real-1-digits-extended", paste0("centroids-", run_label)),
  file.path(getwd(), "real-1-digits-extended", "centroids-core"),
  file.path(getwd(), "real-1-digits-extended", "centroids"),
  file.path(getwd(), paste0("centroids-", run_label)),
  file.path(getwd(), "centroids-core"),
  file.path(getwd(), "centroids")
)
centroid_dir_candidates <- centroid_dir_candidates[dir.exists(centroid_dir_candidates)]
centroid_dir <- if (length(centroid_dir_candidates) > 0L) centroid_dir_candidates[1L] else NA_character_

# load data ----------------------------------------------------------------
df_centroids <- readr::read_csv(path_centroids, show_col_types = FALSE)
df_class <- readr::read_csv(path_classification, show_col_types = FALSE)
df_base <- NULL
if (!is.na(path_baselines)) {
  df_base <- readr::read_csv(path_baselines, show_col_types = FALSE)
}

# output folders
fig_dir <- file.path(getwd(), "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
manuscript_fig_dir <- file.path(dirname(getwd()), "manuscript", "figures")

copy_to_manuscript <- function(path) {
  if (dir.exists(manuscript_fig_dir)) {
    file.copy(path, manuscript_fig_dir, overwrite = TRUE, recursive = TRUE)
  }
}

# plotting parameters -----------------------------------------------------
cb_colors <- c(
  "#0072B2", # blue
  "#E69F00", # orange
  "#009E73", # bluish green
  "#CC79A7", # reddish purple
  "#56B4E9", # sky blue
  "#D55E00"  # vermillion
)

# helper functions --------------------------------------------------------
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

aux_alpha_label <- function(x) {
  xnum <- as.numeric(x)
  factor(paste0("alpha = ", format(xnum, trim = TRUE)),
         levels = c("alpha = 1", "alpha = 0.5", "alpha = 0.75", "alpha = 0.25"))
}

aux_theme <- function(angle_x = 0) {
  theme_classic() +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0),
      legend.position = "bottom",
      axis.text.x = element_text(angle = angle_x, hjust = ifelse(angle_x == 0, 0.5, 1))
    )
}

# Otsu and image helpers, duplicated here so the figure is self-contained ---
aux_otsu <- function(x, num_bins = 256L) {
  x <- as.vector(as.numeric(x))
  rng <- range(x, finite = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) return(rng[1])
  h <- hist(x, breaks = seq(rng[1], rng[2], length.out = num_bins + 1L), plot = FALSE)
  counts <- h$counts
  mids <- h$mids
  probs <- counts / sum(counts)
  best_thresh <- h$breaks[2L]
  max_between_var <- -Inf
  for (i in seq_len(length(probs) - 1L)) {
    w0 <- sum(probs[1:i])
    w1 <- sum(probs[(i + 1L):length(probs)])
    if (w0 == 0 || w1 == 0) next
    mu0 <- sum(mids[1:i] * probs[1:i]) / w0
    mu1 <- sum(mids[(i + 1L):length(probs)] * probs[(i + 1L):length(probs)]) / w1
    between_var <- w0 * w1 * (mu0 - mu1)^2
    if (between_var > max_between_var) {
      max_between_var <- between_var
      best_thresh <- h$breaks[i + 1L]
    }
  }
  best_thresh
}

aux_image_df <- function(mat) {
  mat <- as.matrix(mat)
  nr <- nrow(mat); nc <- ncol(mat)
  out <- expand.grid(row = seq_len(nr), col = seq_len(nc))
  out$value <- as.vector(mat)
  out$x <- out$col
  out$y <- nr - out$row + 1L
  out
}

aux_bin_coords_for_plot <- function(mat, threshold) {
  mat <- as.matrix(mat)
  idx <- which(mat > threshold, arr.ind = TRUE)
  if (nrow(idx) == 0L) return(data.frame(x = 0, y = 0))
  nr <- nrow(mat); nc <- ncol(mat)
  data.frame(
    x = 2 * (idx[, "col"] - 1) / (nc - 1) - 1,
    # negative sign makes the display upright relative to row coordinates
    y = -(2 * (idx[, "row"] - 1) / (nr - 1) - 1)
  )
}

aux_support_for_plot <- function(X) {
  X <- as.matrix(X)
  data.frame(x = X[, 1], y = -X[, 2])
}

aux_make_kde <- function(df, support_label, lim = c(-1.1, 1.1)) {
  if (nrow(df) < 3L || length(unique(df$x)) < 2L || length(unique(df$y)) < 2L) {
    return(data.frame())
  }
  kd <- MASS::kde2d(df$x, df$y, n = 120, lims = c(lim, lim))
  out <- expand.grid(x = kd$x, y = kd$y)
  out$z <- as.vector(kd$z)
  out$support_lab <- support_label
  out
}

# factor labels ------------------------------------------------------------
df_centroids <- df_centroids %>%
  dplyr::mutate(
    rep = as.integer(rep),
    digit = as.integer(digit),
    support = as.integer(support),
    alpha = as.numeric(alpha),
    support_lab = factor(support, levels = sort(unique(support))),
    alpha_lab = aux_alpha_label(alpha)
  )

df_class <- df_class %>%
  dplyr::mutate(
    rep = as.integer(rep),
    support = as.integer(support),
    alpha = as.numeric(alpha),
    support_lab = factor(support, levels = sort(unique(support))),
    alpha_lab = aux_alpha_label(alpha)
  )

if (!is.null(df_base)) {
  df_base <- df_base %>%
    dplyr::mutate(
      rep = as.integer(rep),
      method = as.character(method),
      method_lab = dplyr::recode(method,
                                 euclidean_class_mean = "Euclidean class mean",
                                 euclidean_1nn = "Euclidean 1NN",
                                 .default = method)
    )
}

alpha_cols <- c("alpha = 1" = cb_colors[1],
                "alpha = 0.5" = cb_colors[2],
                "alpha = 0.75" = cb_colors[3],
                "alpha = 0.25" = cb_colors[4])

alpha_ltys <- c("alpha = 1" = "solid",
                "alpha = 0.5" = "dashed",
                "alpha = 0.75" = "dotdash",
                "alpha = 0.25" = "dotted")

# figure 1 : preprocessing example ----------------------------------------
# This reproduces the manuscript's point-cloud construction example:
# original raster, thresholded lattice, intensity histogram with Otsu threshold,
# and foreground point cloud with uniform weights.
save_path <- file.path(fig_dir, "fig-real-digits-1.png")
run_fig1 <- requireNamespace("T4transport", quietly = TRUE)

if (run_fig1) {
  data(digits, package = "T4transport")
  labels <- as.integer(digits$label)
  idx8 <- which(labels == 8L)[1L]
  img <- as.matrix(digits$image[[idx8]])
  thr <- aux_otsu(img, num_bins = 256L)
  img_df <- aux_image_df(img)
  bin_df <- img_df
  bin_df$value <- as.numeric(img > thr)
  pt_df <- aux_bin_coords_for_plot(img, thr)
  hist_df <- data.frame(value = as.vector(img))

  fig1_a <- ggplot(img_df, aes(x = col, y = y, fill = value)) +
    geom_raster() +
    coord_equal(expand = FALSE) +
    scale_fill_gradient(low = "white", high = "black", guide = "none") +
    ggtitle("(A)") +
    theme_void() +
    theme(panel.border = element_rect(color = "gray50", fill = NA, linewidth = 0.5),
          plot.title = element_text(margin = margin(b = 6)))

  fig1_b <- ggplot(bin_df, aes(x = col, y = y, fill = value)) +
    geom_raster() +
    coord_equal(expand = FALSE) +
    scale_fill_gradient(low = "white", high = "black", guide = "none") +
    ggtitle("(B)") +
    theme_void() +
    theme(panel.border = element_rect(color = "gray50", fill = NA, linewidth = 0.5),
          plot.title = element_text(margin = margin(b = 6)))

  fig1_c <- ggplot(hist_df, aes(x = value)) +
    geom_histogram(bins = 30, fill = "gray80", color = "gray40", linewidth = 0.2) +
    geom_vline(xintercept = thr, color = cb_colors[2], linewidth = 0.8, linetype = "dashed") +
    labs(x = "Pixel intensity", y = "Count") +
    ggtitle("(C)") +
    aux_theme()

  fig1_d <- ggplot(pt_df, aes(x = x, y = y)) +
    geom_point(size = 0.9, alpha = 0.85, color = cb_colors[1]) +
    coord_equal(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), expand = FALSE) +
    labs(x = "x", y = "y") +
    ggtitle("(D)") +
    aux_theme()

  fig1_final <- fig1_a + fig1_b + fig1_c + fig1_d +
    plot_layout(nrow = 1)
  plot(fig1_final)
  ggsave(filename = save_path, plot = fig1_final, height = 3.2, width = 9, units = "in")
  copy_to_manuscript(save_path)
} else {
  warning("Skipping fig-real-digits-1.png because T4transport is unavailable.")
}

# figure 2 : digit-8 prototype evolution ----------------------------------
# This uses the saved centroid RData files so that the displayed prototypes
# are exactly those used in the classification study.
save_path <- file.path(fig_dir, "fig-real-digits-2.png")
run_fig2 <- !is.na(centroid_dir)

find_centroid_result <- function(rep_id = 1L, digit_id = 8L, support_size, alpha_value = 1.0) {
  if (is.na(centroid_dir)) return(NULL)
  files <- sort(list.files(centroid_dir, pattern = "^result_.*\\.RData$", full.names = TRUE))
  if (length(files) == 0L) return(NULL)
  for (ff in files) {
    env <- new.env(parent = emptyenv())
    load(ff, envir = env)
    if (exists("metrics", envir = env) && exists("est_bary", envir = env)) {
      met <- get("metrics", envir = env)
      if (isTRUE(as.integer(met$rep[1]) == rep_id) &&
          isTRUE(as.integer(met$digit[1]) == digit_id) &&
          isTRUE(as.integer(met$support[1]) == support_size) &&
          isTRUE(abs(as.numeric(met$alpha[1]) - alpha_value) < 1e-12)) {
        return(get("est_bary", envir = env))
      }
    }
  }
  NULL
}

if (run_fig2) {
  # Representative display: digit 8, first split, full MM step.
  preferred_supports <- sort(unique(df_centroids$support[df_centroids$digit == 8L & df_centroids$alpha == 1]))
  if (length(preferred_supports) > 5L) {
    # Keep the panel readable in heavy mode.
    preferred_supports <- preferred_supports[round(seq(1, length(preferred_supports), length.out = 5L))]
  }
  if (length(preferred_supports) == 0L) preferred_supports <- sort(unique(df_centroids$support))

  support_list <- vector("list", length(preferred_supports))
  names(support_list) <- preferred_supports
  for (ii in seq_along(preferred_supports)) {
    fit <- find_centroid_result(rep_id = 1L, digit_id = 8L,
                                support_size = preferred_supports[ii],
                                alpha_value = 1.0)
    if (!is.null(fit)) {
      tmp <- aux_support_for_plot(fit$support)
      tmp$support <- preferred_supports[ii]
      tmp$support_lab <- paste0("m = ", preferred_supports[ii])
      support_list[[ii]] <- tmp
    }
  }
  df_proto <- dplyr::bind_rows(support_list)

  if (nrow(df_proto) == 0L) {
    warning("Skipping fig-real-digits-2.png because no matching digit-8 centroid RData files were found.")
  } else {
    df_proto$support_lab <- factor(df_proto$support_lab, levels = unique(df_proto$support_lab[order(df_proto$support)]))
    df_kde <- dplyr::bind_rows(lapply(split(df_proto, df_proto$support_lab), function(dd) {
      aux_make_kde(dd, unique(dd$support_lab))
    }))
    if (nrow(df_kde) > 0L) {
      df_kde$support_lab <- factor(df_kde$support_lab, levels = levels(df_proto$support_lab))
    }

    fig2_a <- ggplot(df_proto, aes(x = x, y = y)) +
      geom_point(size = 0.65, alpha = 0.9, color = cb_colors[1]) +
      facet_wrap(~ support_lab, nrow = 1) +
      coord_equal(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), expand = FALSE) +
      labs(x = NULL, y = NULL) +
      ggtitle("(A)") +
      theme_void() +
      theme(strip.text = element_text(size = 10),
            panel.border = element_rect(color = "gray50", fill = NA, linewidth = 0.4),
            plot.title = element_text(margin = margin(b = 6)))

    fig2_b <- ggplot(df_kde, aes(x = x, y = y)) +
      geom_raster(aes(fill = z)) +
      geom_contour(aes(z = z), color = "white", linewidth = 0.25, bins = 6) +
      facet_wrap(~ support_lab, nrow = 1) +
      coord_equal(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), expand = FALSE) +
      scale_fill_viridis_c(guide = "none") +
      labs(x = NULL, y = NULL) +
      ggtitle("(B)") +
      theme_void() +
      theme(strip.text = element_text(size = 10),
            panel.border = element_rect(color = "gray50", fill = NA, linewidth = 0.4),
            plot.title = element_text(margin = margin(b = 6)))

    fig2_final <- fig2_a / fig2_b
    plot(fig2_final)
    ggsave(filename = save_path, plot = fig2_final, height = 4.6, width = 9, units = "in")
    copy_to_manuscript(save_path)
  }
}

# figure 3 : nearest-prototype classification metrics ---------------------
save_path <- file.path(fig_dir, "fig-real-digits-3.png")

make_metric_summary <- function(metric_name) {
  aux_ci(df_class, c("support", "support_lab", "alpha_lab"), metric_name)
}

make_baseline_summary <- function(metric_name) {
  if (is.null(df_base) || !(metric_name %in% names(df_base))) return(data.frame())
  df_base %>%
    dplyr::group_by(method_lab) %>%
    dplyr::summarise(mean = mean(.data[[metric_name]], na.rm = TRUE),
                     sd = stats::sd(.data[[metric_name]], na.rm = TRUE),
                     n = sum(!is.na(.data[[metric_name]])),
                     se = dplyr::if_else(n > 1 & is.finite(sd), sd/sqrt(n), 0),
                     lower = mean - 1.96*se,
                     upper = mean + 1.96*se,
                     .groups = "drop") %>%
    dplyr::filter(n > 0, is.finite(mean))
}

metric_plot <- function(metric_name, panel_label, ylab) {
  dat <- make_metric_summary(metric_name)
  bdat <- make_baseline_summary(metric_name)
  p <- ggplot(dat, aes(x = support_lab, y = mean, color = alpha_lab, group = alpha_lab)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower, ymax = upper), linewidth = 0.45, alpha = 0.7, width = 0) +
    scale_color_manual(values = alpha_cols, name = "Step size", drop = FALSE) +
    labs(x = "Support size", y = ylab, color = NULL) +
    ggtitle(panel_label) +
    aux_theme(angle_x = 45)
  if (nrow(bdat) > 0L) {
    # Overlay simple Euclidean baselines as horizontal reference lines.
    bdat$line_type <- dplyr::recode(bdat$method_lab,
                                    "Euclidean class mean" = "Class mean",
                                    "Euclidean 1NN" = "1NN",
                                    .default = bdat$method_lab)
    p <- p + geom_hline(data = bdat,
                        aes(yintercept = mean, linetype = line_type),
                        color = "black", linewidth = 0.45, inherit.aes = FALSE) +
      scale_linetype_manual(values = c("Class mean" = "dashed", "1NN" = "dotted"),
                            name = "Euclidean baseline")
  }
  p
}

fig3_a <- metric_plot("Accuracy", "(A)", "Accuracy")
fig3_b <- metric_plot("MacroPrecision", "(B)", "Macro precision")
fig3_c <- metric_plot("MacroRecall", "(C)", "Macro recall")
fig3_d <- metric_plot("MacroF1", "(D)", "Macro F1")

fig3_final <- (fig3_a + fig3_b) / (fig3_c + fig3_d) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.key.width = grid::unit(1, "cm")) &
  guides(color = guide_legend(title = "Step size", nrow = 1),
         linetype = guide_legend(title = "Euclidean baseline", nrow = 1))
plot(fig3_final)
ggsave(filename = save_path, plot = fig3_final, height = 6.6, width = 9, units = "in")
copy_to_manuscript(save_path)

# figure 4 : runtime and optimization diagnostics -------------------------
save_path <- file.path(fig_dir, "fig-real-digits-4.png")

centroid_rt <- aux_ci(df_centroids, c("support", "support_lab", "alpha_lab"), "runtime_sec")
centroid_iter <- aux_ci(df_centroids, c("support", "support_lab", "alpha_lab"), "niter")
centroid_var <- aux_ci(df_centroids, c("support", "support_lab", "alpha_lab"), "frechet_variation")
class_rt <- aux_ci(df_class, c("support", "support_lab", "alpha_lab"), "runtime_sec")

fig4_a <- ggplot(centroid_rt, aes(x = support_lab, y = mean, color = alpha_lab, group = alpha_lab)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), linewidth = 0.45, alpha = 0.7, width = 0) +
  scale_color_manual(values = alpha_cols, name = "Step size", drop = FALSE) +
  labs(x = "Support size", y = "Runtime (sec)") +
  ggtitle("(A)") + aux_theme(angle_x = 45)

fig4_b <- ggplot(class_rt, aes(x = support_lab, y = mean, color = alpha_lab, group = alpha_lab)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), linewidth = 0.45, alpha = 0.7, width = 0) +
  scale_color_manual(values = alpha_cols, name = "Step size", drop = FALSE) +
  labs(x = "Support size", y = "Runtime (sec)") +
  ggtitle("(B)") + aux_theme(angle_x = 45)

fig4_c <- ggplot(centroid_iter, aes(x = support_lab, y = mean, color = alpha_lab, group = alpha_lab)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), linewidth = 0.45, alpha = 0.7, width = 0) +
  scale_color_manual(values = alpha_cols, name = "Step size", drop = FALSE) +
  labs(x = "Support size", y = "Iterations") +
  ggtitle("(C)") + aux_theme(angle_x = 45)

fig4_d <- ggplot(centroid_var, aes(x = support_lab, y = mean, color = alpha_lab, group = alpha_lab)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), linewidth = 0.45, alpha = 0.7, width = 0) +
  scale_color_manual(values = alpha_cols, name = "Step size", drop = FALSE) +
  labs(x = "Support size", y = "Within-class OT variation") +
  ggtitle("(D)") + aux_theme(angle_x = 45)

fig4_final <- (fig4_a + fig4_b) / (fig4_c + fig4_d) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.key.width = grid::unit(1, "cm")) &
  guides(color = guide_legend(title = "Step size", nrow = 1))
plot(fig4_final)
ggsave(filename = save_path, plot = fig4_final, height = 6.6, width = 9, units = "in")
copy_to_manuscript(save_path)

# summary tables -----------------------------------------------------------
summary_class <- df_class %>%
  dplyr::group_by(support, alpha) %>%
  dplyr::summarise(
    Accuracy_mean = mean(Accuracy, na.rm = TRUE),
    MacroPrecision_mean = mean(MacroPrecision, na.rm = TRUE),
    MacroRecall_mean = mean(MacroRecall, na.rm = TRUE),
    MacroF1_mean = mean(MacroF1, na.rm = TRUE),
    runtime_sec_mean = mean(runtime_sec, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(summary_class, file.path(fig_dir, "table-real-digits-classification-summary.csv"))

if (!is.null(df_base)) {
  summary_base <- df_base %>%
    dplyr::group_by(method_lab) %>%
    dplyr::summarise(
      Accuracy_mean = mean(Accuracy, na.rm = TRUE),
      MacroPrecision_mean = mean(MacroPrecision, na.rm = TRUE),
      MacroRecall_mean = mean(MacroRecall, na.rm = TRUE),
      MacroF1_mean = mean(MacroF1, na.rm = TRUE),
      runtime_sec_mean = mean(runtime_sec, na.rm = TRUE),
      .groups = "drop"
    )
  readr::write_csv(summary_base, file.path(fig_dir, "table-real-digits-baseline-summary.csv"))
}

summary_centroid <- df_centroids %>%
  dplyr::group_by(support, alpha) %>%
  dplyr::summarise(
    runtime_sec_mean = mean(runtime_sec, na.rm = TRUE),
    niter_mean = mean(niter, na.rm = TRUE),
    frechet_variation_mean = mean(frechet_variation, na.rm = TRUE),
    active_pixels_mean = mean(active_pixels_mean, na.rm = TRUE),
    otsu_threshold_mean = mean(otsu_threshold_mean, na.rm = TRUE),
    monotone_violations_total = sum(monotone_violations, na.rm = TRUE),
    max_history_increase = max(max_history_increase, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(summary_centroid, file.path(fig_dir, "table-real-digits-centroid-summary.csv"))

# Preprocessing table from the raw digit data, if available.
if (run_fig1) {
  all_active <- integer(length(digits$image))
  all_total_mass <- numeric(length(digits$image))
  all_thr <- numeric(length(digits$image))
  for (ii in seq_along(digits$image)) {
    mat <- as.matrix(digits$image[[ii]])
    th <- aux_otsu(mat, num_bins = 256L)
    all_thr[ii] <- th
    all_active[ii] <- sum(mat > th)
    all_total_mass[ii] <- sum(mat)
  }
  preprocessing_summary <- data.frame(
    n_images = length(digits$image),
    raw_positive_pixels_mean = mean(vapply(digits$image, function(z) sum(as.matrix(z) > 0), numeric(1))),
    otsu_active_pixels_mean = mean(all_active),
    otsu_active_pixels_sd = stats::sd(all_active),
    otsu_threshold_mean = mean(all_thr),
    otsu_threshold_sd = stats::sd(all_thr),
    raw_total_intensity_mean = mean(all_total_mass)
  )
  readr::write_csv(preprocessing_summary, file.path(fig_dir, "table-real-digits-preprocessing-summary.csv"))
}

message("Saved digit revision figures and tables under: ", fig_dir)
