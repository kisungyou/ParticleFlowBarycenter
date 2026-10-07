# final visualization for revised Bayesian posterior aggregation experiment
# This script follows the style of the original vis-simulation-normal.R:
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
               mvtnorm,
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
path_wasp <- find_file(
  file.path(getwd(), "simulation-normal-extended", paste0("assembled_wasp_alpha_metrics-", run_label, ".csv")),
  file.path(getwd(), "simulation-normal-extended", "assembled_wasp_alpha_metrics-core.csv"),
  file.path(getwd(), "simulation-normal-extended", "assembled_wasp_alpha_metrics.csv"),
  file.path(getwd(), paste0("assembled_wasp_alpha_metrics-", run_label, ".csv")),
  file.path(getwd(), "assembled_wasp_alpha_metrics-core.csv"),
  file.path(getwd(), "assembled_wasp_alpha_metrics.csv")
)

path_full_ref <- find_optional_file(
  file.path(getwd(), "simulation-normal-extended", paste0("assembled_wasp_full_reference-", run_label, ".csv")),
  file.path(getwd(), "simulation-normal-extended", "assembled_wasp_full_reference-core.csv"),
  file.path(getwd(), "simulation-normal-extended", "assembled_wasp_full_reference.csv"),
  file.path(getwd(), paste0("assembled_wasp_full_reference-", run_label, ".csv")),
  file.path(getwd(), "assembled_wasp_full_reference-core.csv"),
  file.path(getwd(), "assembled_wasp_full_reference.csv")
)

# original reduced data, needed only for Fig. 1 density visualization
path_data <- find_optional_file(
  Sys.getenv("WASP_DATA_FILE", unset = ""),
  file.path(getwd(), "simulation-normal-extended", "normal-simplified-reduced.RData"),
  file.path(getwd(), "simulation-normal", "normal-simplified-reduced.RData"),
  file.path(getwd(), "simulation-normal", "assemble-normal-data.RData"),
  file.path(dirname(getwd()), "simulation-normal", "normal-simplified-reduced.RData"),
  file.path(dirname(getwd()), "simulation-normal", "assemble-normal-data.RData")
)

# load data ----------------------------------------------------------------
df_wasp <- readr::read_csv(path_wasp, show_col_types = FALSE)
full_reference <- NULL
if (!is.na(path_full_ref)) {
  full_reference <- readr::read_csv(path_full_ref, show_col_types = FALSE)
}

# output folders
fig_dir <- file.path(getwd(), "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
manuscript_fig_dir <- file.path(dirname(getwd()), "manuscript", "figures")

# plotting parameters -----------------------------------------------------
# 6 colorblind-friendly colors, matching the original normal visualization
cb_colors <- c(
  "#0072B2", # blue
  "#E69F00", # orange
  "#009E73", # bluish green
  "#CC79A7", # reddish purple
  "#56B4E9", # sky blue
  "#D55E00"  # vermillion
)

par_save_height = 4

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

copy_to_manuscript <- function(path) {
  if (dir.exists(manuscript_fig_dir)) {
    file.copy(path, manuscript_fig_dir, overwrite = TRUE, recursive = TRUE)
  }
}

# factor and labels used in all metric plots ------------------------------
# Keep only the step sizes used in the revised core experiment.  Some older
# assembled files may contain unused factor levels or NA rows; removing them
# here prevents spurious legend entries such as alpha = 0.75, alpha = 0.25,
# or NA from appearing in the final figures.
df_wasp <- df_wasp %>%
  dplyr::mutate(
    nsplit = as.integer(nsplit),
    support = as.integer(support),
    alpha = as.numeric(alpha)
  ) %>%
  dplyr::filter(!is.na(nsplit), !is.na(support), !is.na(alpha)) %>%
  dplyr::filter(alpha %in% c(1.0, 0.5)) %>%
  dplyr::mutate(
    support_lab = factor(support, levels = sort(unique(support))),
    alpha_lab = factor(ifelse(alpha == 1.0, "alpha = 1", "alpha = 0.5"),
                       levels = c("alpha = 1", "alpha = 0.5"))
  ) %>%
  droplevels()

support_levels <- levels(droplevels(df_wasp$support_lab))
# Use first several colors for support sizes.
support_cols <- setNames(cb_colors[seq_along(support_levels)], support_levels)
alpha_ltys <- c("alpha = 1" = "solid",
                "alpha = 0.5" = "dashed")
alpha_ltys <- alpha_ltys[names(alpha_ltys) %in% levels(df_wasp$alpha_lab)]

# figure 1 : representative posterior distributions -----------------------
# This retains the original qualitative display: analytic posterior,
# full-data MCMC, powered subset posteriors, and the Wasserstein barycenter.
# It is skipped if the reduced data file or matching run output is absent.
run_fig1 <- !is.na(path_data)

if (run_fig1) {
  save_path = file.path(fig_dir, "fig-sim-normal-1.png")

  load(path_data)

  # Accept both old and new object names.
  if (exists("common_posterior_mean")) analytic_mean <- common_posterior_mean
  if (exists("common_posterior_cov"))  analytic_cov  <- common_posterior_cov
  if (exists("common_full_mcmc"))      mcmc_full     <- common_full_mcmc
  if (exists("individual_subfits"))    mcmc_subsets  <- individual_subfits
  if (exists("individual_nsplits"))    mcmc_nsplits  <- individual_nsplits

  needed_fig1_objects <- c("analytic_mean", "analytic_cov", "mcmc_full", "mcmc_subsets")
  missing_fig1 <- needed_fig1_objects[!vapply(needed_fig1_objects, exists, logical(1))]
  if (length(missing_fig1) > 0L) {
    warning("Skipping Fig. 1 because the data file is missing objects: ", paste(missing_fig1, collapse = ", "))
    run_fig1 <- FALSE
  }
}

if (run_fig1) {
  # Representative setting. Use a small number of subsets for a readable contour plot.
  preferred_rep     <- 1L
  preferred_nsplit  <- 5L
  preferred_support <- 100L
  preferred_alpha   <- 1.0

  # If the preferred setting is unavailable, choose the closest available setting.
  avail <- df_wasp %>%
    dplyr::mutate(score = abs(rep - preferred_rep) +
                    10*abs(nsplit - preferred_nsplit) +
                    abs(support - preferred_support)/100 +
                    100*abs(alpha - preferred_alpha)) %>%
    dplyr::arrange(score)
  sel_cfg <- avail[1, ]
  sel_rep <- as.integer(sel_cfg$rep)
  sel_nsplit <- as.integer(sel_cfg$nsplit)
  sel_support <- as.integer(sel_cfg$support)
  sel_alpha <- as.numeric(sel_cfg$alpha)

  # Locate the subset posterior list in the reduced data.
  # New reduced data were indexed by expand.grid(rep = 1:50, nsplits = 2:20).
  if (exists("individual_subfits")) {
    original_grid <- expand.grid(rep = seq_len(50L), nsplits = seq(from = 2L, to = 20L, by = 1L))
    target_idx <- which(original_grid$rep == sel_rep & original_grid$nsplits == sel_nsplit)
    if (length(target_idx) != 1L) target_idx <- which(mcmc_nsplits == sel_nsplit)[1L]
    sel_subfit <- mcmc_subsets[[target_idx]]
  } else {
    target_idx <- which(mcmc_nsplits == sel_nsplit)[1L]
    sel_subfit <- mcmc_subsets[[target_idx]]
  }

  # Locate the corresponding barycenter support from the new result files.
  run_dir_candidates <- c(
    file.path(getwd(), "simulation-normal-extended", paste0("runs-", run_label)),
    file.path(getwd(), "simulation-normal-extended", "runs-core"),
    file.path(getwd(), "simulation-normal-extended", "runs"),
    file.path(getwd(), paste0("runs-", run_label)),
    file.path(getwd(), "runs-core"),
    file.path(getwd(), "runs")
  )
  run_dir_candidates <- run_dir_candidates[dir.exists(run_dir_candidates)]

  sel_bary <- NULL
  if (length(run_dir_candidates) > 0L) {
    result_files <- sort(unique(unlist(lapply(run_dir_candidates, function(dd) {
      list.files(dd, pattern = "^result_.*\\.RData$", full.names = TRUE)
    }))))

    # First pass: exact match.
    for (ff in result_files) {
      env <- new.env(parent = emptyenv())
      load(ff, envir = env)
      if (exists("metadata", envir = env) && exists("est_bary", envir = env)) {
        md <- get("metadata", envir = env)
        if (isTRUE(as.integer(md$rep) == sel_rep) &&
            isTRUE(as.integer(md$nsplit) == sel_nsplit) &&
            isTRUE(as.integer(md$support) == sel_support) &&
            isTRUE(abs(as.numeric(md$alpha) - sel_alpha) < 1e-12)) {
          eb <- get("est_bary", envir = env)
          sel_bary <- eb$support
          break
        }
      }
    }
  }

  if (is.null(sel_bary)) {
    warning("Skipping Fig. 1 because no matching barycenter result file was found. Metric figures will still be generated.")
  } else {
    # common parametric setting: preserve original visual window.
    par_xlim = c(0.63, 1.13)
    par_ylim = c(0.77, 1.27)

    # common grid
    grid_x <- seq(from=par_xlim[1], to=par_xlim[2], length.out=400)
    grid_y <- seq(from=par_ylim[1], to=par_ylim[2], length.out=400)
    grid <- expand.grid(x = grid_x, y = grid_y)

    # (a) analytic posterior distribution
    grid_1a <- grid
    grid_1a$z_analytic <- mvtnorm::dmvnorm(grid_1a,
                                           mean=analytic_mean,
                                           sigma=analytic_cov)

    # (b) full-data MCMC posterior KDE
    kde_full <- MASS::kde2d(mcmc_full[,1], mcmc_full[,2],
                            lims=c(range(grid_x), range(grid_y)))
    kde_b <- expand.grid(x=kde_full$x, y=kde_full$y)
    kde_b$z_kde <- as.vector(kde_full$z); rm(kde_full)

    # (c) powered subset posterior contours
    if (is.null(names(sel_subfit))) {
      names(sel_subfit) <- paste0("Subset ", seq_along(sel_subfit))
    }
    df_sub <- do.call(rbind, lapply(seq_along(sel_subfit), function(i) {
      m <- as.matrix(sel_subfit[[i]])
      data.frame(x = m[,1], y = m[,2], subset = names(sel_subfit)[i], row.names = NULL)
    }))
    df_sub$subset <- factor(df_sub$subset, levels = names(sel_subfit))

    # (d) barycenter KDE
    kde_bary <- MASS::kde2d(sel_bary[,1], sel_bary[,2],
                            lims=c(range(grid_x), range(grid_y)))
    kde_d <- expand.grid(x=kde_bary$x, y=kde_bary$y)
    kde_d$z_kde <- as.vector(kde_bary$z); rm(kde_bary)

    # common density scale for raster panels
    zlim = range(c(grid_1a$z_analytic, kde_b$z_kde, kde_d$z_kde), finite = TRUE)

    fig_1_a <- ggplot(grid_1a, aes(x, y)) +
      geom_raster(aes(fill = z_analytic)) +
      geom_contour(aes(z = z_analytic), color = "white", linewidth = 0.3) +
      scale_fill_viridis_c(name = "Density", limits = zlim) +
      theme_void() +
      theme(legend.position = "bottom",
            panel.border = element_rect(color = "gray50", fill = NA, linewidth = 0.5),
            plot.title = element_text(margin = margin(b = 6)))

    fig_1_b <- ggplot(kde_b, aes(x, y)) +
      geom_raster(aes(fill = z_kde)) +
      geom_contour(aes(z = z_kde), color = "white", linewidth = 0.3) +
      scale_fill_viridis_c(name = "Density", limits = zlim) +
      theme_void() +
      theme(legend.position = "bottom",
            panel.border = element_rect(color = "gray50", fill = NA, linewidth = 0.5),
            plot.title = element_text(margin = margin(b = 6)))

    # color-blind friendly set for subsets; repeat if needed.
    okabe <- rep(cb_colors, length.out = length(levels(df_sub$subset)))
    names(okabe) <- levels(df_sub$subset)
    subset_ltys <- rep(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),
                       length.out = length(levels(df_sub$subset)))
    names(subset_ltys) <- levels(df_sub$subset)

    fig_1_c <- ggplot(df_sub, aes(x = x, y = y, color = subset, linetype = subset)) +
      stat_density_2d(linewidth = 0.9, bins = 8) +
      scale_color_manual(values = okabe, name = NULL) +
      scale_linetype_manual(values = subset_ltys, name = NULL) +
      theme_void() +
      theme(legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "horizontal",
            panel.border = element_rect(color = "gray50", fill = NA, linewidth = 0.5),
            plot.title = element_text(margin = margin(b = 6))) +
      xlim(par_xlim) + ylim(par_ylim)

    fig_1_d <- ggplot(kde_d, aes(x, y)) +
      geom_raster(aes(fill = z_kde)) +
      geom_contour(aes(z = z_kde), color = "white", linewidth = 0.3) +
      scale_fill_viridis_c(name = "Density", limits = zlim) +
      theme_void() +
      theme(legend.position = "bottom",
            panel.border = element_rect(color = "gray50", fill = NA, linewidth = 0.5),
            plot.title = element_text(margin = margin(b = 6)))

    # compute half-cell padding from grid
    dx <- diff(range(grid_x)) / (length(unique(grid_x)) - 1)
    dy <- diff(range(grid_y)) / (length(unique(grid_y)) - 1)
    xlim_pad <- c(min(grid_x) - dx/2, max(grid_x) + dx/2)
    ylim_pad <- c(min(grid_y) - dy/2, max(grid_y) + dy/2)
    common_margin <- theme(plot.margin = margin(r = 5, l = 5))

    fig_1_a <- fig_1_a + coord_equal(xlim = xlim_pad, ylim = ylim_pad, expand = FALSE) + ggtitle("(A)") + common_margin
    fig_1_b <- fig_1_b + coord_equal(xlim = xlim_pad, ylim = ylim_pad, expand = FALSE) + ggtitle("(B)") + common_margin
    fig_1_c <- fig_1_c + coord_equal(xlim = xlim_pad, ylim = ylim_pad, expand = FALSE) + ggtitle("(C)") + common_margin
    fig_1_d <- fig_1_d + coord_equal(xlim = xlim_pad, ylim = ylim_pad, expand = FALSE) + ggtitle("(D)") + common_margin

    fig_1_final <- fig_1_a + fig_1_b + fig_1_c + fig_1_d +
      plot_layout(nrow = 1, guides = "collect") &
      theme(legend.position = "bottom")
    plot(fig_1_final)

    ggsave(filename=save_path,
           plot=fig_1_final,
           height=3.75,
           width=8,
           units="in")
    copy_to_manuscript(save_path)
  }
}

# figure 2 : accuracy, support size, and step size ------------------------
# Path to save
save_path = file.path(fig_dir, "fig-sim-normal-2.png")

# Summary data
error_semi  <- aux_ci(df_wasp, c("nsplit", "support_lab", "alpha_lab"), "semi_w2")
error_mmean <- aux_ci(df_wasp, c("nsplit", "support_lab", "alpha_lab"), "mean_error")
error_mcov  <- aux_ci(df_wasp, c("nsplit", "support_lab", "alpha_lab"), "cov_error")

full_semi  <- if (!is.null(full_reference) && "semi_w2" %in% names(full_reference)) full_reference$semi_w2[1] else NA_real_
full_mmean <- if (!is.null(full_reference) && "mean_error" %in% names(full_reference)) full_reference$mean_error[1] else NA_real_
full_mcov  <- if (!is.null(full_reference) && "cov_error" %in% names(full_reference)) full_reference$cov_error[1] else NA_real_

# Main accuracy plots. Colors encode support size; linetypes encode alpha.
fig_2_a <- ggplot(error_semi, aes(x=nsplit, y=mean,
                                  color=support_lab, linetype=alpha_lab,
                                  group=interaction(support_lab, alpha_lab))) +
  geom_line(linewidth=0.8) +
  geom_point(size=1.8) +
  geom_errorbar(aes(ymin=lower, ymax=upper), linewidth=0.4, alpha=0.75, width=0.15) +
  scale_color_manual(values = support_cols, limits = names(support_cols), breaks = names(support_cols), name = "Support size", drop = TRUE) +
  scale_linetype_manual(values = alpha_ltys, limits = names(alpha_ltys), breaks = names(alpha_ltys), name = "Step size", drop = TRUE) +
  labs(x = "Number of subsets", y="Error")
if (is.finite(full_semi)) {
  fig_2_a <- fig_2_a + geom_hline(yintercept = full_semi, color="black", linetype="dotted")
}

fig_2_b <- ggplot(error_mmean, aes(x=nsplit, y=mean,
                                   color=support_lab, linetype=alpha_lab,
                                   group=interaction(support_lab, alpha_lab))) +
  geom_line(linewidth=0.8) +
  geom_point(size=1.8) +
  geom_errorbar(aes(ymin=lower, ymax=upper), linewidth=0.4, alpha=0.75, width=0.15) +
  scale_color_manual(values = support_cols, limits = names(support_cols), breaks = names(support_cols), name = "Support size", drop = TRUE) +
  scale_linetype_manual(values = alpha_ltys, limits = names(alpha_ltys), breaks = names(alpha_ltys), name = "Step size", drop = TRUE) +
  labs(x = "Number of subsets", y="Error")
if (is.finite(full_mmean)) {
  fig_2_b <- fig_2_b + geom_hline(yintercept = full_mmean, color="black", linetype="dotted")
}

fig_2_c <- ggplot(error_mcov, aes(x=nsplit, y=mean,
                                  color=support_lab, linetype=alpha_lab,
                                  group=interaction(support_lab, alpha_lab))) +
  geom_line(linewidth=0.8) +
  geom_point(size=1.8) +
  geom_errorbar(aes(ymin=lower, ymax=upper), linewidth=0.4, alpha=0.75, width=0.15) +
  scale_color_manual(values = support_cols, limits = names(support_cols), breaks = names(support_cols), name = "Support size", drop = TRUE) +
  scale_linetype_manual(values = alpha_ltys, limits = names(alpha_ltys), breaks = names(alpha_ltys), name = "Step size", drop = TRUE) +
  labs(x = "Number of subsets", y="Error")
if (is.finite(full_mcov)) {
  fig_2_c <- fig_2_c + geom_hline(yintercept = full_mcov, color="black", linetype="dotted")
}

fig_2_a <- fig_2_a + ggtitle("(A)") + aux_theme()
fig_2_b <- fig_2_b + ggtitle("(B)") + aux_theme()
fig_2_c <- fig_2_c + ggtitle("(C)") + aux_theme()

fig_2_final <- fig_2_a + fig_2_b + fig_2_c +
  plot_layout(guides = "collect") &
  theme_classic() &
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.key.width = grid::unit(1, "cm"),
        aspect.ratio=1) &
  guides(
    color    = guide_legend(title = "Barycenter support size", nrow = 1, order = 1),
    linetype = guide_legend(title = "Step size", nrow = 1, order = 2)
  )
plot(fig_2_final)

ggsave(filename=save_path,
       plot=fig_2_final,
       height=3.5,
       width=9,
       units="in")
copy_to_manuscript(save_path)

# figure 3 : conditional-mean, interval-inclusion, and runtime diagnostics -
# Legacy pred_* CSV fields concern x^T beta and omit future-response noise.
# Inclusion is descriptive, not predictive or repeated-data calibration.
save_path = file.path(fig_dir, "fig-sim-normal-3.png")

pred_rmse <- aux_ci(df_wasp, c("nsplit", "support_lab", "alpha_lab"), "pred_rmse_mean")
pred_cov  <- aux_ci(df_wasp, c("nsplit", "support_lab", "alpha_lab"), "pred_interval_coverage")
post_cov  <- aux_ci(df_wasp, c("nsplit", "support_lab", "alpha_lab"), "marginal_coverage")
runtime_supp <- aux_ci(df_wasp, c("support", "support_lab", "alpha_lab"), "runtime_sec")

full_pred_rmse <- if (!is.null(full_reference) && "pred_rmse_mean" %in% names(full_reference)) full_reference$pred_rmse_mean[1] else NA_real_
full_pred_cov  <- if (!is.null(full_reference) && "pred_interval_coverage" %in% names(full_reference)) full_reference$pred_interval_coverage[1] else NA_real_
full_post_cov  <- if (!is.null(full_reference) && "marginal_coverage" %in% names(full_reference)) full_reference$marginal_coverage[1] else NA_real_

fig_3_a <- ggplot(pred_rmse, aes(x=nsplit, y=mean,
                                 color=support_lab, linetype=alpha_lab,
                                 group=interaction(support_lab, alpha_lab))) +
  geom_line(linewidth=0.8) + geom_point(size=1.8) +
  geom_errorbar(aes(ymin=lower, ymax=upper), linewidth=0.4, alpha=0.75, width=0.15) +
  scale_color_manual(values = support_cols, limits = names(support_cols), breaks = names(support_cols), name = "Support size", drop = TRUE) +
  scale_linetype_manual(values = alpha_ltys, limits = names(alpha_ltys), breaks = names(alpha_ltys), name = "Step size", drop = TRUE) +
  labs(x = "Number of subsets", y = "Conditional-mean RMSE")
if (is.finite(full_pred_rmse)) fig_3_a <- fig_3_a + geom_hline(yintercept = full_pred_rmse, color="black", linetype="dotted")

fig_3_b <- ggplot(pred_cov, aes(x=nsplit, y=mean,
                                color=support_lab, linetype=alpha_lab,
                                group=interaction(support_lab, alpha_lab))) +
  geom_line(linewidth=0.8) + geom_point(size=1.8) +
  geom_errorbar(aes(ymin=lower, ymax=upper), linewidth=0.4, alpha=0.75, width=0.15) +
  geom_hline(yintercept = 0.95, color="gray40", linetype="dashed") +
  scale_color_manual(values = support_cols, limits = names(support_cols), breaks = names(support_cols), name = "Support size", drop = TRUE) +
  scale_linetype_manual(values = alpha_ltys, limits = names(alpha_ltys), breaks = names(alpha_ltys), name = "Step size", drop = TRUE) +
  labs(x = "Number of subsets", y = "Conditional-mean inclusion")
if (is.finite(full_pred_cov)) fig_3_b <- fig_3_b + geom_hline(yintercept = full_pred_cov, color="black", linetype="dotted")

fig_3_c <- ggplot(post_cov, aes(x=nsplit, y=mean,
                                color=support_lab, linetype=alpha_lab,
                                group=interaction(support_lab, alpha_lab))) +
  geom_line(linewidth=0.8) + geom_point(size=1.8) +
  geom_errorbar(aes(ymin=lower, ymax=upper), linewidth=0.4, alpha=0.75, width=0.15) +
  geom_hline(yintercept = 0.95, color="gray40", linetype="dashed") +
  scale_color_manual(values = support_cols, limits = names(support_cols), breaks = names(support_cols), name = "Support size", drop = TRUE) +
  scale_linetype_manual(values = alpha_ltys, limits = names(alpha_ltys), breaks = names(alpha_ltys), name = "Step size", drop = TRUE) +
  labs(x = "Number of subsets", y = "Coefficient inclusion")
if (is.finite(full_post_cov)) fig_3_c <- fig_3_c + geom_hline(yintercept = full_post_cov, color="black", linetype="dotted")

fig_3_d <- ggplot(runtime_supp, aes(x=support_lab, y=mean,
                                    linetype=alpha_lab, group=alpha_lab)) +
  geom_line(color=cb_colors[1], linewidth=0.8) +
  geom_point(color=cb_colors[1], size=1.8) +
  geom_errorbar(aes(ymin=lower, ymax=upper), color=cb_colors[1], linewidth=0.4, alpha=0.75, width=0) +
  scale_linetype_manual(values = alpha_ltys, limits = names(alpha_ltys), breaks = names(alpha_ltys), name = "Step size", drop = TRUE) +
  labs(x = "Support size", y = "Runtime (sec)")

fig_3_a <- fig_3_a + ggtitle("(A)") + aux_theme()
fig_3_b <- fig_3_b + ggtitle("(B)") + aux_theme()
fig_3_c <- fig_3_c + ggtitle("(C)") + aux_theme()
fig_3_d <- fig_3_d + ggtitle("(D)") + aux_theme(angle_x = 45)

# Use a single-row layout for the diagnostic figure.
# The legends are collected once at the bottom.
fig_3_final <- fig_3_a + fig_3_b + fig_3_c + fig_3_d +
  plot_layout(nrow = 1, guides = "collect") &
  theme_classic() &
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.key.width = grid::unit(1, "cm"),
        aspect.ratio = 1) &
  guides(
    color    = guide_legend(title = "Barycenter support size", nrow = 1, order = 1),
    linetype = guide_legend(title = "Step size", nrow = 1, order = 2)
  )
plot(fig_3_final)

ggsave(filename=save_path,
       plot=fig_3_final,
       height=3.6,
       width=13.2,
       units="in")
copy_to_manuscript(save_path)

# summary tables ----------------------------------------------------------
# Keep compact tables that are useful for manuscript text and response letter.
table_wasp_accuracy <- df_wasp %>%
  dplyr::group_by(nsplit, support, alpha) %>%
  dplyr::summarise(
    semi_w2_mean = mean(semi_w2, na.rm = TRUE),
    semi_w2_sd = stats::sd(semi_w2, na.rm = TRUE),
    mean_error_mean = mean(mean_error, na.rm = TRUE),
    cov_error_mean = mean(cov_error, na.rm = TRUE),
    pred_rmse_mean = mean(pred_rmse_mean, na.rm = TRUE),
    pred_coverage_mean = mean(pred_interval_coverage, na.rm = TRUE),
    marginal_coverage_mean = mean(marginal_coverage, na.rm = TRUE),
    runtime_mean = mean(runtime_sec, na.rm = TRUE),
    niter_mean = mean(niter, na.rm = TRUE),
    monotone_violations_total = sum(monotone_violations, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(table_wasp_accuracy,
          file = file.path(fig_dir, "table-sim-normal-wasp-summary.csv"),
          row.names = FALSE)

if (!is.null(full_reference)) {
  write.csv(full_reference,
            file = file.path(fig_dir, "table-sim-normal-full-reference.csv"),
            row.names = FALSE)
}

message("Saved WASP figures and tables under: ", fig_dir)
