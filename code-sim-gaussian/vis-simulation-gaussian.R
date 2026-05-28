# final visualization for revised Gaussian experiment
# This script follows the style of the original vis-simulation-gaussian.R:
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
               tidyr)

# Set working directory to the script location when run in RStudio; otherwise
# keep the current working directory, which should be the project root.
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# helper for paths: first look in the revision experiment folders, then in cwd
find_file <- function(...) {
  candidates <- c(...)
  for (p in candidates) {
    if (file.exists(p)) return(p)
  }
  stop("None of the candidate files exists: ", paste(candidates, collapse = ", "))
}

path_g1 <- find_file(
  file.path(getwd(), "simulation-gaussian-extended", "assembled_gaussian_resolution_alpha-core.csv"),
  file.path(getwd(), "assembled_gaussian_resolution_alpha-core.csv")
)
path_g2 <- find_file(
  file.path(getwd(), "simulation-gaussian-extended", "assembled_runtime_scaling-core.csv"),
  file.path(getwd(), "assembled_runtime_scaling-core.csv")
)
path_pot <- find_file(
  file.path(getwd(), "pot-baselines", "assembled_pot_gaussian_baselines-core.csv"),
  file.path(getwd(), "assembled_pot_gaussian_baselines-core.csv")
)

# load the data for plotting
df_g1  <- readr::read_csv(path_g1, show_col_types = FALSE)
df_g2  <- readr::read_csv(path_g2, show_col_types = FALSE)
df_pot <- readr::read_csv(path_pot, show_col_types = FALSE)

# output folders
fig_dir <- file.path(getwd(), "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# plotting parameters -----------------------------------------------------
# 4 colorblind-friendly colors
cb_colors <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7")

# size to save
par_save_height = 4

# small helper functions --------------------------------------------------
aux_ci <- function(data, group_vars, value_var) {
  data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      mean = mean(.data[[value_var]], na.rm = TRUE),
      sd = stats::sd(.data[[value_var]], na.rm = TRUE),
      n = sum(!is.na(.data[[value_var]])),
      se = sd/sqrt(n),
      lower = mean - 1.96*se,
      upper = mean + 1.96*se,
      .groups = "drop"
    )
}

aux_alpha_label <- function(x) {
  factor(paste0("alpha = ", x),
         levels = c("alpha = 1", "alpha = 0.5"))
}

aux_theme <- function() {
  theme_classic() +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# figure 2 : accuracy and step size ---------------------------------------
# This replaces/extends the original fig-sim-gauss-2 by showing both alpha=1
# and alpha=0.5. The style follows the original figure: line + point + errorbar.
save_path = file.path(fig_dir, "fig-sim-gauss-2.png")

support_breaks <- sort(unique(df_g1$support))

df_g1 <- df_g1 %>%
  mutate(alpha_lab = aux_alpha_label(alpha),
         support_lab = factor(support, levels = support_breaks))

df1_semi <- aux_ci(df_g1, c("support", "support_lab", "alpha_lab"), "semi_w2")
df2_moment_mean <- aux_ci(df_g1, c("support", "support_lab", "alpha_lab"), "mean_error")
df3_moment_cov <- aux_ci(df_g1, c("support", "support_lab", "alpha_lab"), "cov_error")

fig2_a <- ggplot(df1_semi, aes(x=support_lab, y=mean, color=alpha_lab, group=alpha_lab)) +
  geom_line(linewidth=1) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lower, ymax=upper),
                linewidth=0.5,
                alpha=0.5,
                width=0) +
  scale_color_manual(values = cb_colors[1:2], drop = FALSE) +
  labs(x = "Support size", y = "Error", color = NULL)

fig2_b <- ggplot(df2_moment_mean, aes(x=support_lab, y=mean, color=alpha_lab, group=alpha_lab)) +
  geom_line(linewidth=1) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lower, ymax=upper),
                linewidth=0.5,
                alpha=0.5,
                width=0) +
  scale_color_manual(values = cb_colors[1:2], drop = FALSE) +
  labs(x = "Support size", y = "Error", color = NULL)

fig2_c <- ggplot(df3_moment_cov, aes(x=support_lab, y=mean, color=alpha_lab, group=alpha_lab)) +
  geom_line(linewidth=1) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lower, ymax=upper),
                linewidth=0.5,
                alpha=0.5,
                width=0) +
  scale_color_manual(values = cb_colors[1:2], drop = FALSE) +
  labs(x = "Support size", y = "Error", color = NULL)

fig2_a <- fig2_a + ggtitle("(A)") + aux_theme()
fig2_b <- fig2_b + ggtitle("(B)") + aux_theme()
fig2_c <- fig2_c + ggtitle("(C)") + aux_theme()

fig2_final <- fig2_a + fig2_b + fig2_c +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
plot(fig2_final)

ggsave(filename=save_path,
       plot=fig2_final,
       height=3.4,
       width=9,
       unit="in")

# copy to manuscript/figures if that directory exists
manuscript_fig_dir <- file.path(dirname(getwd()), "manuscript", "figures")
if (dir.exists(manuscript_fig_dir)) {
  file.copy(save_path, manuscript_fig_dir, overwrite = TRUE, recursive=TRUE)
}

# figure 3 : runtime and scaling ------------------------------------------
save_path = file.path(fig_dir, "fig-sim-gauss-3.png")

make_runtime_summary <- function(sweep_name, xvar) {
  out <- df_g2 %>%
    filter(sweep == sweep_name) %>%
    aux_ci(c(xvar), "runtime_sec") %>%
    mutate(x = .data[[xvar]])
  out$x_lab <- factor(out$x, levels = sort(unique(out$x)))
  out
}

df_rt_m      <- make_runtime_summary("m", "support")
df_rt_target <- make_runtime_summary("target_n", "target_n")
df_rt_N      <- make_runtime_summary("N", "N")
df_rt_d      <- make_runtime_summary("d", "d")
df_rt_alpha  <- make_runtime_summary("alpha", "alpha")

fig3_a <- ggplot(df_rt_m, aes(x=x_lab, y=mean, group=1)) +
  geom_line(color=cb_colors[1], linewidth=1) +
  geom_point(color=cb_colors[1], size=2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), color=cb_colors[1], alpha=0.5, width=0) +
  labs(x = "Barycenter support size", y = "Runtime (sec)")

fig3_b <- ggplot(df_rt_target, aes(x=x_lab, y=mean, group=1)) +
  geom_line(color=cb_colors[1], linewidth=1) +
  geom_point(color=cb_colors[1], size=2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), color=cb_colors[1], alpha=0.5, width=0) +
  labs(x = "Target support size", y = "Runtime (sec)")

fig3_c <- ggplot(df_rt_N, aes(x=x_lab, y=mean, group=1)) +
  geom_line(color=cb_colors[1], linewidth=1) +
  geom_point(color=cb_colors[1], size=2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), color=cb_colors[1], alpha=0.5, width=0) +
  labs(x = "Number of measures", y = "Runtime (sec)")

fig3_d <- ggplot(df_rt_d, aes(x=x_lab, y=mean, group=1)) +
  geom_line(color=cb_colors[1], linewidth=1) +
  geom_point(color=cb_colors[1], size=2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), color=cb_colors[1], alpha=0.5, width=0) +
  labs(x = "Dimension", y = "Runtime (sec)")

fig3_e <- ggplot(df_rt_alpha, aes(x=x_lab, y=mean, group=1)) +
  geom_line(color=cb_colors[1], linewidth=1) +
  geom_point(color=cb_colors[1], size=2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), color=cb_colors[1], alpha=0.5, width=0.05) +
  labs(x = "Step size", y = "Runtime (sec)")

fig3_a <- fig3_a + ggtitle("(A)") + aux_theme()
fig3_b <- fig3_b + ggtitle("(B)") + aux_theme()
fig3_c <- fig3_c + ggtitle("(C)") + aux_theme()
fig3_d <- fig3_d + ggtitle("(D)") + aux_theme()
fig3_e <- fig3_e + ggtitle("(E)") + aux_theme()

fig3_final <- (fig3_a + fig3_b + fig3_c) / (fig3_d + fig3_e + plot_spacer())
plot(fig3_final)

ggsave(filename=save_path,
       plot=fig3_final,
       height=6.5,
       width=9,
       unit="in")
if (dir.exists(manuscript_fig_dir)) {
  file.copy(save_path, manuscript_fig_dir, overwrite = TRUE, recursive=TRUE)
}

# figure 4 : representative POT baselines ---------------------------------
# The proposed method values use alpha=1 from the G1 experiment at common
# support sizes. The POT values use the representative POT benchmark file.
save_path = file.path(fig_dir, "fig-sim-gauss-4.png")

common_supp <- c(50L, 100L, 200L)

df_ours <- df_g1 %>%
  filter(alpha == 1, support %in% common_supp) %>%
  transmute(rep, support,
            method_lab = "Proposed (alpha = 1)",
            semi_w2, mean_error, cov_error, runtime_sec)

df_pot2 <- df_pot %>%
  filter(status == "ok", support %in% common_supp) %>%
  mutate(method_lab = ifelse(method == "pot_exact",
                             "POT exact",
                             paste0("Sinkhorn, reg = ", reg))) %>%
  transmute(rep, support, method_lab, semi_w2, mean_error, cov_error, runtime_sec)

df_base <- bind_rows(df_ours, df_pot2) %>%
  mutate(support_lab = factor(support, levels = common_supp),
         method_lab = factor(method_lab,
                             levels = c("Proposed (alpha = 1)",
                                        "POT exact",
                                        "Sinkhorn, reg = 1",
                                        "Sinkhorn, reg = 5")))

method_cols <- c("Proposed (alpha = 1)" = cb_colors[1],
                 "POT exact" = cb_colors[2],
                 "Sinkhorn, reg = 1" = cb_colors[3],
                 "Sinkhorn, reg = 5" = cb_colors[4])

df_base_w2  <- aux_ci(df_base, c("support", "support_lab", "method_lab"), "semi_w2")
df_base_cov <- aux_ci(df_base, c("support", "support_lab", "method_lab"), "cov_error")
df_base_rt  <- aux_ci(df_base, c("support", "support_lab", "method_lab"), "runtime_sec")

fig4_a <- ggplot(df_base_w2, aes(x=support_lab, y=mean, color=method_lab, group=method_lab)) +
  geom_line(linewidth=1) + geom_point(size=2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), linewidth=0.5, alpha=0.5, width=0) +
  scale_color_manual(values=method_cols, drop=FALSE) +
  labs(x="Support size", y="Error", color=NULL)

fig4_b <- ggplot(df_base_cov, aes(x=support_lab, y=mean, color=method_lab, group=method_lab)) +
  geom_line(linewidth=1) + geom_point(size=2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), linewidth=0.5, alpha=0.5, width=0) +
  scale_color_manual(values=method_cols, drop=FALSE) +
  labs(x="Support size", y="Error", color=NULL)

fig4_c <- ggplot(df_base_rt, aes(x=support_lab, y=mean, color=method_lab, group=method_lab)) +
  geom_line(linewidth=1) + geom_point(size=2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), linewidth=0.5, alpha=0.5, width=0) +
  scale_color_manual(values=method_cols, drop=FALSE) +
  labs(x="Support size", y="Runtime (sec)", color=NULL)

fig4_a <- fig4_a + ggtitle("(A)") + aux_theme()
fig4_b <- fig4_b + ggtitle("(B)") + aux_theme()
fig4_c <- fig4_c + ggtitle("(C)") + aux_theme()

fig4_final <- fig4_a + fig4_b + fig4_c +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
plot(fig4_final)

ggsave(filename=save_path,
       plot=fig4_final,
       height=3.8,
       width=10,
       unit="in")
if (dir.exists(manuscript_fig_dir)) {
  file.copy(save_path, manuscript_fig_dir, overwrite = TRUE, recursive=TRUE)
}

# summary tables for manuscript text --------------------------------------
summary_g1 <- df_g1 %>%
  group_by(support, alpha) %>%
  summarise(
    semi_w2_mean = mean(semi_w2),
    cov_error_mean = mean(cov_error),
    runtime_mean = mean(runtime_sec),
    niter_mean = mean(niter),
    monotone_violations_total = sum(monotone_violations),
    .groups = "drop"
  )
readr::write_csv(summary_g1, file.path(fig_dir, "table-sim-gauss-resolution-alpha-summary.csv"))

summary_pot <- df_base %>%
  group_by(support, method_lab) %>%
  summarise(
    semi_w2_mean = mean(semi_w2),
    cov_error_mean = mean(cov_error),
    runtime_mean = mean(runtime_sec),
    .groups = "drop"
  )
readr::write_csv(summary_pot, file.path(fig_dir, "table-sim-gauss-pot-summary.csv"))
