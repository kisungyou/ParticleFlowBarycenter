# Install dependencies for the current replication scripts and original notebooks.
# Archived R results used R 4.5.1 and T4transport 0.1.9.

cran_packages <- c(
  "pacman", "ggplot2", "patchwork", "dplyr", "readr", "tidyr",
  "mvtnorm", "MASS", "abind", "skmeans", "cluster", "clusterCrit",
  "mclustcomp", "remotes", "jsonlite", "rstudioapi"
)
missing <- cran_packages[!vapply(cran_packages, requireNamespace,
                               quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing) > 0L) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}
if (!requireNamespace("T4transport", quietly = TRUE)) {
  remotes::install_github("kisungyou/T4transport")
}
if (as.character(packageVersion("T4transport")) != "0.1.9") {
  warning("Archived computations used T4transport 0.1.9. Review version differences before rerunning.")
}
if (!"alpha" %in% names(formals(T4transport::rbaryGD))) {
  stop("The original damped workflows require an rbaryGD implementation with an alpha argument. See revision-work/refresh/README.md.")
}
message("Dependencies are available. Saved-result checks and figure builds do not require rerunning the native barycenter solver.")
