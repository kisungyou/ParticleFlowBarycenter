# Install R packages used across the replication notebooks.

cran_packages <- c(
  "pacman", "ggplot2", "patchwork", "dplyr", "readr", "tidyr",
  "mvtnorm", "MASS", "abind", "skmeans", "cluster", "clusterCrit",
  "mclustcomp", "remotes"
)

missing <- cran_packages[!vapply(cran_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("T4transport", quietly = TRUE)) {
  remotes::install_github("kisungyou/T4transport")
}

message("Package installation/check complete.")
