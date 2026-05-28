# Helper routines for Quarto replication notebooks.

ensure_pacman <- function(install = FALSE) {
  if (!requireNamespace("pacman", quietly = TRUE)) {
    if (!install) {
      stop("The 'pacman' package is required. Install it or set install_packages <- TRUE.", call. = FALSE)
    }
    install.packages("pacman", repos = "https://cloud.r-project.org")
  }
}

load_or_install <- function(pkgs, install = FALSE) {
  ensure_pacman(install = install)
  if (install) {
    github_pkgs <- intersect(pkgs, "T4transport")
    cran_pkgs <- setdiff(pkgs, github_pkgs)
    if (length(cran_pkgs) > 0) {
      pacman::p_load(char = cran_pkgs, install = TRUE, update = FALSE, character.only = TRUE)
    }
    if (length(github_pkgs) > 0 && !requireNamespace("T4transport", quietly = TRUE)) {
      if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes", repos = "https://cloud.r-project.org")
      }
      remotes::install_github("kisungyou/T4transport")
    }
  } else {
    missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
    if (length(missing) > 0) {
      stop(
        "Missing required R packages: ", paste(missing, collapse = ", "), "\n",
        "Install them manually or set install_packages <- TRUE in the first code cell.",
        call. = FALSE
      )
    }
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

set_run_environment <- function(run_mode, run_label = NULL, experiment_scale = "core") {
  run_mode <- match.arg(tolower(run_mode), c("smoke", "full"))
  if (is.null(run_label) || !nzchar(run_label)) {
    run_label <- if (run_mode == "smoke") "smoke" else "core"
  }
  Sys.setenv(
    REV_RUN_LABEL = run_label,
    REV_EXPERIMENT_SCALE = experiment_scale,
    REV_SMOKE = if (run_mode == "smoke") "1" else "0"
  )
  invisible(run_label)
}

make_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

append_log_tail <- function(log_file, n = 80L) {
  if (file.exists(log_file)) {
    x <- readLines(log_file, warn = FALSE)
    tail_x <- tail(x, n)
    cat(paste(tail_x, collapse = "\n"), "\n")
  }
}

run_rscript_task <- function(script, task_id, log_dir = "logs") {
  script <- normalizePath(script, mustWork = TRUE)
  script_dir <- dirname(script)
  script_base <- basename(script)
  make_dir(file.path(script_dir, log_dir))
  log_file <- file.path(script_dir, log_dir, sprintf("%s_%05d.log", tools::file_path_sans_ext(script_base), as.integer(task_id)))
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(script_dir)
  status <- system2(file.path(R.home("bin"), "Rscript"), c(script_base, as.character(task_id)), stdout = log_file, stderr = log_file)
  if (!identical(status, 0L)) {
    cat("\nTask failed: ", script, " id=", task_id, "\n", sep = "")
    cat("Log file: ", log_file, "\n", sep = "")
    append_log_tail(log_file)
    stop("Rscript task failed.", call. = FALSE)
  }
  invisible(log_file)
}

run_rscript_grid <- function(script, task_ids, jobs = 1L, label = NULL) {
  if (length(task_ids) == 0) return(invisible(NULL))
  jobs <- max(1L, as.integer(jobs))
  label <- label %||% basename(script)
  message("Running ", length(task_ids), " task(s) for ", label, " with jobs=", jobs, ".")
  if (.Platform$OS.type == "unix" && jobs > 1L) {
    parallel::mclapply(task_ids, function(id) run_rscript_task(script, id), mc.cores = jobs)
  } else {
    for (id in task_ids) run_rscript_task(script, id)
  }
  invisible(NULL)
}

run_rscript_once <- function(script, log_name = NULL) {
  script <- normalizePath(script, mustWork = TRUE)
  script_dir <- dirname(script)
  script_base <- basename(script)
  log_dir <- file.path(script_dir, "logs")
  make_dir(log_dir)
  if (is.null(log_name)) log_name <- paste0(tools::file_path_sans_ext(script_base), ".log")
  log_file <- file.path(log_dir, log_name)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(script_dir)
  status <- system2(file.path(R.home("bin"), "Rscript"), script_base, stdout = log_file, stderr = log_file)
  if (!identical(status, 0L)) {
    cat("\nScript failed: ", script, "\n", sep = "")
    cat("Log file: ", log_file, "\n", sep = "")
    append_log_tail(log_file)
    stop("Rscript failed.", call. = FALSE)
  }
  invisible(log_file)
}

run_python_task <- function(script, task_id, log_dir = "logs") {
  script <- normalizePath(script, mustWork = TRUE)
  script_dir <- dirname(script)
  script_base <- basename(script)
  make_dir(file.path(script_dir, log_dir))
  log_file <- file.path(script_dir, log_dir, sprintf("%s_%05d.log", tools::file_path_sans_ext(script_base), as.integer(task_id)))
  python <- Sys.which("python3")
  if (!nzchar(python)) python <- Sys.which("python")
  if (!nzchar(python)) stop("Python was not found on PATH.", call. = FALSE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(script_dir)
  status <- system2(python, c(script_base, as.character(task_id)), stdout = log_file, stderr = log_file)
  if (!identical(status, 0L)) {
    cat("\nPython task failed: ", script, " id=", task_id, "\n", sep = "")
    cat("Log file: ", log_file, "\n", sep = "")
    append_log_tail(log_file)
    stop("Python task failed. Check that POT is installed: python3 -m pip install -r ../requirements.txt", call. = FALSE)
  }
  invisible(log_file)
}

run_python_grid <- function(script, task_ids, jobs = 1L, label = NULL) {
  if (length(task_ids) == 0) return(invisible(NULL))
  jobs <- max(1L, as.integer(jobs))
  label <- label %||% basename(script)
  message("Running ", length(task_ids), " Python task(s) for ", label, " with jobs=", jobs, ".")
  if (.Platform$OS.type == "unix" && jobs > 1L) {
    parallel::mclapply(task_ids, function(id) run_python_task(script, id), mc.cores = jobs)
  } else {
    for (id in task_ids) run_python_task(script, id)
  }
  invisible(NULL)
}

require_file <- function(path, label = "required file") {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
  normalizePath(path)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
