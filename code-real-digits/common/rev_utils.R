# Print a traceback in Rscript logs before exiting on uncaught errors.
options(error = function() {
  traceback(2)
  q(save = "no", status = 1, runLast = FALSE)
})

# Shared utilities for revision experiments.
# These functions are intentionally dependency-light. Experiment scripts load
# heavier packages only when needed.

rev_arg_id <- function(default = 1L) {
  if (interactive()) return(as.integer(default))
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0L) return(as.integer(default))
  as.integer(as.numeric(args[1L]))
}

rev_ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}


rev_is_smoke <- function() {
  tolower(Sys.getenv("REV_SMOKE", unset = "0")) %in% c("1", "true", "yes", "y", "smoke")
}

rev_run_label <- function() {
  Sys.getenv("REV_RUN_LABEL", unset = "")
}

rev_label_suffix <- function() {
  lab <- rev_run_label()
  if (nzchar(lab)) paste0("-", lab) else ""
}

rev_output_dir <- function(script_dir, base) {
  file.path(script_dir, paste0(base, rev_label_suffix()))
}

rev_output_file <- function(script_dir, stem, ext = ".RData") {
  file.path(script_dir, paste0(stem, rev_label_suffix(), ext))
}

rev_smoke_int <- function(full, smoke) {
  if (rev_is_smoke()) as.integer(smoke) else as.integer(full)
}

rev_smoke_num <- function(full, smoke) {
  if (rev_is_smoke()) as.numeric(smoke) else as.numeric(full)
}

rev_require <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "),
         ". Install them before running this experiment.", call. = FALSE)
  }
  invisible(TRUE)
}

rev_elapsed <- function(expr) {
  t0 <- proc.time()
  value <- eval.parent(substitute(expr))
  elapsed <- as.numeric((proc.time() - t0)["elapsed"])
  list(value = value, elapsed = elapsed)
}

rev_history_summary <- function(history, tol = 1e-10) {
  history <- as.numeric(history)
  if (length(history) == 0L || all(!is.finite(history))) {
    return(list(first = NA_real_, last = NA_real_, min = NA_real_,
                n = 0L, monotone_violations = NA_integer_,
                max_increase = NA_real_))
  }
  diffs <- diff(history)
  allowed <- tol * (1 + abs(head(history, -1L)))
  list(
    first = history[1L],
    last = history[length(history)],
    min = min(history, na.rm = TRUE),
    n = length(history),
    monotone_violations = sum(diffs > allowed, na.rm = TRUE),
    max_increase = if (length(diffs) == 0L) 0 else max(c(0, diffs), na.rm = TRUE)
  )
}

rev_bures <- function(A, B) {
  A <- (A + t(A)) / 2
  B <- (B + t(B)) / 2
  ea <- eigen(A, symmetric = TRUE)
  Asqrt <- ea$vectors %*% (sqrt(pmax(ea$values, 0)) * t(ea$vectors))
  M <- Asqrt %*% B %*% Asqrt
  M <- (M + t(M)) / 2
  em <- eigen(M, symmetric = TRUE)
  tr_Msqrt <- sum(sqrt(pmax(em$values, 0)))
  out <- sqrt(max(0, sum(diag(A)) + sum(diag(B)) - 2 * tr_Msqrt))
  as.numeric(Re(out))
}

rev_empirical_objective <- function(support, measures, weights = NULL) {
  rev_require(c("T4transport"))
  if (is.null(weights)) weights <- rep(1 / length(measures), length(measures))
  vals <- numeric(length(measures))
  for (i in seq_along(measures)) {
    vals[i] <- T4transport::wasserstein(support, measures[[i]], p = 2)$distance^2
  }
  sum(weights * vals)
}

rev_make_gaussian_measures <- function(seed, num_measures = 4L, n_per_measure = 500L,
                                       d = 2L, loc_scale = 10, wishart_df = NULL,
                                       loc_sd = 1) {
  rev_require(c("mvtnorm", "abind"))
  set.seed(seed)
  if (is.null(wishart_df)) wishart_df <- max(d + 2L, 4L)
  means <- matrix(0, nrow = num_measures, ncol = d)
  # Symmetric locations for d >= 2; random signs for additional measures/dims.
  base_signs <- expand.grid(rep(list(c(-1, 1)), min(d, ceiling(log2(max(2, num_measures))))))
  base_signs <- as.matrix(base_signs)
  if (ncol(base_signs) < d) {
    extra <- matrix(sample(c(-1, 1), nrow(base_signs) * (d - ncol(base_signs)), replace = TRUE),
                    nrow = nrow(base_signs))
    base_signs <- cbind(base_signs, extra)
  }
  for (i in seq_len(num_measures)) {
    s <- base_signs[((i - 1L) %% nrow(base_signs)) + 1L, seq_len(d)]
    means[i, ] <- loc_scale * s + rnorm(d, sd = loc_sd)
  }
  covs <- array(0, dim = c(d, d, num_measures))
  measures <- vector("list", num_measures)
  for (i in seq_len(num_measures)) {
    covs[, , i] <- stats::rWishart(1, wishart_df, diag(d))[, , 1]
    measures[[i]] <- mvtnorm::rmvnorm(n = n_per_measure, mean = means[i, ], sigma = covs[, , i])
  }
  list(means = means, covs = covs, measures = measures)
}

rev_classification_metrics <- function(predicted, actual) {
  predicted <- as.integer(predicted)
  actual <- as.integer(actual)
  labs <- sort(unique(c(predicted, actual)))
  cm <- table(factor(actual, levels = labs), factor(predicted, levels = labs))
  acc <- sum(diag(cm)) / sum(cm)
  precision <- recall <- f1 <- rep(NA_real_, length(labs))
  for (i in seq_along(labs)) {
    tp <- cm[i, i]
    fp <- sum(cm[, i]) - tp
    fn <- sum(cm[i, ]) - tp
    precision[i] <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
    recall[i] <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
    f1[i] <- if (!is.finite(precision[i] + recall[i]) || (precision[i] + recall[i]) == 0) {
      NA_real_
    } else {
      2 * precision[i] * recall[i] / (precision[i] + recall[i])
    }
  }
  c(Accuracy = acc,
    MacroPrecision = mean(precision, na.rm = TRUE),
    MacroRecall = mean(recall, na.rm = TRUE),
    MacroF1 = mean(f1, na.rm = TRUE))
}

rev_otsu <- function(x, num_bins = 256L) {
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

rev_bin2coords <- function(mat) {
  dims <- dim(mat)
  p <- dims[1L]
  q <- dims[2L]
  idx <- which(abs(mat) > sqrt(.Machine$double.eps), arr.ind = TRUE)
  if (nrow(idx) == 0L) return(matrix(c(0, 0), ncol = 2, dimnames = list(NULL, c("x", "y"))))
  x_centered <- 2 * (idx[, "col"] - 1) / (q - 1) - 1
  y_centered <- 2 * (idx[, "row"] - 1) / (p - 1) - 1
  cbind(x = x_centered, y = y_centered)
}

rev_img2coords <- function(mat, num_bins = 256L) {
  image_thr <- rev_otsu(mat, num_bins = num_bins)
  image_bin <- round(as.matrix(mat > image_thr))
  rev_bin2coords(image_bin)
}

rev_normalized_image_vector <- function(mat) {
  x <- as.vector(as.numeric(mat))
  sx <- sum(x)
  if (!is.finite(sx) || sx <= 0) return(x)
  x / sx
}

rev_balanced_indices <- function(labels, per_class, seed = 1L) {
  set.seed(seed)
  labs <- sort(unique(labels))
  idx <- integer(0)
  for (lab in labs) {
    pool <- which(labels == lab)
    if (length(pool) < per_class) {
      warning("Class ", lab, " has only ", length(pool), " observations; using all of them.")
      idx <- c(idx, pool)
    } else {
      idx <- c(idx, sample(pool, per_class, replace = FALSE))
    }
  }
  sample(idx, length(idx), replace = FALSE)
}

rev_kmeans_centers <- function(X, k, seed = 1L, nstart = 10L, max_n = 5000L) {
  set.seed(seed)
  if (nrow(X) > max_n) X_fit <- X[sample(seq_len(nrow(X)), max_n, replace = FALSE), , drop = FALSE] else X_fit <- X
  fit <- stats::kmeans(X_fit, centers = k, nstart = nstart, iter.max = 100)
  fit$centers
}

rev_assign_to_centers <- function(X, centers) {
  # Avoid constructing an enormous full distance matrix when possible.
  out <- integer(nrow(X))
  chunk <- 2000L
  for (start in seq(1L, nrow(X), by = chunk)) {
    end <- min(nrow(X), start + chunk - 1L)
    XX <- X[start:end, , drop = FALSE]
    dmat <- as.matrix(stats::dist(rbind(XX, centers)))[seq_len(nrow(XX)),
                                                        nrow(XX) + seq_len(nrow(centers)), drop = FALSE]
    out[start:end] <- max.col(-dmat, ties.method = "first")
  }
  out
}

rev_script_dir <- function() {
  # Works for Rscript; falls back to getwd() in interactive use.
  frames <- sys.frames()
  files <- vapply(frames, function(x) {
    f <- attr(x, "ofile")
    if (is.null(f)) "" else f
  }, character(1))
  files <- files[nzchar(files)]
  if (length(files) > 0L) return(dirname(normalizePath(files[length(files)])))
  getwd()
}
