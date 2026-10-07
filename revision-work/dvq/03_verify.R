# Verify retained run histories, assignments, metrics, and historical audit.
# Recomputes only final-state OT and evaluation, not local fits or the full study.
script_arg <- grep('^--file=', commandArgs(), value = TRUE)
script_arg <- gsub("~+~", " ", script_arg, fixed=TRUE) # Decode spaces encoded by Rscript.
stopifnot(length(script_arg) == 1L)
root <- dirname(normalizePath(sub('^--file=', '', script_arg)))
repo_root <- normalizePath(file.path(root, '../..'))
for (pkg in c('T4transport', 'mclustcomp', 'jsonlite')) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop('Required R package: ', pkg)
}
source(file.path(root, 'common.R'))
metrics <- read.csv(file.path(root, 'paired_metrics.csv'))
files <- sort(list.files(file.path(root, 'runs'), '\\.rds$', full.names = TRUE))
stopifnot(length(files) == 12L, nrow(metrics) == 24L, all(metrics$converged))
max_feasibility <- max_residual <- max_native_difference <- 0
updates <- integer(0)
for (file in files) {
  key <- sub('\\.rds$', '', basename(file)); bits <- strsplit(key, '_', fixed = TRUE)[[1]]
  name <- bits[1]; S <- as.integer(sub('S', '', bits[2])); repetition <- as.integer(sub('rep', '', bits[3]))
  result <- readRDS(file); input <- readRDS(file.path(root, 'inputs', basename(file)))
  stopifnot(!grepl('^/', input$data_file), identical(input$initial, result$dvq$initial))
  data_file <- file.path(repo_root, input$data_file)
  stopifnot(unname(tools::md5sum(data_file)) == input$data_md5)
  e <- new.env(); load(data_file, e)
  h <- result$dvq$history; terminal <- tail(h, 1L)
  stopifnot(result$dvq$converged,
            terminal$relative_change <= result$dvq$objective_tol,
            terminal$normalized_residual <= result$dvq$residual_tol,
            all(diff(h$objective) <= 1e-10 * (1 + abs(head(h$objective, -1L)))))
  final <- exact_state(result$dvq$centers, input$summaries)
  stopifnot(abs(final$objective - terminal$objective) <= 1e-10 * max(1, terminal$objective))
  max_feasibility <- max(max_feasibility, final$feasible_error)
  max_residual <- max(max_residual, sqrt(final$residual2 / result$dvq$scale2))
  max_native_difference <- max(max_native_difference, abs(result$dvq$centers - result$native$support))
  updates <- c(updates, result$dvq$iterations)
  for (method in c('DVQ', 'Summary k-means')) {
    centers <- if (method == 'DVQ') result$dvq$centers else result$kmeans$centers
    score <- evaluate(as.matrix(e$X_reduced), e$Y, centers)
    expected <- metrics[metrics$data == name & metrics$S == S & metrics$repeat_id == repetition & metrics$method == method, ]
    stopifnot(nrow(expected) == 1L, identical(score$pred, result$scores[[method]]$pred),
              expected$transmitted_coordinates == S * 10L * 20L)
    for (measure in c('distortion', 'ARI', 'NMI')) {
      stopifnot(abs(score[[measure]] - expected[[measure]]) <= 1e-10 * max(1, abs(expected[[measure]])))
    }
  }
}
audit <- read.csv(file.path(root, 'historical_cap_audit.csv'))
reproduction <- read.csv(file.path(root, 'reproduction_checks.csv'))
stopifnot(sum(audit$n_available) == 64L, sum(audit$n_at_cap) == 0L,
          max(audit$max_niter) == 22L, sum(reproduction$historical_record_available) == 8L,
          all(reproduction$historical_assignments_identical[reproduction$historical_record_available]),
          all(reproduction$max_historical_center_abs_diff[reproduction$historical_record_available] == 0),
          min(updates) == 3L, max(updates) == 5L, max_feasibility < 1e-12, max_residual <= 1e-6)
# Independently compare the eight archived matching historical records when present.
historical <- read.csv(file.path(root, 'available_historical_metrics.csv'))
legacy_dir <- file.path(repo_root, 'code-real-clustering/real-2-clustering-extended/outcome-core')
n_legacy_checked <- 0L
for (key in reproduction$key[reproduction$historical_record_available]) {
  bits <- strsplit(key, '_', fixed = TRUE)[[1]]
  input <- readRDS(file.path(root, 'inputs', paste0(key, '.rds')))
  matching <- historical[historical$data == bits[1] & historical$seed == input$seed & historical$base_algorithm == 'DVQ', ]
  stopifnot(nrow(matching) == 1L)
  file <- file.path(legacy_dir, matching$file)
  if (file.exists(file)) {
    e <- new.env(); load(file, e)
    run <- readRDS(file.path(root, 'runs', paste0(key, '.rds')))
    stopifnot(identical(e$centers, run$native$support),
              identical(as.integer(e$pred_label), run$scores$DVQ$pred))
    n_legacy_checked <- n_legacy_checked + 1L
  }
}
cat(jsonlite::toJSON(list(passed = TRUE, n_runs = length(files),
  update_range = range(updates), max_plan_feasibility_error = max_feasibility,
  max_normalized_terminal_residual = max_residual,
  max_difference_from_native_support = max_native_difference,
  historical_records_independently_checked = n_legacy_checked), pretty = TRUE, auto_unbox = TRUE), '\n')
