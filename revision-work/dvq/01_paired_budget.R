# Replay the controlled identical-summary comparison; retained evidence is read-only.
# Rscript revision-work/dvq/01_paired_budget.R [--key pbmc_S5_rep1]
#   [--output-dir PATH] [--rebuild-inputs]
options(warn = 1)
script_arg <- grep('^--file=', commandArgs(), value = TRUE)
script_arg <- gsub("~+~", " ", script_arg, fixed=TRUE) # Decode spaces encoded by Rscript.
stopifnot(length(script_arg) == 1L)
script_dir <- dirname(normalizePath(sub('^--file=', '', script_arg)))
repo_root <- normalizePath(file.path(script_dir, '../..'))
args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  index <- match(flag, args)
  if (is.na(index)) return(default)
  if (index == length(args)) stop('Missing value for ', flag)
  args[index + 1L]
}
out <- arg_value('--output-dir', file.path(script_dir, 'recomputed'))
rebuild <- '--rebuild-inputs' %in% args
key_filter <- arg_value('--key', '')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
if (out == script_dir) stop('Choose a separate output directory to preserve retained evidence.')
for (pkg in c('T4transport', 'mclustcomp', 'jsonlite')) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop('Required R package: ', pkg)
}
source(file.path(repo_root, 'code-real-clustering/common/rev_utils.R'))
source(file.path(script_dir, 'common.R'))
dir.create(file.path(out, 'runs'), showWarnings = FALSE)
dir.create(file.path(out, 'inputs'), showWarnings = FALSE)
files <- sort(list.files(file.path(script_dir, 'inputs'), '\\.rds$', full.names = TRUE))
if (nzchar(key_filter)) files <- files[basename(files) == paste0(key_filter, '.rds')]
if (!length(files)) stop('No bundled input matches requested key.')
metrics <- list()
for (file in files) {
  key <- sub('\\.rds$', '', basename(file))
  bits <- strsplit(key, '_', fixed = TRUE)[[1]]
  name <- bits[1]; S <- as.integer(sub('S', '', bits[2]))
  repetition <- as.integer(sub('rep', '', bits[3]))
  input <- readRDS(file)
  data_path <- file.path(repo_root, input$data_file)
  if (!file.exists(data_path)) stop('Missing repository data: ', input$data_file)
  if (unname(tools::md5sum(data_path)) != input$data_md5) stop('Data checksum mismatch: ', input$data_file)
  e <- new.env(); load(data_path, e)
  X <- as.matrix(e$X_reduced); Y <- as.integer(e$Y)
  stopifnot(nrow(X) == length(Y), all(is.finite(X)))
  summaries <- input$summaries; k <- nrow(summaries[[1]])
  seed <- input$seed
  local_time <- NA_real_
  if (rebuild) {
    if (!exists('aux_ginit', asNamespace('T4transport'), inherits = FALSE)) {
      stop('--rebuild-inputs requires the retained T4transport Gaussian initializer aux_ginit. Use default replay with bundled initial atoms on other builds.')
    }
    local_time <- system.time({
      set.seed(seed)
      partition <- sample(rep(seq_len(S), length.out = nrow(X)))
      summaries <- lapply(seq_len(S), function(i) rev_kmeans_centers(
        X[partition == i, , drop = FALSE], k, seed = seed + 100L * i,
        nstart = 10L, max_n = 5000L))
    })[['elapsed']]
    set.seed(seed + 999L)
    input$initial <- get('aux_ginit', asNamespace('T4transport'))(summaries, k)
    input$partition <- partition; input$summaries <- summaries
  }
  saveRDS(input, file.path(out, 'inputs', basename(file)))
  dvq_time <- system.time({dvq <- exact_dvq(summaries, input$initial)})[['elapsed']]
  km_time <- system.time({
    set.seed(seed + 999L)
    km <- stats::kmeans(do.call(rbind, summaries), centers = k, nstart = 10L, iter.max = 100L)
  })[['elapsed']]
  results <- list(DVQ = list(centers = dvq$centers, runtime = dvq_time),
                  'Summary k-means' = list(centers = km$centers, runtime = km_time))
  scores <- list()
  for (method in names(results)) {
    result <- results[[method]]
    eval_time <- system.time({score <- evaluate(X, Y, result$centers)})[['elapsed']]
    scores[[method]] <- score
    metrics[[length(metrics) + 1L]] <- data.frame(
      data = name, S = S, k = k, repeat_id = repetition, seed = seed, method = method,
      n = nrow(X), d = ncol(X), transmitted_coordinates = S * k * ncol(X),
      distortion = score$distortion, ARI = score$ARI, NMI = score$NMI,
      summary_barycenter_objective = exact_state(result$centers, summaries)$objective,
      local_summary_sec = local_time, aggregation_sec = result$runtime,
      construction_sec = local_time + result$runtime, evaluation_sec = eval_time,
      converged = if (method == 'DVQ') dvq$converged else km$ifault == 0L,
      iterations = if (method == 'DVQ') dvq$iterations else km$iter,
      terminal_residual = if (method == 'DVQ') tail(dvq$history$normalized_residual, 1L) else NA_real_)
  }
  saveRDS(list(dvq = dvq, kmeans = km, scores = scores), file.path(out, 'runs', basename(file)))
  write.csv(do.call(rbind, metrics), file.path(out, 'paired_metrics.csv'), row.names = FALSE)
  cat('Completed', key, ':', dvq$iterations, 'updates, joint stopping rule', dvq$converged, '\n')
}
jsonlite::write_json(list(
  experiment = 'central-k center-only identical-summary paired check',
  mode = if (rebuild) 'rebuild local summaries and initializer' else 'replay saved summaries and initial atoms',
  keys = sub('\\.rds$', '', basename(files)),
  summary_weights = 'uniform local, uniform outer; no mass counts transmitted',
  communication = 'Skd real coordinates, excluding identifiers',
  timing = 'Aggregation and evaluation timed separately. Local/construction times are unavailable in replay mode; rebuild mode recomputes them. Retained timing evidence is not overwritten.',
  system = as.list(Sys.info()[c('sysname', 'release', 'version', 'machine')]),
  package_versions = lapply(c('T4transport', 'mclustcomp', 'jsonlite'), function(pkg)
    list(package = pkg, version = as.character(packageVersion(pkg)))),
  source_helpers = 'code-real-clustering/common/rev_utils.R'),
  file.path(out, 'metadata.json'), pretty = TRUE, auto_unbox = TRUE)
