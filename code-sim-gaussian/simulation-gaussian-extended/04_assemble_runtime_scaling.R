# Assemble G2 runtime scaling results.
script_dir <- getwd()
source(file.path(dirname(script_dir), "common", "rev_utils.R"))
run_dir <- rev_output_dir(script_dir, "runs-scaling")
files <- sort(list.files(run_dir, pattern = "^result_.*\\.RData$", full.names = TRUE))
if (length(files) == 0L) stop("No result files found in ", run_dir, call. = FALSE)
rows <- vector("list", length(files))
for (i in seq_along(files)) {
  load(files[i])
  rows[[i]] <- metrics
}
df_runtime_scaling <- do.call(rbind, rows)
df_runtime_scaling <- df_runtime_scaling[order(df_runtime_scaling$sweep,
                                               df_runtime_scaling$rep,
                                               df_runtime_scaling$N,
                                               df_runtime_scaling$target_n,
                                               df_runtime_scaling$d,
                                               df_runtime_scaling$support,
                                               df_runtime_scaling$alpha), ]
save(df_runtime_scaling, file = rev_output_file(script_dir, "assembled_runtime_scaling"))
write.csv(df_runtime_scaling, file = rev_output_file(script_dir, "assembled_runtime_scaling", ext = ".csv"), row.names = FALSE)
