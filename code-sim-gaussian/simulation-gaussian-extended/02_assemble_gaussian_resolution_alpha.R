# Assemble G1 results into a compact data frame for later visualization/tables.
script_dir <- getwd()
source(file.path(dirname(script_dir), "common", "rev_utils.R"))

run_dir <- rev_output_dir(script_dir, "runs")
files <- sort(list.files(run_dir, pattern = "^result_.*\\.RData$", full.names = TRUE))
if (length(files) == 0L) stop("No result files found in ", run_dir, call. = FALSE)

rows <- vector("list", length(files))
for (i in seq_along(files)) {
  load(files[i])
  rows[[i]] <- metrics
}
df_gaussian_resolution_alpha <- do.call(rbind, rows)
df_gaussian_resolution_alpha <- df_gaussian_resolution_alpha[order(df_gaussian_resolution_alpha$rep,
                                                                   df_gaussian_resolution_alpha$support,
                                                                   df_gaussian_resolution_alpha$alpha), ]
save(df_gaussian_resolution_alpha, file = rev_output_file(script_dir, "assembled_gaussian_resolution_alpha"))
write.csv(df_gaussian_resolution_alpha,
          file = rev_output_file(script_dir, "assembled_gaussian_resolution_alpha", ext = ".csv"),
          row.names = FALSE)
