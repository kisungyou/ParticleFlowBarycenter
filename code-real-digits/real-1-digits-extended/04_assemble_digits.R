# Assemble digit centroid, classification, and Euclidean baseline results.
script_dir <- getwd()
source(file.path(dirname(script_dir), "common", "rev_utils.R"))

assemble_dir <- function(subdir) {
  files <- sort(list.files(rev_output_dir(script_dir, subdir), pattern = "^result_.*\\.RData$", full.names = TRUE))
  if (length(files) == 0L) return(data.frame())
  rows <- vector("list", length(files))
  for (i in seq_along(files)) {
    load(files[i])
    rows[[i]] <- metrics
  }
  do.call(rbind, rows)
}

df_digit_centroids <- assemble_dir("centroids")
df_digit_classification <- assemble_dir("classification")
df_digit_baselines <- assemble_dir("baselines")

save(df_digit_centroids, df_digit_classification, df_digit_baselines,
     file = rev_output_file(script_dir, "assembled_digits_extended"))
if (nrow(df_digit_centroids) > 0) write.csv(df_digit_centroids, rev_output_file(script_dir, "assembled_digit_centroids", ext = ".csv"), row.names = FALSE)
if (nrow(df_digit_classification) > 0) write.csv(df_digit_classification, rev_output_file(script_dir, "assembled_digit_classification", ext = ".csv"), row.names = FALSE)
if (nrow(df_digit_baselines) > 0) write.csv(df_digit_baselines, rev_output_file(script_dir, "assembled_digit_baselines", ext = ".csv"), row.names = FALSE)
