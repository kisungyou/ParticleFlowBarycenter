# Locate inputs/helpers relative to this script, from any working directory.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_arg <- gsub("~+~", " ", script_arg, fixed=TRUE) # Decode spaces encoded by Rscript.
stopifnot(length(script_arg) == 1L)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
repo_root <- normalizePath(file.path(script_dir, "../.."))
work_dir <- Sys.getenv("MNIST_WORK_DIR", unset = script_dir)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)
source(file.path(repo_root, "code-real-digits/common/rev_utils.R"))
options(stringsAsFactors = FALSE)
stopifnot(requireNamespace("T4transport", quietly = TRUE),
          requireNamespace("jsonlite", quietly = TRUE))
if (file.exists("point_clouds.csv.gz") && file.exists("splits.json") &&
    !("--force" %in% commandArgs(trailingOnly = TRUE))) {
  cat("Prepared inputs already exist. Use MNIST_WORK_DIR for a fresh run or --force to regenerate.\n")
  quit(save = "no")
}
data(digits, package = "T4transport")
labels <- as.integer(digits$label)
make_digit_split <- function(labels, rep_id, n_train_per_digit = 100L,
                             n_test_per_digit = 250L) {
  set.seed(800000L + rep_id)
  labs <- sort(unique(labels)); train_idx <- test_idx <- integer(0)
  for (lab in labs) {
    pool <- which(labels == lab)
    pool <- sample(pool, length(pool), replace = FALSE)
    ntr <- min(n_train_per_digit, length(pool))
    train_lab <- pool[seq_len(ntr)]
    rem <- setdiff(pool, train_lab)
    nts <- min(n_test_per_digit, length(rem))
    train_idx <- c(train_idx, train_lab)
    test_idx <- c(test_idx, rem[seq_len(nts)])
  }
  list(train_idx = sample(train_idx, length(train_idx), replace = FALSE),
       test_idx = sample(test_idx, length(test_idx), replace = FALSE))
}
t0 <- proc.time()[["elapsed"]]
clouds <- lapply(digits$image, rev_img2coords, num_bins = 256L)
preprocessing_sec <- proc.time()[["elapsed"]] - t0
coords <- do.call(rbind, lapply(seq_along(clouds), function(i)
  data.frame(image_id = i, x = clouds[[i]][,1], y = clouds[[i]][,2])))
con <- gzfile("point_clouds.csv.gz", "wt")
write.csv(coords, con, row.names = FALSE)
close(con)
write.csv(data.frame(image_id = seq_along(labels), label = labels,
                     atoms = vapply(clouds, nrow, integer(1))),
          "images.csv", row.names = FALSE)
splits <- lapply(1:3, function(rep_id) {
  sp <- make_digit_split(labels, rep_id)
  # The shuffled training order is fixed by the published split seed.
  # Selection is restricted to its first ten members in each class.
  sp$candidates <- lapply(0:9, function(lab)
    head(sp$train_idx[labels[sp$train_idx] == lab], 10L))
  names(sp$candidates) <- as.character(0:9)
  sp$rep <- rep_id; sp$seed <- 800000L + rep_id
  sp
})
jsonlite::write_json(splits, "splits.json", pretty = TRUE, auto_unbox = TRUE)
# Independently validate Python's exact W2 implementation against the solver
# used by the manuscript on twenty pairs spanning all classes.
pair_ids <- do.call(rbind, lapply(0:9, function(lab) {
  ix <- which(labels == lab)
  rbind(c(ix[1], ix[2]), c(ix[1], which(labels == ((lab + 1) %% 10))[3]))
}))
reference <- data.frame(i = pair_ids[,1], j = pair_ids[,2],
  W2_squared = apply(pair_ids, 1, function(ij)
    T4transport::wasserstein(clouds[[ij[1]]], clouds[[ij[2]]], p = 2)$distance^2))
write.csv(reference, "r_reference_distances.csv", row.names = FALSE)
metadata <- list(R_version = R.version.string,
  T4transport_version = as.character(packageVersion("T4transport")),
  platform = R.version$platform, preprocessing_all_5000_sec = preprocessing_sec,
  data = "T4transport::digits: 5000 images, 500 per digit",
  preprocessing = "Same rev_img2coords as published: Otsu, 256 bins; foreground coordinates in [-1,1]^2; uniform masses",
  criterion = "Mean squared exact W2 to all100 class-training images, minimized among first10 in shuffled training order",
  old_helper_md5 = unname(tools::md5sum(file.path(repo_root, "code-real-digits/common/rev_utils.R"))))
jsonlite::write_json(metadata, "r_metadata.json", pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), "r_session_info.txt")
cat("Prepared", nrow(coords), "foreground points; preprocessing", preprocessing_sec, "seconds\n")
