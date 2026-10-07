args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=",commandArgs())][1])
script <- gsub("~+~", " ", script, fixed=TRUE) # Decode spaces encoded by Rscript.
here <- dirname(normalizePath(script))
out <- if (length(args)) args[1] else file.path(here,"data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
data(digits, package = "T4transport")
idx <- which(as.integer(digits$label) == 8L)[1L]
img <- as.matrix(digits$image[[idx]])
# Same 256-bin Otsu rule as the submitted visualization script.
xx <- as.vector(img)
h <- hist(xx, breaks = seq(min(xx), max(xx), length.out = 257L), plot = FALSE)
probs <- h$counts / sum(h$counts)
threshold <- h$breaks[2L]
best <- -Inf
for (i in seq_len(255L)) {
  w0 <- sum(probs[1:i])
  w1 <- sum(probs[(i+1L):256L])
  if (w0 == 0 || w1 == 0) next
  mu0 <- sum(h$mids[1:i] * probs[1:i]) / w0
  mu1 <- sum(h$mids[(i+1L):256L] * probs[(i+1L):256L]) / w1
  between <- w0 * w1 * (mu0-mu1)^2
  if (between > best) {
    best <- between
    threshold <- h$breaks[i+1L]
  }
}
write.table(img, file.path(out, "digit8_image.csv"), sep=",", row.names=FALSE, col.names=FALSE)
writeLines(format(threshold, digits=17), file.path(out, "digit8_threshold.txt"))
message("Exported original first digit-8 image; threshold = ", threshold)
