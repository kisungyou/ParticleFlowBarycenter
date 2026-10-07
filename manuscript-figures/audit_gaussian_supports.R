# Read-only comparison of saved fits; no optimization or new data generation.
suppressPackageStartupMessages(library(T4transport))
this_file <- sub("^--file=", "", commandArgs()[grepl("^--file=",commandArgs())][1])
this_file <- gsub("~+~", " ", this_file, fixed=TRUE) # Decode spaces encoded by Rscript.
here <- dirname(normalizePath(this_file))
out <- file.path(here,"data")
base <- file.path(dirname(here),"revision-work","refresh")
rows <- list(); j <- 0
for (rep in 1:3) for (m in c(10,25,50,100,200,500)) {
  half <- readRDS(file.path(base,sprintf("gaussian_%d_%d_0.5.rds",rep,m)))
  full <- readRDS(file.path(base,sprintf("gaussian_%d_%d_1.rds",rep,m)))
  target_mean <- Reduce("+",lapply(half$dat$measures,colMeans))/length(half$dat$measures)
  j <- j+1
  rows[[j]] <- data.frame(rep=rep,support=m,
    identical_inputs=identical(half$dat,full$dat),
    identical_output=identical(half$fit$support,full$fit$support),
    max_pointwise_difference=max(abs(half$fit$support-full$fit$support)),
    w2_between_outputs=wasserstein(half$fit$support,full$fit$support,p=2)$distance,
    full_mean_deviation=sqrt(sum((colMeans(full$fit$support)-target_mean)^2)),
    half_mean_deviation=sqrt(sum((colMeans(half$fit$support)-target_mean)^2)),
    full_iterations=full$fit$niter,half_iterations=half$fit$niter)
}
write.csv(do.call(rbind,rows),file.path(out,"gaussian_support_audit.csv"),row.names=FALSE)
