#!/usr/bin/env Rscript
# Reproduce published numerical summaries from bundled inputs and saved fits.
all_args <- commandArgs()
script <- sub("^--file=","",all_args[grepl("^--file=",all_args)][1])
script <- gsub("~+~"," ",script,fixed=TRUE) # Rscript encodes spaces in --file.
base <- dirname(normalizePath(script))
source(file.path(base,"helpers.R"))
args <- commandArgs(trailingOnly=TRUE)
mode <- if (length(args)) args[1] else "--check"
stopifnot(mode %in% c("--check","--assemble","--run","--recompute","--check-native"))
gi <- readRDS(file.path(base,"gaussian_inputs.rds"))
wi <- readRDS(file.path(base,"wasp_inputs.rds"))
initial <- readRDS(file.path(base,"initial_supports.rds"))
gaussian_row <- function(fit,rep,m,alpha) {
  input <- gi[[as.character(rep)]]
  data.frame(rep=rep,support=m,alpha=alpha,
    semi_w2=T4transport::wasserstein(fit$support,input$truth,p=2)$distance,
    mean_error=sqrt(sum((colMeans(fit$support)-input$oracle$mean)^2)),
    cov_error=bures(atomcov(fit$support),input$oracle$var),runtime=fit$runtime,
    niter=fit$niter,F=fit$F,residual=fit$residual,delta=fit$delta,status=fit$status,feas=fit$feas)
}
wasp_row <- function(fit,rep,K,m,alpha) {
  cbind(data.frame(rep=rep,nsplit=K,support=m,alpha=alpha,runtime=fit$runtime,
    niter=fit$niter,F=fit$F,residual=fit$residual,delta=fit$delta,status=fit$status,feas=fit$feas),
    as.data.frame(as.list(posterior_metrics(fit$support,wi))))
}
cached_rows <- function(study) {
  pattern <- if(study=="gaussian") "^gaussian_[0-9]+_[0-9]+_(0.5|1)\\.rds$" else "^wasp_[0-9]+_[0-9]+_[0-9]+_(0.5|1)\\.rds$"
  rows <- lapply(list.files(base,pattern=pattern,full.names=TRUE),function(f) readRDS(f)$row)
  do.call(rbind,rows)
}
if (mode=="--check" || mode=="--check-native") {
  require_transport()
  published_g <- read.csv(file.path(base,"gaussian_resolution.csv"))
  published_w <- read.csv(file.path(base,"wasp_metrics.csv"))
  for (study in c("gaussian","wasp")) {
    raw <- cached_rows(study); ref <- if(study=="gaussian")published_g else published_w
    keys <- if(study=="gaussian")c("rep","support","alpha") else c("rep","nsplit","support","alpha")
    raw <- raw[do.call(order,raw[keys]),];ref <- ref[do.call(order,ref[keys]),]
    rownames(raw)<-rownames(ref)<-NULL
    stopifnot(isTRUE(all.equal(raw,ref,tolerance=1e-12)))
  }
  stopifnot(nrow(published_g)==36,nrow(published_w)==54,
            all(published_g$status=="converged"),all(published_w$status=="converged"))
  g <- readRDS(file.path(base,"gaussian_1_50_0.5.rds"))
  w <- readRDS(file.path(base,"wasp_1_5_100_1.rds"))
  stopifnot(isTRUE(all.equal(gaussian_row(g$fit,1,50,.5),g$row,tolerance=1e-10)),
            isTRUE(all.equal(wasp_row(w$fit,1,5,100,1),w$row,tolerance=1e-10)))
  reference <- as.data.frame(as.list(posterior_metrics(wi$full_mcmc,wi)))
  stopifnot(isTRUE(all.equal(reference,read.csv(file.path(base,"wasp_full_reference.csv")),tolerance=1e-10)))
  cat("Verified all 36 Gaussian and 54 WASP records, representative accuracy metrics, and full-data reference.\n")
  if(mode=="--check-native") {
    saved <- readRDS(file.path(base,"gaussian_1_10_1.rds"))$fit
    reproduced <- solve_checked(gi[["1"]]$dat$measures,10,1,initial[["gaussian_1_10"]])
    discrepancy <- max(abs(saved$support-reproduced$support))
    stopifnot(discrepancy<1e-10,saved$niter==reproduced$niter)
    cat("Native small-case rerun: support max difference",discrepancy,"; iterations",reproduced$niter,"\n")
  }
  quit(save="no",status=0)
}
if(mode %in% c("--run","--recompute")) {
  require_transport(); native_solver()
  for(rep in 1:3)for(m in c(10,25,50,100,200,500))for(alpha in c(1,.5)) {
    fn <- file.path(base,sprintf("gaussian_%d_%d_%s.rds",rep,m,alpha))
    if(file.exists(fn) && mode!="--recompute")next
    input <- gi[[as.character(rep)]]
    fit <- solve_checked(input$dat$measures,m,alpha,initial[[sprintf("gaussian_%d_%d",rep,m)]])
    saveRDS(list(row=gaussian_row(fit,rep,m,alpha),fit=fit,dat=input$dat,oracle=input$oracle),fn)
    cat(basename(fn),fit$status,"\n")
  }
  for(rep in 1:3)for(K in c(5,10,20))for(m in c(50,100,200))for(alpha in c(1,.5)) {
    fn <- file.path(base,sprintf("wasp_%d_%d_%d_%s.rds",rep,K,m,alpha))
    if(file.exists(fn) && mode!="--recompute")next
    fit <- solve_checked(wi$subsets[[sprintf("%d_%d",rep,K)]],m,alpha,initial[[sprintf("wasp_%d_%d_%d",rep,K,m)]])
    saveRDS(list(row=wasp_row(fit,rep,K,m,alpha),fit=fit),fn)
    cat(basename(fn),fit$status,"\n")
  }
}
g <- cached_rows("gaussian");g<-g[order(g$rep,g$support,-g$alpha),]
w <- cached_rows("wasp");w<-w[order(w$rep,w$nsplit,w$support,-w$alpha),]
stopifnot(nrow(g)==36,nrow(w)==54)
write.csv(g,file.path(base,"gaussian_resolution.csv"),row.names=FALSE)
write.csv(w,file.path(base,"wasp_metrics.csv"),row.names=FALSE)
cat("Rebuilt Gaussian and WASP summaries from saved fit records.\n")
