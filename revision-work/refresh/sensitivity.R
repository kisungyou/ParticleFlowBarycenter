#!/usr/bin/env Rscript
# Nested input-size sensitivity from the bundled selected 200-draw inputs.
script <- sub("^--file=","",commandArgs()[grepl("^--file=",commandArgs())][1])
script <- gsub("~+~"," ",script,fixed=TRUE) # Rscript encodes spaces in --file.
base <- dirname(normalizePath(script));source(file.path(base,"helpers.R"))
args <- commandArgs(trailingOnly=TRUE);mode <- if(length(args))args[1] else "--check"
stopifnot(mode %in% c("--check","--assemble","--run","--recompute"))
input <- readRDS(file.path(base,"sensitivity_inputs.rds"));require_transport()
rows <- list();i <- 0L
for(rep in 1:3) {
  case <- input$cases[[as.character(rep)]];fits <- list()
  for(n in c(25,50,100,200)) {
    fn <- file.path(base,sprintf("sensitivity_%d_%d.rds",rep,n))
    if((!file.exists(fn) && mode=="--run") || mode=="--recompute") {
      ys <- lapply(seq_along(case$measures),function(j)case$measures[[j]][case$permutations[[j]][seq_len(n)],,drop=FALSE])
      fit <- solve_checked(ys,100,1,case$initial_support);saveRDS(fit,fn)
    }
    if(!file.exists(fn))stop("Missing saved fit ",basename(fn),"; use --run with the native solver available.")
    fits[[as.character(n)]] <- readRDS(fn)
  }
  for(n in c(25,50,100,200)) {
    fit <- fits[[as.character(n)]];i <- i+1L
    rows[[i]] <- data.frame(rep=rep,draws=n,
      oracle_w2=T4transport::wasserstein(fit$support,input$truth,p=2)$distance,
      reference_w2=T4transport::wasserstein(fit$support,fits[["200"]]$support,p=2)$distance,
      full_input_objective=diag_ot(fit$support,case$measures)$F,
      runtime=fit$runtime,niter=fit$niter,residual=fit$residual,delta=fit$delta,status=fit$status)
  }
}
result <- do.call(rbind,rows)
if(mode=="--check") {
  stopifnot(isTRUE(all.equal(result,read.csv(file.path(base,"draw_sensitivity.csv")),tolerance=1e-10)),all(result$status=="converged"))
  cat("Verified all 12 sensitivity fits by rescoring their saved supports against bundled inputs.\n")
} else {
  write.csv(result,file.path(base,"draw_sensitivity.csv"),row.names=FALSE)
  cat("Rebuilt retained-draw sensitivity summaries.\n")
}
