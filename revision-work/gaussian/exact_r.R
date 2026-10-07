args <- commandArgs(trailingOnly=TRUE)
base <- args[1]
maxiter <- 300L
objective_tol <- 1e-8
residual_tol <- 1e-7
cases <- read.csv(file.path(base, "cases.csv"), stringsAsFactors=FALSE)
stopifnot(requireNamespace("T4transport",quietly=TRUE))
writeLines(c(paste("R",getRversion()),paste("T4transport",packageVersion("T4transport"))),file.path(base,"rerun_R_versions.txt"))
sqdist <- function(X,Y) pmax(outer(rowSums(X*X),rowSums(Y*Y),"+") - 2*tcrossprod(X,Y),0)
for (j in seq_len(nrow(cases))) {
  case <- cases[j,]
  d <- file.path(base,"inputs",case$id)
  out <- file.path(base,"runs",paste0(case$id,"_r_exact"))
  if (file.exists(paste0(out,".csv"))) next
  X <- as.matrix(read.csv(file.path(d,"init.csv"),header=FALSE))
  measures <- lapply(seq_len(4),function(i) as.matrix(read.csv(file.path(d,paste0("measure",i,".csv")),header=FALSE)))
  b <- rep(1/nrow(X),nrow(X)); a <- lapply(measures,function(Y) rep(1/nrow(Y),nrow(Y)))
  previous <- NA_real_; scale2 <- case$scale2; hist <- list()
  start <- proc.time()[3]
  for (it in 0:maxiter) {
    M <- X*0; obj <- 0; marginal <- 0
    for (i in seq_len(4)) {
      C <- sqdist(X,measures[[i]])
      P <- T4transport::wasserstein(X,measures[[i]],p=2,wx=b,wy=a[[i]])$plan
      marginal <- max(marginal,max(abs(rowSums(P)-b)),max(abs(colSums(P)-a[[i]])))
      obj <- obj + sum(P*C)/4
      M <- M + (P %*% measures[[i]])/b/4
    }
    residual <- sum(b*rowSums((M-X)^2))/scale2
    change <- if (is.na(previous)) Inf else abs(previous-obj)/max(1,abs(previous),abs(obj))
    hist[[it+1]] <- data.frame(iter=it,objective=obj,residual=residual,relative_change=change,marginal_error=marginal)
    if (change<=objective_tol && residual<=residual_tol && marginal<=1e-8) {status<-"converged"; break}
    if (it==maxiter) {status<-"iteration_limit"; break}
    previous <- obj; X <- M
  }
  elapsed <- unname(proc.time()[3]-start)
  write.table(X,paste0(out,"_support.csv"),sep=",",row.names=FALSE,col.names=FALSE)
  write.csv(do.call(rbind,hist),paste0(out,"_history.csv"),row.names=FALSE)
  hh <- do.call(rbind,hist)
  write.csv(data.frame(id=case$id,rep=case$rep,support=case$support,method="r_exact",reg_ratio=0,
    runtime_sec=elapsed,niter=it,objective=obj,residual=residual,relative_change=change,
    marginal_error=marginal,status=status,scale2=scale2,max_marginal_error=max(hh$marginal_error),
    max_objective_increase=max(diff(hh$objective))),paste0(out,".csv"),row.names=FALSE)
  cat(case$id,"R exact",it,status,"sec",elapsed,"\n"); flush.console()
}
