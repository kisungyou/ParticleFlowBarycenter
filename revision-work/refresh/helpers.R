# Numerical definitions used by the published Gaussian-resolution and WASP fits.
# No source-directory or private-machine paths are required.
require_transport <- function() {
  if (!requireNamespace("T4transport", quietly=TRUE))
    stop("Install T4transport to score saved supports or compute fits.",call.=FALSE)
}
atomcov <- function(z) crossprod(sweep(z,2,colMeans(z)))/nrow(z)
bures <- function(A,B) {
  A <- (A+t(A))/2; B <- (B+t(B))/2
  ea <- eigen(A,symmetric=TRUE)
  As <- ea$vectors %*% (sqrt(pmax(ea$values,0))*t(ea$vectors))
  M <- As %*% B %*% As
  vals <- eigen((M+t(M))/2,symmetric=TRUE,only.values=TRUE)$values
  sqrt(max(0,sum(diag(A))+sum(diag(B))-2*sum(sqrt(pmax(vals,0)))))
}
diag_ot <- function(z,ys) {
  require_transport()
  m <- nrow(z); K <- length(ys); M <- matrix(0,m,ncol(z)); F <- 0; feas <- 0
  for (y in ys) {
    ot <- T4transport::wasserstein(z,y,p=2)
    F <- F+ot$distance^2/K; M <- M+(ot$plan%*%y)*m/K
    feas <- max(feas,max(abs(rowSums(ot$plan)-1/m)),max(abs(colSums(ot$plan)-1/nrow(y))))
  }
  list(F=F,R=mean(rowSums((z-M)^2)),M=M,feas=feas)
}
native_solver <- function() {
  require_transport()
  fn <- get0("cpp_free_bary_gradient_damped_init",envir=asNamespace("T4transport"),inherits=FALSE)
  if (is.null(fn)) stop(paste(
    "A full rerun requires a T4transport installation exposing",
    "cpp_free_bary_gradient_damped_init (present in the recorded 0.1.9 build).",
    "The version number alone does not establish this capability.",
    "Saved-fit verification, summary rebuilding, and plotting do not require it."),call.=FALSE)
  fn
}
solve_checked <- function(ys,m,alpha,initial_support) {
  fn <- native_solver(); z <- initial_support
  stopifnot(is.matrix(z),nrow(z)==m,ncol(z)==ncol(ys[[1]]))
  s2 <- max(mean(vapply(ys,function(y) sum(diag(atomcov(y))),numeric(1))),1e-12)
  mar <- lapply(ys,function(y) rep(1/nrow(y),nrow(y))); w <- rep(1/length(ys),length(ys))
  nt <- 0; hist <- numeric(); status <- "iteration_limit"; start <- proc.time()[3]
  repeat {
    fit <- fn(ys,mar,w,min(50L,500L-nt),1e-12,alpha,z)
    z <- as.matrix(fit$support); nt <- nt+as.integer(fit$niter)
    hist <- c(hist,as.vector(fit$cost_history))
    d <- diag_ot(z,ys); zz <- (1-alpha)*z+alpha*d$M; dn <- diag_ot(zz,ys)
    delta <- abs(dn$F-d$F)/max(s2,d$F)
    if (delta<=1e-8 && sqrt(d$R/s2)<=1e-4) {status <- "converged";break}
    if (nt>=500L) break
    z <- zz; nt <- nt+1L
  }
  list(support=z,runtime=unname(proc.time()[3]-start),niter=nt,F=d$F,R=d$R,
       residual=sqrt(d$R/s2),delta=delta,status=status,feas=d$feas,history=hist,s2=s2)
}
posterior_metrics <- function(z,inputs) {
  require_transport()
  q <- apply(z,2,quantile,c(.025,.975)); ep <- inputs$test_covariates%*%t(z)
  lo <- apply(ep,1,quantile,.025); hi <- apply(ep,1,quantile,.975)
  eta <- inputs$true_conditional_mean
  c(semi_w2=T4transport::wasserstein(z,inputs$truth,p=2)$distance,
    mean_error=sqrt(sum((colMeans(z)-inputs$posterior_mean)^2)),
    cov_error=bures(atomcov(z),inputs$posterior_cov),
    mean_rmse=sqrt(mean((rowMeans(ep)-eta)^2)),mean_inclusion=mean(eta>=lo & eta<=hi),
    coef_inclusion=mean(c(1,1)>=q[1,] & c(1,1)<=q[2,]),coef_width=mean(q[2,]-q[1,]))
}
