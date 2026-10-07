# Exact-OT MM updates and evaluation used by the controlled DVQ check.
exact_state <- function(Z, summaries) {
  k <- nrow(Z); S <- length(summaries)
  plans <- lapply(summaries, function(A) T4transport::wasserstein(Z, A, p = 2))
  M <- Reduce(`+`, Map(function(p, A) k * p$plan %*% A, plans, summaries)) / S
  objective <- mean(vapply(plans, function(p) p$distance^2, numeric(1)))
  residual2 <- mean(rowSums((Z - M)^2))
  feasible_error <- max(vapply(plans, function(p) max(abs(rowSums(p$plan) - 1 / k),
                                                       abs(colSums(p$plan) - 1 / k)), numeric(1)))
  list(objective = objective, residual2 = residual2, M = M, feasible_error = feasible_error)
}

exact_dvq <- function(summaries, initial, maxiter = 100L, objective_tol = 1e-8, residual_tol = 1e-6) {
  Z <- as.matrix(initial)
  A <- do.call(rbind, summaries)
  scale2 <- mean(rowSums(sweep(A, 2, colMeans(A))^2))
  if (scale2 <= 0) scale2 <- 1
  previous <- NA_real_; rows <- list(); converged <- FALSE
  for (iteration in 0:maxiter) {
    st <- exact_state(Z, summaries)
    relative_change <- if (is.na(previous)) Inf else abs(st$objective - previous) / max(scale2, abs(previous))
    normalized_residual <- sqrt(st$residual2 / scale2)
    rows[[length(rows) + 1L]] <- data.frame(iteration = iteration, objective = st$objective,
      residual2 = st$residual2, normalized_residual = normalized_residual,
      relative_change = relative_change, feasibility_error = st$feasible_error)
    if (iteration > 0L && relative_change <= objective_tol && normalized_residual <= residual_tol) {
      converged <- TRUE; break
    }
    if (iteration == maxiter) break
    previous <- st$objective; Z <- st$M
  }
  list(centers = Z, initial = initial, history = do.call(rbind, rows),
       converged = converged, iterations = iteration, scale2 = scale2,
       objective_tol = objective_tol, residual_tol = residual_tol, maxiter = maxiter)
}

evaluate <- function(X, Y, centers) {
  d2 <- vapply(seq_len(nrow(centers)), function(j) rowSums(sweep(X, 2, centers[j, ])^2), numeric(nrow(X)))
  pred <- max.col(-d2, ties.method = "first")
  scores <- mclustcomp::mclustcomp(as.integer(Y), pred, types = c("adjrand", "nmi1"))
  list(pred = pred, distortion = mean(d2[cbind(seq_len(nrow(X)), pred)]),
       ARI = scores$scores[1], NMI = scores$scores[2])
}
