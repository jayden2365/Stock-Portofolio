portfolio_fun <- function(data) {
  
  # Robust GMVP
  portfolioGMVPRobustSphereX <- function(X_hat, delta) {
    N <- ncol(X_hat)
    w <- Variable(N)
    X_ <- scale(X_hat, center = TRUE, scale = FALSE)  # demean
    prob <- Problem(Minimize(norm2(X_ %*% w) + delta*norm2(w)),
                    constraints = list(w >= 0, sum(w) == 1))
    result <- solve(prob)
    return(as.vector(result$getValue(w)))
  }
  
  ## Run Code and paramaters
  X <- as.matrix(diff(log(data$adjusted))[-1])
  
  # Robust GMVP
  delta <- 0.1
  w <- portfolioGMVPRobustSphereX(X, delta)
  return(w)
}
